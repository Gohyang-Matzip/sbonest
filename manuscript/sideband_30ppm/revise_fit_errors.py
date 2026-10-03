"""Add all fit uncertainties as Donghan Lee redlines.

Use the bundled document runtime Python. Preserve the accepted previous version.
"""
from copy import deepcopy
from datetime import datetime, timezone
from difflib import SequenceMatcher
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from zipfile import ZipFile
from lxml import etree as E

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'Sideband_CEST_1p2GHz_parameters_13C.docx'
CLEAN = HERE/'Sideband_CEST_1p2GHz_parameters_errors.docx'
TRACKED = HERE/'Sideband_CEST_1p2GHz_parameters_errors_tracked.docx'
QA = HERE/'qa/fit_errors_20261003'
W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
NS = {'w': W, 'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'}
AUTHOR = 'Donghan Lee'
WHEN = datetime.now(timezone.utc).isoformat(timespec='seconds').replace('+00:00','Z')
assert not CLEAN.exists() and not TRACKED.exists()
QA.mkdir(exist_ok=True)
source_bytes = SOURCE.read_bytes()
with ZipFile(SOURCE) as z:
    parts = {i.filename:(i,z.read(i.filename)) for i in z.infolist()}
root = E.fromstring(parts['word/document.xml'][1])
original = deepcopy(root)
body = root.find('w:body',NS)
paras = body.findall('w:p',NS)
assert not root.xpath('//w:ins|//w:del',namespaces=NS)
next_id = max([int(x) for x in root.xpath('//@w:id',namespaces=NS)]+[0])+1
ledger = []


def tag(s):
    return f'{{{W}}}{s}'


def txt(p):
    return ''.join(p.xpath('.//w:t/text()',namespaces=NS))


def rev(kind):
    global next_id
    x=E.Element(tag(kind),{tag('id'):str(next_id),tag('author'):AUTHOR,tag('date'):WHEN})
    next_id+=1
    return x


def run(value, rpr=None, deleted=False):
    r=E.Element(tag('r'))
    if rpr is not None:
        r.append(deepcopy(rpr))
    for i,line in enumerate(value.split('\n')):
        if i:E.SubElement(r,tag('br'))
        t=E.SubElement(r,tag('delText' if deleted else 't'))
        t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
        t.text=line
    return r


def replace(i,new):
    p=paras[i] if isinstance(i,int) else i; old=txt(p)
    assert not p.xpath('.//w:drawing|.//m:oMath',namespaces=NS)
    rpr=p.find('w:r/w:rPr',NS)
    for child in list(p):
        if child.tag!=tag('pPr'):p.remove(child)
    if isinstance(i,int):
        a,b=re.findall(r'\S+\s*',old),re.findall(r'\S+\s*',new)
        suffix=''
    else:
        # Keep unchanged digits visible once in numeric table redlines.
        point,separator,error=new.partition('\n')
        a,b=list(old),list(point)
        suffix=separator+error
    for op,l,r,s,t in SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
        if op=='equal':p.append(run(''.join(a[l:r]),rpr))
        else:
            if l!=r:
                x=rev('del');x.append(run(''.join(a[l:r]),rpr,True));p.append(x)
            if s!=t:
                x=rev('ins');x.append(run(''.join(b[s:t]),rpr));p.append(x)
    if suffix:
        x=rev('ins');x.append(run(suffix,rpr));p.append(x)
    ledger.append(dict(kind='replace',paragraph=i if isinstance(i,int) else 'table cell',old=old,new=new))


def insert(anchor,value,template=6):
    p=E.Element(tag('p'))
    pp=deepcopy(paras[template].find('w:pPr',NS))
    if pp is None:pp=E.Element(tag('pPr'))
    for x in pp.findall('w:pageBreakBefore',NS):pp.remove(x)
    rp=pp.find('w:rPr',NS)
    if rp is None:rp=E.SubElement(pp,tag('rPr'))
    rp.append(rev('ins'));p.append(pp)
    x=rev('ins');x.append(run(value,paras[template].find('w:r/w:rPr',NS)));p.append(x)
    anchor.addnext(p)
    ledger.append(dict(kind='insert',new=value))
    return p


errors=json.loads((HERE/'results/fit_errors_20261003/summary.json').read_text())
fits=errors['fits']
tables=body.findall('w:tbl',NS)


def estimate(key,name,multiline=False):
    v=fits[key]['parameters'][name]
    digits=max(0,1-math.floor(math.log10(v['stderr'])))
    return f"{v['value']:.{digits}f}"+('\n± ' if multiline else ' ± ')+f"{v['stderr']:.{digits}f}"


def cells(table):
    return [[cell.find('w:p',NS) for cell in row.findall('w:tc',NS)] for row in table.findall('w:tr',NS)]


for table_index,keys in ((2,('NH_full','NH_masked','NH_fixed')),
                        (4,('NH_full','N_noisy_full','NH_masked','N_noisy_masked'))):
    rows=cells(tables[table_index])
    for row,key in zip(rows[1 if table_index==2 else 2:],keys,strict=True):
        for col,name in ((2,'kex'),(3,'pB_percent'),(4,'v1n_scale')):
            if key=='NH_fixed' and name=='v1n_scale':continue
            replace(row[col],estimate(key,name,True))

param_names=('kab','kba','v1n_scale','A1.peak_ppm','A1.dw_ppm','A1.R1','A1.R2a','A1.R2b','delta_B')
for row,name in zip(cells(tables[5])[1:],param_names,strict=True):
    for cell,key in zip(row[2:],('NH_full','N_noisy_full','NH_masked','N_noisy_masked'),strict=True):
        replace(cell,estimate(key,name,True))

# Keep all repeated numerical estimates in prose consistent with the tables.
replacements={
    3:[('kex=321.4 s⁻¹',f"kex={estimate('N_noisy_full','kex')} s⁻¹"),
       ('300.2 s⁻¹ for',f"{estimate('NH_full','kex')} s⁻¹ for"),
       ('to 7.80 s⁻¹',f"to {estimate('N_noisy_full','A1.R2b')} s⁻¹"),
       ('with 15.22 s⁻¹',f"with {estimate('NH_full','A1.R2b')} s⁻¹")],
    84:[('kex=300.224 s⁻¹',f"kex={estimate('NH_full','kex')} s⁻¹"),
        ('pB=4.9946%',f"pB=({estimate('NH_full','pB_percent')})%"),
        ('scale=1.08138',f"scale={estimate('NH_full','v1n_scale')}"),
        ('300.316 s⁻¹, 4.9896% and 1.08205',f"{estimate('NH_masked','kex')} s⁻¹, ({estimate('NH_masked','pB_percent')})% and {estimate('NH_masked','v1n_scale')}")],
    89:[('5.6654%',f"({estimate('NH_fixed','pB_percent')})%")],
    95:[('to 321.4 s⁻¹',f"to {estimate('N_noisy_full','kex')} s⁻¹"),
        ('to 4.812%',f"to ({estimate('N_noisy_full','pB_percent')})%"),
        ('to 1.0958',f"to {estimate('N_noisy_full','v1n_scale')}"),
        ('to 7.80 s⁻¹',f"to {estimate('N_noisy_full','A1.R2b')} s⁻¹"),
        ('kex=300.2 s⁻¹',f"kex={estimate('NH_full','kex')} s⁻¹"),
        ('pB=4.995%',f"pB=({estimate('NH_full','pB_percent')})%"),
        ('scale=1.0814',f"scale={estimate('NH_full','v1n_scale')}"),
        ('R2N=15.22 s⁻¹',f"R2N={estimate('NH_full','A1.R2b')} s⁻¹")],
    96:[('kex remains 321.2 s⁻¹',f"kex remains {estimate('N_noisy_masked','kex')} s⁻¹"),
        ('pB becomes 4.710%',f"pB becomes ({estimate('N_noisy_masked','pB_percent')})%"),
        ('kex=300.3 s⁻¹',f"kex={estimate('NH_masked','kex')} s⁻¹"),
        ('kex=321.1 s⁻¹',f"kex={estimate('N_exact_full','kex')} s⁻¹"),
        ('and 320.8 s⁻¹',f"and {estimate('N_exact_masked','kex')} s⁻¹")],
    102:[('R2NA=11.971 s⁻¹',f"R2NA={estimate('NH_full','A1.R2a')} s⁻¹"),
         ('R2NB=15.218 s⁻¹',f"R2NB={estimate('NH_full','A1.R2b')} s⁻¹"),
         ('R2NB=7.554 s⁻¹',f"R2NB={estimate('N_exact_full','A1.R2b')} s⁻¹")],
    152:[('kex=321.4 s⁻¹',f"kex={estimate('N_noisy_full','kex')} s⁻¹"),
         ('pB=4.812%',f"pB=({estimate('N_noisy_full','pB_percent')})%"),
         ('300.2 s⁻¹, 4.995%',f"{estimate('NH_full','kex')} s⁻¹, ({estimate('NH_full','pB_percent')})%"),
         ('kex=321.2 s⁻¹',f"kex={estimate('N_noisy_masked','kex')} s⁻¹"),
         ('321.1 s⁻¹(전체)',f"{estimate('N_exact_full','kex')} s⁻¹(전체)"),
         ('320.8 s⁻¹(masked)',f"{estimate('N_exact_masked','kex')} s⁻¹(masked)")],
    153:[('kAB=15.466',f"kAB={estimate('N_noisy_full','kab')}"),
         ('kBA=305.931',f"kBA={estimate('N_noisy_full','kba')}"),
         ('scale=1.09581',f"scale={estimate('N_noisy_full','v1n_scale')}"),
         ('δNA=119.98824',f"δNA={estimate('N_noisy_full','A1.peak_ppm')}"),
         ('ΔδN=3.01667',f"ΔδN={estimate('N_noisy_full','A1.dw_ppm')}"),
         ('R1N=1.5096',f"R1N={estimate('N_noisy_full','A1.R1')}"),
         ('R2NA=11.506',f"R2NA={estimate('N_noisy_full','A1.R2a')}"),
         ('R2NB=7.796',f"R2NB={estimate('N_noisy_full','A1.R2b')}"),
         ('R2NB=15.218',f"R2NB={estimate('NH_full','A1.R2b')}"),
         ('R2NB는 8.465',f"R2NB는 {estimate('N_noisy_masked','A1.R2b')}")],
}
for i,pairs in replacements.items():
    value=txt(paras[i])
    for old,new in pairs:
        assert old in value,(i,old)
        value=value.replace(old,new)
    if i==3:value=value.replace('The rate displacement persists without noise.',
        'Quoted uncertainties are formal local standard errors (1 SE), conditional on each model. The rate displacement persists without noise.')
    replace(i,value)

insert(paras[65],
    'All fitted estimates are reported with one local standard error (estimate ± 1 SE). At each fitted parameter vector, J is the Jacobian of residuals divided by their supplied intensity errors, and C=(JᵀJ)⁻¹ is evaluated by singular-value decomposition. SEᵢ=√Cᵢᵢ; C is not multiplied by χ²/dof. For kex=kAB+kBA, pB=kAB/kex and δB=δA+Δδ, we use SE(g)=√(∇gᵀC∇g), including cross-covariances. Population SE values accompanying percentages are in percentage points. These fitted-point SE values differ from the generating-point information estimates in Table 4. For misspecified nitrogen-only models, they quantify within-model local curvature and exclude model discrepancy; they are not validated confidence intervals. Noiseless checks retain the assigned σ=0.001 weights, so their nonzero SE describes hypothetical measurement noise, not observed replicate scatter. Known generating values and fixed inputs have no fitted SE.',65)
replace(86,txt(paras[86])+' Fitted kex, pB and RF scale are estimate ± 1 SE from the fitted-point covariance without χ²/dof rescaling; pB errors are percentage points. The fixed RF scale has no fitted error.')
replace(97,txt(paras[97]).replace('Values are point estimates; the known generating values are shown for reference.',
    'Fitted values are estimate ± 1 SE from the fitted-point covariance without χ²/dof rescaling; pB SE values are percentage points. N-only SE is formal and excludes model discrepancy. The known generating values are shown without fitted errors.'))
replace(103,txt(paras[103]).replace('These are point estimates from one noise realization.',
    'All fitted entries show the estimate on the first line and ± 1 SE on the second line. SE is from the fitted-point covariance without χ²/dof rescaling; the δNB SE includes Cov(δNA,ΔδN). N-only SE is formal and excludes model discrepancy. These results use one noise realization.'))
replace(113,txt(paras[113])+' Table 7 gives every carbon fit parameter and the derived kex, pB and δCB as estimate ± 1 SE under the assigned σ=0.001 weights; these are hypothetical-noise sensitivities for this noiseless check.')
replace(131,txt(paras[131])+' collect_fit_errors.py reconstructs the full covariance at the unchanged saved optima, verifies all stored parameter SE values and fit objectives, and propagates derived errors. The full covariance matrices and source hashes are in results/fit_errors_20261003/summary.json; revise_fit_errors.py adds the corresponding uncertainties to Tables 3, 5, 6 and 7 and all repeated numerical fit estimates in the text.')
insert(paras[153],
    '모든 fitted parameter는 추정값 ± 1 SE로 표시했다. kex·pB·δB의 SE에는 parameter 간 공분산을 포함했다. pB(%)의 SE 단위는 percentage point이다. χ²/dof에 따른 오차 재조정은 하지 않았다. 일반 fitting의 형식적 SE는 모델 불일치를 포함하지 않는다. ¹³C noiseless 검증의 모든 parameter와 SE는 Table 7에 있으며, 이 SE는 실제 반복 측정 오차가 아니라 가정한 intensity σ=0.001에 대한 민감도이다. 참값·고정값에는 fitted SE를 붙이지 않았다.',153)

# Carbon fit table: all eight free parameters and the three derived quantities.
carbon_specs=[('kAB / s⁻¹','kab',15),('kBA / s⁻¹','kba',285),('ν₁C scale','v1n_scale',1.08),
    ('δCA / ppm','A1.peak_ppm',55),('ΔδC / ppm','A1.dw_ppm',1.2),
    ('R1C / s⁻¹','A1.R1',1.2),('R2CA / s⁻¹','A1.R2a',15),('R2CB / s⁻¹','A1.R2b',20),
    ('kex / s⁻¹ (derived)','kex',300),('pB / % (derived)','pB_percent',5),('δCB / ppm (derived)','delta_B',56.2)]
table=deepcopy(tables[5])
templates=table.findall('w:tr',NS)
for row in templates:table.remove(row)
grid=table.find('w:tblGrid',NS)
for col in list(grid):grid.remove(col)
widths=(2900,2100,4115)
for width in widths:E.SubElement(grid,tag('gridCol'),{tag('w'):str(width)})
carbon_rows=[['Parameter','Generating','¹³C fit ± 1 SE']]+[[label,f'{true:g}',estimate('CH_exact',name)] for label,name,true in carbon_specs]
for i,values in enumerate(carbon_rows):
    row=deepcopy(templates[0 if i==0 else 1])
    rp=row.find('w:trPr',NS)
    if rp is None:rp=E.SubElement(row,tag('trPr'))
    rp.append(rev('ins'))
    tc=row.findall('w:tc',NS)
    for c in tc[3:]:row.remove(c)
    for cell,value,width in zip(tc[:3],values,widths,strict=True):
        cell.find('w:tcPr/w:tcW',NS).set(tag('w'),str(width))
        p=cell.find('w:p',NS);rpr=p.find('.//w:rPr',NS)
        for child in list(p):
            if child.tag!=tag('pPr'):p.remove(child)
        x=rev('ins');x.append(run(value,rpr));p.append(x)
    table.append(row)
paras[113].addnext(table)
insert(table,
    'Table 7. Carbon parameter recovery and local SE for the noiseless isolated ¹³C–¹H implementation check. All eight free parameters and three derived quantities are listed. Errors are 1 SE from the full fitted-point covariance with assigned absolute intensity σ=0.001, without χ²/dof rescaling. They express hypothetical-noise sensitivity, not scatter observed in these noiseless data or experimental uncertainty. Derived errors include parameter cross-covariances; pB SE is in percentage points. Fixed inputs (including JCH=140 Hz, R1H=2 s⁻¹ and R2H=25 s⁻¹) have no fitted SE.',103)


def resolve(doc,accept):
    r=deepcopy(doc)
    if not accept:
        for table in list(r.xpath('//w:tbl',namespaces=NS)):
            rows=table.findall('w:tr',NS)
            if rows and all(x.find('w:trPr/w:ins',NS) is not None for x in rows):table.getparent().remove(table)
        for p in list(r.xpath('//w:p[w:pPr/w:rPr/w:ins]',namespaces=NS)):p.getparent().remove(p)
    for x in reversed(r.xpath('//w:ins|//w:del',namespaces=NS)):
        parent=x.getparent()
        if (x.tag==tag('ins'))==accept:
            position=parent.index(x)
            for child in list(x):
                for t in child.iter(tag('delText')):t.tag=tag('t')
                parent.insert(position,child);position+=1
        parent.remove(x)
    return r


def xml(x):
    return E.tostring(x,xml_declaration=True,encoding='UTF-8',standalone=True)


accepted=resolve(root,True)
rejected=resolve(root,False)
assert [txt(p) for p in original.xpath('//w:p',namespaces=NS)]==[txt(p) for p in rejected.xpath('//w:p',namespaces=NS)]
for query in ('//m:oMath','//w:drawing'):
    for target in (accepted,rejected):
        assert [E.tostring(x) for x in original.xpath(query,namespaces=NS)]==[E.tostring(x) for x in target.xpath(query,namespaces=NS)]
assert len(accepted.xpath('//w:tbl',namespaces=NS))==8
assert not re.search(r'1200\s*MHz',txt(accepted))
assert set(root.xpath('//w:ins/@w:author|//w:del/@w:author',namespaces=NS))=={AUTHOR}
core=E.fromstring(parts['docProps/core.xml'][1])
core.find('{http://purl.org/dc/terms/}modified').text=WHEN
core.find('{http://schemas.openxmlformats.org/package/2006/metadata/core-properties}lastModifiedBy').text=AUTHOR
for path,doc,tracking in ((TRACKED,root,True),(CLEAN,accepted,False)):
    settings=E.fromstring(parts['word/settings.xml'][1])
    if tracking:
        settings.insert(0,E.Element(tag('trackRevisions')))
        settings.insert(1,E.Element(tag('revisionView'),{tag('markup'):'1',tag('insDel'):'1'}))
    replacements={'word/document.xml':xml(doc),'docProps/core.xml':xml(core),'word/settings.xml':xml(settings)}
    with ZipFile(path,'w') as z:
        for name,(info,data) in parts.items():z.writestr(info,replacements.get(name,data))
assert SOURCE.read_bytes()==source_bytes
(QA/'revision_audit.json').write_text(json.dumps(dict(source=str(SOURCE),source_sha256=sha256(source_bytes).hexdigest(),
    outputs={str(p):sha256(p.read_bytes()).hexdigest() for p in (CLEAN,TRACKED)},author=AUTHOR,
    operations=ledger,carbon_table=carbon_rows,source_results='results/fit_errors_20261003/summary.json',
    reject_restores_original_text=True,original_figures_and_equations_preserved=True,visual_qa='pending'),ensure_ascii=False,indent=2)+'\n')
print(CLEAN)
print(TRACKED)
