# SBONEST

Exact NH sideband CEST fitting with checkpoints, uncertainty analysis and design
tools. The code, examples and archived evidence are in the
[repository](https://github.com/Gohyang-Matzip/sbonest).

| Start here | |
|---|---|
| [First fit, step by step](guides/DUMMY_GUIDE.md) ([한국어](guides/DUMMY_GUIDE.ko.md)) | Install, run the bundled demo, read the outputs |
| [SBONEST manual](manual/SBONEST_MANUAL.md) ([한국어](manual/SBONEST_MANUAL.ko.md)) | Inputs, configuration, analyses, command line, web runner |
| [API reference](API_REFERENCE.md) | Public functions of every module (generated) |
| [Validation](reports/VALIDATION.md) | What has been checked, on synthetic data, and how |
| [Inherited ONEST manual](manual/ONEST_MANUAL.md) ([한국어](manual/ONEST_MANUAL.ko.md)) | The original ONEST models and workflows |

Install the command line with `pip install .` from a checkout, or use the
container image:

```bash
docker run --rm -v "$PWD:/work" ghcr.io/gohyang-matzip/sbonest init-demo --out demo
docker run --rm -v "$PWD:/work" ghcr.io/gohyang-matzip/sbonest fit demo/fit.json --no-pdf
```

All uncertainty statements are within-model statements under the supplied
absolute sigma; the validation evidence is synthetic unless a report says otherwise.
