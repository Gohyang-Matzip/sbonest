"""MkDocs hook: links that leave docs/ (sources, results, README) point at the repository."""
BLOB = "https://github.com/Gohyang-Matzip/sbonest/blob/main/"


def on_page_markdown(markdown, page, **_):
    """Rewrite relative links that climb out of the docs folder to GitHub URLs."""
    outside = "](" + "../" * (page.file.src_uri.count("/") + 1)
    return markdown.replace(outside, "](" + BLOB)
