#!/usr/bin/env bash
# Render the manuscript to .docx, for co-authors and for journal submission.
#
# The master is the Markdown file: it is what versions cleanly, what the figures
# and tables are cross-checked against, and what the audit trail refers to. The
# .docx is a derived product and is regenerated rather than edited, so that a
# correction never has to be made twice.
#
# Two details make the conversion work. Pandoc resolves image paths against the
# working directory rather than against the input file, so the resource path has
# to include the repository root for `../figures/*.png` to resolve. And the input
# is read as pandoc's own Markdown flavour rather than as CommonMark, because the
# abstract carries superscripts in `10^8^` form that CommonMark would render
# literally.
#
# Pandoc then leaves the tables unreadable on a portrait page, and for a reason
# worth recording: it gives every column the same width and sets the cell text at
# body size, so a nine-column table of percentages wraps every heading onto three
# lines while a one-word stratum column keeps a tenth of the page. The step below
# fixes both — cell text goes to TABLE_PT and the table is switched to autofit, so
# each column takes the width its content needs. Neither can be expressed in
# Markdown, and doing it by hand in Word would have to be redone on every
# regeneration.
#
# The figures travel as the 300 dpi PNGs, which is what a .docx can embed. The
# vector masters in `figures/*.svg` and `figures/*.pdf` are what a journal should
# receive alongside; the README says so and the manuscript repeats it.

set -euo pipefail
source "$(dirname "$0")/../config.sh"

# The manuscript by default; any other Markdown of the same directory when named,
# so that the addenda handed to a co-author go through this same conversion — the
# table fitting and the figure embedding are not worth reimplementing twice, and
# a document that reaches a reader by a different path is a document with
# different formatting.
SRC=$(realpath "${1:-$MANUSCRIPT/mash_boundaries_paper.md}" 2>/dev/null || echo "${1:-}")
OUT=${SRC%.md}.docx
TABLE_PT=${TABLE_PT:-7}

[[ -f $SRC ]] || { echo "no está el documento: $SRC"; exit 1; }
command -v pandoc >/dev/null || { echo "hace falta pandoc"; exit 1; }

cd "$MANUSCRIPT"
pandoc "$SRC" \
    --from markdown \
    --to docx \
    --resource-path "$MANUSCRIPT:$WORK" \
    --output "$OUT"

python3 - "$OUT" "$TABLE_PT" <<'PY'
"""Make the tables fit the page, and report what the file ended up carrying."""
import os
import re
import shutil
import sys
import zipfile

path, pt = sys.argv[1], int(sys.argv[2])
half = pt * 2                       # Word measures font size in half-points
size = '<w:sz w:val="%d"/><w:szCs w:val="%d"/>' % (half, half)

src = zipfile.ZipFile(path)
doc = src.read("word/document.xml").decode("utf8")


def shrink(run):
    if "<w:rPr>" in run:            # keep the bold or italic already there
        return run.replace("<w:rPr>", "<w:rPr>" + size, 1)
    return run.replace("<w:r>", "<w:r><w:rPr>" + size + "</w:rPr>", 1)


# Word pads every cell by 0.075 in on each side. Over a ten-column table that
# is a quarter of the printable width spent on whitespace, which is what forces
# a figure like 456,276,736 to break across two lines. The padding goes down to
# 0.03 in, where the columns still read as separate and the numbers stay whole.
MARGIN = ('<w:tblCellMar>'
          '<w:top w:w="0" w:type="dxa"/><w:left w:w="43" w:type="dxa"/>'
          '<w:bottom w:w="0" w:type="dxa"/><w:right w:w="43" w:type="dxa"/>'
          '</w:tblCellMar>')


def fit(table):
    table = table.replace(
        '<w:tblW w:type="pct" w:w="5000" />',
        '<w:tblW w:type="pct" w:w="5000" />'
        '<w:tblLayout w:type="autofit" />' + MARGIN, 1)
    return re.sub(r"<w:r>.*?</w:r>", lambda m: shrink(m.group(0)), table,
                  flags=re.S)


doc, n = re.subn(r"<w:tbl>.*?</w:tbl>", lambda m: fit(m.group(0)), doc,
                 flags=re.S)

# The document has to arrive editable, and that is not automatic. Word refuses to
# edit a file that carries document protection or that opens with revision
# tracking enforced, and neither would be visible in the Markdown. Nothing here
# adds them, but a reference document could, so they are removed rather than
# assumed absent.
GUARDED = ("documentProtection", "writeProtection", "trackChanges",
           "readOnlyRecommended", "revisionView")
settings = "word/settings.xml"
sett = src.read(settings).decode("utf8") if settings in src.namelist() else None
stripped = []
if sett:
    for tag in GUARDED:
        # las dos formas: <w:tag .../> y <w:tag ...></w:tag>
        sett, k = re.subn(r"<w:%s\b[^>]*?(?:/>|>.*?</w:%s>)" % (tag, tag), "",
                          sett, flags=re.S)
        if k:
            stripped.append("%s x%d" % (tag, k))

# ZIP entries need an explicit permission too. Left at zero, which is what pandoc
# writes, zipfile substitutes 0o600 — owner-only — so every part of the package
# ends up unreadable to anyone else the moment it is unpacked or copied through a
# share. 0o644 plus the DOS archive bit is what an ordinary Word file carries.
ATTR = (0o644 << 16) | 0x20

tmp = path + ".tmp"
with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
    for item in src.infolist():
        if item.filename == "word/document.xml":
            data = doc.encode("utf8")
        elif sett is not None and item.filename == settings:
            data = sett.encode("utf8")
        else:
            data = src.read(item.filename)
        item.external_attr = ATTR
        out.writestr(item, data)
src.close()
shutil.move(tmp, path)
os.chmod(path, 0o644)

print("escrito %s" % path)
z = zipfile.ZipFile(path)
d = z.read("word/document.xml").decode("utf8")
print("  figuras incrustadas: %d" % sum(1 for f in z.namelist()
                                        if f.startswith("word/media/")))
print("  tablas ajustadas:    %d de %d, a %d pt"
      % (n, d.count("<w:tbl>"), pt))
print("  superíndices:        %d" % d.count('vertAlign w:val="superscript"'))

# se comprueba, no se da por hecho
s2 = z.read(settings).decode("utf8") if settings in z.namelist() else ""
bad = [t for t in GUARDED if t in s2]
marks = sum(d.count("<w:%s " % t) for t in ("ins", "del", "rPrChange", "pPrChange"))
attrs = {i.external_attr for i in z.infolist()}
print("  editable:            %s"
      % ("sí" if not bad and not marks and attrs == {ATTR} else
         "NO — %s" % (bad or "%d marcas de revisión" % marks or attrs)))
if stripped:
    print("  se retiró de settings.xml: %s" % ", ".join(stripped))
PY
