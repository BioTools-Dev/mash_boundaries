"""Lector mínimo de .xlsx, sin dependencias externas."""
import re, zipfile
import xml.etree.ElementTree as ET
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

def read(path, sheet=1, limit=None):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")):
            shared.append("".join(t.text or "" for t in si.iter(NS + "t")))
    name = "xl/worksheets/sheet%d.xml" % sheet
    rows = []
    for r in ET.fromstring(z.read(name)).iter(NS + "row"):
        cells = {}
        for c in r.iter(NS + "c"):
            ref = c.get("r") or ""
            col = re.match(r"[A-Z]+", ref).group(0) if ref else str(len(cells))
            v = c.find(NS + "v")
            if c.get("t") == "s" and v is not None:
                val = shared[int(v.text)]
            elif c.get("t") == "inlineStr":
                val = "".join(t.text or "" for t in c.iter(NS + "t"))
            else:
                val = v.text if v is not None else ""
            cells[col] = val
        order = sorted(cells, key=lambda k: (len(k), k))
        rows.append([cells[k] for k in order])
        if limit and len(rows) >= limit:
            break
    return rows
