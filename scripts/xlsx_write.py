"""Escritor mínimo de .xlsx, sin dependencias externas.

El proyecto no exige más biblioteca que matplotlib, y añadir una para escribir
una hoja de cálculo sería añadirla al lector que quiera reproducir. Un .xlsx es
un zip de XML, así que se escribe directamente.

Se usan cadenas en línea en lugar de la tabla de cadenas compartidas: cuesta algo
de tamaño y ahorra una pasada previa sobre datos que aquí llegan a un millón de
filas. Las hojas se escriben en flujo, sin construir el XML en memoria.

Cada entrada del zip se declara con fecha fija y permisos explícitos, y ninguna
de las dos cosas es cosmética. Con la fecha por defecto —la de generación— dos
ejecuciones sobre los mismos datos producen ficheros distintos byte a byte, de
modo que regenerar el directorio de suplementarios ensucia el repositorio con
megabytes de diferencia binaria que no corresponden a ningún cambio; y un
suplementario que no se puede reconstruir idéntico no es reproducible. Con los
permisos a cero, que es lo que `zipfile` escribe si no se le dice otra cosa, al
desempaquetar sustituye `0o600` y las piezas salen ilegibles para cualquiera que
no sea el dueño — el mismo defecto que hubo que corregir en el paso 25.
"""
import re, zipfile

INVALID = re.compile(r"[\[\]:*?/\\]")
NUM = re.compile(r"^-?\d+(\.\d+)?([eE][-+]?\d+)?$")
FECHA = (1980, 1, 1, 0, 0, 0)                    # el origen del formato zip
ATTR = (0o644 << 16) | 0x20                      # rw-r--r-- más el bit de archivo


def _entrada(nombre):
    zi = zipfile.ZipInfo(nombre, FECHA)
    zi.compress_type = zipfile.ZIP_DEFLATED
    zi.external_attr = ATTR
    return zi


def sheet_name(raw, taken):
    """Excel: 31 caracteres, sin []:*?/\\ y sin repetir."""
    n = INVALID.sub("_", raw)[:31] or "sheet"
    if n not in taken:
        taken.add(n)
        return n
    for i in range(2, 1000):                     # desambiguar sin pasarse de 31
        cand = "%s_%d" % (n[:31 - len(str(i)) - 1], i)
        if cand not in taken:
            taken.add(cand)
            return cand
    raise ValueError("no hay nombre libre para %r" % raw)


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;"))


def _col(i):
    s = ""
    while True:
        s = chr(ord("A") + i % 26) + s
        i = i // 26 - 1
        if i < 0:
            return s


def write(path, sheets):
    """sheets: lista de (nombre, iterable de filas, cada fila lista de cadenas)."""
    cols = [_col(i) for i in range(1024)]
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(_entrada("[Content_Types].xml"),
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            + "".join('<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % (i + 1)
                      for i in range(len(sheets))) + '</Types>')
        z.writestr(_entrada("_rels/.rels"),
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            '</Relationships>')
        z.writestr(_entrada("xl/workbook.xml"),
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
            + "".join('<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (_esc(n), i + 1, i + 1)
                      for i, (n, _) in enumerate(sheets)) + '</sheets></workbook>')
        z.writestr(_entrada("xl/_rels/workbook.xml.rels"),
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + "".join('<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>' % (i + 1, i + 1)
                      for i in range(len(sheets))) + '</Relationships>')
        for i, (name, rows) in enumerate(sheets):
            with z.open(_entrada("xl/worksheets/sheet%d.xml" % (i + 1)), "w") as fh:
                fh.write(b'<?xml version="1.0" encoding="UTF-8"?>'
                         b'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>')
                for r, row in enumerate(rows, 1):
                    buf = ['<row r="%d">' % r]
                    for c, v in enumerate(row):
                        ref = "%s%d" % (cols[c], r)
                        if v == "":
                            continue                       # celda vacía, no se escribe
                        if NUM.match(v):
                            buf.append('<c r="%s"><v>%s</v></c>' % (ref, v))
                        else:
                            buf.append('<c r="%s" t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>'
                                       % (ref, _esc(v)))
                    buf.append("</row>")
                    fh.write("".join(buf).encode("utf8"))
                fh.write(b"</sheetData></worksheet>")
