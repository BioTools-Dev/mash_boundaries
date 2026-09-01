#!/usr/bin/env python3
"""Arma el directorio de suplementarios con los nombres que usa el manuscrito.

Las tablas del artículo se llaman S1 a S12, pero los ficheros que las contienen
llevan el nombre funcional con que los escribió el paso que los produjo:
`conflicts_s10000_genera.tsv` y no `Table_S1`. Quien llega desde el artículo
busca lo primero y encuentra lo segundo, así que aquí se construye la vista que
el lector espera.

No es una copia mantenida a mano. La correspondencia se lee del propio
manuscrito —de las rutas que cada entrada de la sección de material suplementario
declara— de modo que añadir un fichero a una tabla, o una tabla nueva, se refleja
aquí sin tocar este script y sin que las dos vistas puedan divergir.

Una tabla respaldada por un solo fichero se escribe como un fichero, `Table_S2.tsv`,
porque es lo que un lector quiere descargar. Una respaldada por varios se escribe
como un directorio con esos ficheros dentro, conservando sus nombres de origen:
fundirlos exigiría un esquema común que no tienen, y renombrarlos perdería la
trazabilidad con el paso que los generó.

Una tabla de varios ficheros se entrega además como un `.xlsx` con una hoja por
fichero, que es el formato en que una revista espera un suplementario. El `.tsv`
se conserva al lado: es lo que lee una máquina, y no todo cabe cómodamente en una
hoja de cálculo.

Los ficheros se duplican en lugar de enlazarse. Los enlaces simbólicos viajan mal
entre sistemas y en las descargas comprimidas de las plataformas de código; y como
git guarda el contenido por su huella, dos rutas con el mismo contenido comparten
objeto y la duplicación no pesa en el repositorio.
"""
import argparse, glob, io, os, re, shutil, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xlsx_write


def sources(catalogo):
    """Las rutas que cada entrada del catálogo declara."""
    blk = io.open(catalogo, encoding="utf8").read()
    out = []
    for item in re.split(r"\n- \*\*Table S", blk)[1:]:
        num = item.split(".")[0]
        titulo = item.partition("**")[2].strip().split("\n")[0]
        # Solo las rutas que siguen al marcador de procedencia, y solo mientras
        # sigan encadenadas: una entrada puede nombrar después otros ficheros
        # -por ejemplo para decir que NO forman parte de la tabla- y recogerlos
        # metería en el suplementario justo lo que el texto excluye.
        m = re.search(r"assembled\s+by\s+the\s+pipeline\s+from\s+"
                      r"((?:`[^`]+`(?:,| and|\s)*)+)", item)
        files = []
        for p in re.findall(r"`([^`]+)`", m.group(1) if m else ""):
            if not p.startswith(("results/", "data/")):
                continue
            m = re.search(r"\{([^}]*)\}", p)
            cands = ([p[:m.start()] + o + p[m.end():] for o in m.group(1).split(",")]
                     if m else [p])
            for c in cands:
                files += sorted(glob.glob(c))
        out.append((num, titulo, sorted(set(files))))
    return out


def base_names(files):
    """Nombres de hoja: se quita el prefijo común, que no distingue nada."""
    stems = [os.path.splitext(os.path.basename(f))[0] for f in files]
    pre = os.path.commonprefix(stems) if len(stems) > 1 else ""
    if pre and any(len(s) == len(pre) for s in stems):
        pre = ""                                   # un nombre es prefijo de otro
    taken = set()
    return [xlsx_write.sheet_name(s[len(pre):] or s, taken) for s in stems]


def write_book(path, files):
    """Un .xlsx con una hoja por fichero, leídas en flujo."""
    def rows(f):
        with open(f) as fh:
            for line in fh:
                yield line.rstrip("\n").split("\t")
    xlsx_write.write(path, list(zip(base_names(files), (rows(f) for f in files))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalogue", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir,
        "supplementary", "tables.md"),
        help="el catálogo que declara qué ficheros respaldan cada tabla")
    ap.add_argument("--out", required=True, help="directorio de suplementarios")
    ap.add_argument("--xlsx-max-rows", type=int, default=200000,
                    help="por encima de esto no se genera el .xlsx: una hoja de "
                         "cálculo con cientos de miles de filas es lenta de abrir "
                         "y pesa más que el dato que contiene. 0 lo desactiva")
    a = ap.parse_args()

    # Se reconstruye entero, nunca se parchea, pero solo lo que este paso escribe:
    # el catálogo del que se lee vive en el mismo directorio que describe, y borrarlo
    # dejaría la vista sin su fuente.
    for viejo in glob.glob(os.path.join(a.out, "Table_S*")) + \
                 [os.path.join(a.out, "README.md")]:
        shutil.rmtree(viejo) if os.path.isdir(viejo) else os.path.exists(viejo) and os.remove(viejo)
    os.makedirs(a.out, exist_ok=True)

    index, total = [], 0
    for num, titulo, files in sources(a.catalogue):
        if not files:
            # Un aviso aquí dejaba el directorio incompleto y el paso terminaba
            # bien, que es la peor combinación posible: la vista que el artículo
            # promete salía con una tabla de menos y nada lo señalaba.
            sys.exit("la tabla S%s no resuelve ningún fichero: revisa el catálogo"
                     % num)
        if len(files) == 1:
            dest = os.path.join(a.out, "Table_S%s%s" % (num, os.path.splitext(files[0])[1]))
            shutil.copy2(files[0], dest)
            donde = os.path.basename(dest)
        else:
            d = os.path.join(a.out, "Table_S%s" % num)
            os.makedirs(d)
            for f in files:
                shutil.copy2(f, os.path.join(d, os.path.basename(f)))
            donde = "Table_S%s/ (%d files)" % (num, len(files))
            if a.xlsx_max_rows and sum(1 for f in files for _ in open(f)) <= a.xlsx_max_rows:
                book = os.path.join(a.out, "Table_S%s.xlsx" % num)
                write_book(book, files)
                donde += " and Table_S%s.xlsx" % num
        total += len(files)
        index.append((num, titulo, donde, files))
        print("  S%-3s %2d fichero(s) -> %s" % (num, len(files), donde))

    with io.open(os.path.join(a.out, "README.md"), "w", encoding="utf8") as fh:
        fh.write("# Supplementary tables\n\n"
                 "The supplementary tables of the manuscript, under the names the paper\n"
                 "uses. This directory is assembled from `results/` and `data/` by\n"
                 "`scripts/33_supplementary.py`, which reads the catalogue in `tables.md`,\n"
                 "so the two views cannot disagree. Each entry lists the files it was built\n"
                 "from, which is where the pipeline writes them.\n\n")
        for num, titulo, donde, files in index:
            fh.write("## Table S%s\n\n%s\n\n**Here:** `%s`\n\n**Built from:** %s\n\n"
                     % (num, titulo, donde,
                        ", ".join("`%s`" % f for f in files)))
    print("\n%d ficheros en %d tablas" % (total, len(index)))


if __name__ == "__main__":
    main()
