#!/usr/bin/env python3
"""Los suplementarios de Bi et al. (2026), descargados y puestos en formato usable.

Este paso no depende de que nadie haya bajado nada a mano: pide el paquete de
suplementarios a Europe PMC por el identificador del artículo, que es open
access, y de él saca las dos piezas que el estudio necesita. La alternativa
—pedirle al lector que descargue un .xlsx y lo coloque— rompe la cadena en el
primer paso.

De Supplementary Data 1 salen las etiquetas: accesión, nombre en NCBI y especie
revisada, que es el patrón de oro contra el que se contrasta. De Supplementary
Data 2 sale la matriz ANI.

La matriz viene como 540 x 540 en una hoja de cálculo, con las etiquetas de
especie revisada en las columnas de la izquierda. Se escribe un par por fila,
una sola vez por par, para poder unirla con cualquier otro eje por la clave de
accesiones ordenada.

Sus valores son ANIb (pyani con BLASTN), que no es la misma escala que skani ni
que la conversión de la distancia Mash: el objeto de este paso es tener la
columna con la que medir ese desfase en lugar de suponerlo.
"""
import argparse, csv, io, os, re, shutil, sys, urllib.request, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xlsx import read


def fetch(pmcid, work):
    """El paquete de suplementarios del artículo, de Europe PMC."""
    os.makedirs(work, exist_ok=True)
    zp = os.path.join(work, "supp.zip")
    if not os.path.exists(zp):
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/supplementaryFiles"
               % pmcid)
        with urllib.request.urlopen(url, timeout=600) as r, open(zp, "wb") as f:
            shutil.copyfileobj(r, f)
    zipfile.ZipFile(zp).extractall(work)
    got = sorted(f for f in os.listdir(work) if f.endswith(".xlsx"))
    if not got:
        sys.exit("ABORT: el paquete de %s no trae hojas de cálculo" % pmcid)
    print("suplementarios en %s: %d hojas" % (work, len(got)))
    return work


def pick(work, n):
    """MOESM<n> es la numeración con que Springer nombra los suplementarios."""
    for f in sorted(os.listdir(work)):
        if f.endswith(".xlsx") and "MOESM%d_" % n in f:
            return os.path.join(work, f)
    sys.exit("ABORT: no está MOESM%d en %s" % (n, work))


def labels(path, out):
    """Supplementary Data 1: accesión, nombre en NCBI y especie revisada."""
    rows = read(path)
    hdr = rows[1]
    ia = hdr.index("Genome assembly accession")
    io_ = [i for i, h in enumerate(hdr) if h.startswith("Original organism name")][0]
    ir = [i for i, h in enumerate(hdr) if h.startswith("Revised species")][0]
    n = 0
    with open(out, "w") as fh:
        fh.write("accession\tncbi_name\trevised_species\n")
        for r in rows[2:]:
            if len(r) > ia and r[ia].startswith("GC"):
                fh.write("%s\t%s\t%s\n" % (r[ia], r[io_].strip(), r[ir].strip()))
                n += 1
    print("etiquetas escritas: %d genomas" % n)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmcid", default="PMC13106818",
                    help="el artículo en Europe PMC; open access")
    ap.add_argument("--work", required=True, help="dónde dejar los suplementarios")
    ap.add_argument("--labels-out", required=True)
    ap.add_argument("--out", required=True, help="la matriz ANI en formato largo")
    a = ap.parse_args()

    work = fetch(a.pmcid, a.work)
    labels(pick(work, 3), a.labels_out)
    rows = read(pick(work, 4))
    hdr = rows[1]
    ia = hdr.index("Genome assembly accession")
    ir = [i for i, h in enumerate(hdr) if h.startswith("Revised species")][0]
    cols = hdr[ia + 1:]
    data = [r for r in rows[2:] if len(r) > ia and r[ia].startswith("GC")]
    lab = {r[ia]: r[ir].strip() for r in data}

    n = 0
    with open(a.out, "w") as fh:
        fh.write("acc_a\tacc_b\tsp_a\tsp_b\tsame_species\tanib\n")
        for r in data:
            x = r[ia]
            for j, y in enumerate(cols):
                if y not in lab or y <= x:
                    continue
                try:
                    v = float(r[ia + 1 + j])
                except (ValueError, IndexError):
                    continue
                # la hoja mezcla fracción y porcentaje según la celda
                v = v * 100 if v <= 1.0 else v
                fh.write("%s\t%s\t%s\t%s\t%d\t%.4f\n"
                         % (x, y, lab[x], lab[y], lab[x] == lab[y], v))
                n += 1
    print("pares de ANI escritos: %d (esperados %d)"
          % (n, len(data) * (len(data) - 1) // 2))


if __name__ == "__main__":
    main()
