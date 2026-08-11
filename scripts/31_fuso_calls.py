#!/usr/bin/env python3
"""Tabla S10: la clasificación de cada genoma del caso de Fusobacterium.

Se escribe donde el manuscrito la cita, no en este proyecto, para que la ruta de
la tabla suplementaria resuelva a un fichero real. Mientras el caso siga sin
decidirse el fichero está gitignorado junto con la variante del manuscrito; si
entra, se quita una línea del .gitignore y viaja con el repositorio.

Cada fila es una consulta: su cepa tipo más cercana y la distancia a ella, la
especie que este trabajo llama con el corte de 0.043, el nombre que llevaba en
NCBI y el que le asigna la revisión contra la que se contrasta.
"""
import argparse, collections, csv, re

ACC = re.compile(r"(GC[AF]_\d+\.\d+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True, help="genome_labels.tsv del estudio")
    ap.add_argument("--bi", required=True, help="data/bi2026_labels.tsv")
    ap.add_argument("--dist", required=True, help="salida de mash dist tipos x consultas")
    ap.add_argument("--cutoff", type=float, default=0.043)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tipo = {r["accession"]: r["species"]
            for r in csv.DictReader(open(a.labels), delimiter="\t")
            if r["genus"] == "Fusobacterium"}
    bi = {r["accession"]: r for r in csv.DictReader(open(a.bi), delimiter="\t")}
    norm = lambda s: s.split("(")[0].strip()

    best = collections.defaultdict(lambda: (9.0, None))
    for line in open(a.dist):
        f = line.split("\t")
        t, q = ACC.search(f[0]).group(1), ACC.search(f[1]).group(1)
        if t not in tipo or q not in bi or f[4].split("/")[0] == "0":
            continue
        d = float(f[2])
        if d < best[q][0]:
            best[q] = (d, t)

    n = agree = called = 0
    with open(a.out, "w") as fh:
        fh.write("query\tnearest_type_strain\tnearest_type_species\tmash_dist\t"
                 "call\tncbi_name\trevised_species\tagrees_with_revision\n")
        for q, (d, t) in sorted(best.items()):
            call = tipo[t] if d <= a.cutoff else ""
            rev = norm(bi[q]["revised_species"])
            ok = "" if not call else ("1" if call == rev else "0")
            fh.write("%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\n"
                     % (q, t, tipo[t], d, call, norm(bi[q]["ncbi_name"]), rev, ok))
            n += 1
            called += bool(call)
            agree += ok == "1"
    print("escrito %s: %d consultas, %d con llamada a d <= %.3f, %d coinciden (%.2f %%)"
          % (a.out, n, called, a.cutoff, agree, 100.0 * agree / called))


if __name__ == "__main__":
    main()
