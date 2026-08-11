#!/usr/bin/env python3
"""¿Es el hueco una propiedad del género o de cuántos genomas se muestrearon?

El hueco que separa las especies tiene dos bordes y se mueven por razones
distintas. El inferior es el par interespecie más parecido: al añadir genomas
solo puede subir, porque aparecen pares nuevos que pueden ser más parecidos que
los que había. El superior es el par intraespecie más divergente: al añadir
genomas solo puede bajar, porque se muestrea más diversidad dentro de cada
especie. Los dos se mueven en la misma dirección —hacia distancias mayores, hacia
identidades menores— así que un hueco medido sobre una colección profunda cae por
fuerza más abajo que uno medido sobre material tipo.

Si eso es lo que ocurre, el umbral que Bi et al. proponen para Fusobacterium no
es una constante del género sino una lectura de su colección, y no se puede
transferir a otro conjunto de genomas ni a otro género.

El experimento no necesita recalcular nada: submuestrea la matriz ya medida a
distintas profundidades por especie y sigue los dos bordes. Como el submuestreo
es aleatorio, cada profundidad se repite y se reportan los cuartiles.
"""
import argparse, collections, csv, random, statistics


def load(path, xcol):
    rows = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if not r[xcol]:
            continue
        rows.append((r["acc_a"], r["acc_b"], r["sp_a"], r["sp_b"], float(r[xcol])))
    sp = {}
    for a, b, sa, sb, _ in rows:
        sp[a] = sa
        sp[b] = sb
    return rows, sp


def edges(rows, keep, higher_is_closer):
    """(borde inferior, borde superior) del hueco sobre el subconjunto `keep`."""
    intra, inter = [], []
    for a, b, sa, sb, v in rows:
        if a in keep and b in keep:
            (intra if sa == sb else inter).append(v)
    if not intra or not inter:
        return None
    if higher_is_closer:                 # identidad: intra alto, inter bajo
        return max(inter), min(intra)
    return max(intra), min(inter)        # distancia: intra bajo, inter alto


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--column", required=True, help="anib, skani o mash_s10000")
    ap.add_argument("--identity", action="store_true",
                    help="la columna es identidad y no distancia")
    ap.add_argument("--depths", default="1,2,3,4,5,7,10,15,20,30,50,0",
                    help="genomas por especie; 0 significa todos")
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows, sp = load(a.pairs, a.column)
    by_sp = collections.defaultdict(list)
    for acc, s in sp.items():
        by_sp[s].append(acc)
    for v in by_sp.values():
        v.sort()
    rng = random.Random(a.seed)
    print("genomas %d, especies %d, pares con %s: %d"
          % (len(sp), len(by_sp), a.column, len(rows)))

    with open(a.out, "w") as fh:
        fh.write("depth\trep\tgenomes\tintra_pairs\tlow_edge\thigh_edge\tgap\n")
        for d in [int(x) for x in a.depths.split(",")]:
            res = []
            reps = 1 if d == 0 else a.reps
            for rep in range(reps):
                keep = set()
                for s, accs in by_sp.items():
                    keep |= set(accs if d == 0 or d >= len(accs)
                                else rng.sample(accs, d))
                e = edges(rows, keep, a.identity)
                if e is None:
                    continue
                lo, hi = e
                ip = sum(1 for x, y, sa, sb, _ in rows
                         if x in keep and y in keep and sa == sb)
                gap = (hi - lo) if a.identity else (hi - lo)
                fh.write("%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n"
                         % (d, rep, len(keep), ip, lo, hi, gap))
                res.append((lo, hi, gap, len(keep), ip))
            if not res:
                print("  profundidad %-4s sin pares intraespecie" % (d or "todos"))
                continue
            f = lambda i: statistics.median(x[i] for x in res)
            print("  profundidad %-5s genomas %4d  pares intra %6d   borde bajo %8.4f   "
                  "borde alto %8.4f   hueco %+.4f%s"
                  % (d or "todos", f(3), f(4), f(0), f(1), f(2),
                     "" if f(2) > 0 else "   SIN HUECO"))


if __name__ == "__main__":
    main()
