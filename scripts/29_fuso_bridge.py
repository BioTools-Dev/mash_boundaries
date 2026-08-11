#!/usr/bin/env python3
"""El puente entre escalas: ANIb de Bi et al., skani y distancia Mash sobre los mismos pares.

Sin esto, trasladar su umbral de 93.38 % a distancia Mash es una suposición: su
ANI es ANIb (pyani con BLASTN) y nuestra conversión ANI = 1 - 1.12 d se calibró
contra skani. Aquí se mide el desfase entre las tres escalas par a par, y con él
se sitúa su umbral en el eje de Mash en vez de darlo por equivalente.
"""
import argparse, collections, csv, math, os, re, statistics

ACC = re.compile(r"(GC[AF]_\d+\.\d+)")


def load_mash(path, keep):
    d = {}
    for line in open(path):
        f = line.split("\t")
        a, b = ACC.search(f[0]), ACC.search(f[1])
        if not (a and b):
            continue
        a, b = a.group(1), b.group(1)
        if a == b:
            continue
        k = tuple(sorted((a, b)))
        if k in keep and f[4].split("/")[0] != "0":   # sin hash compartido no es medida
            d[k] = float(f[2])
    return d


def load_skani(path, keep):
    d = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, b = ACC.search(f[0]), ACC.search(f[1])
            if not (a and b):
                continue
            k = tuple(sorted((a.group(1), b.group(1))))
            if k in keep:
                d[k] = (float(f[2]), min(float(f[3]), float(f[4])))
    return d


def fit_through_origin(pts):
    """1 - ANI = k d, forzada por el origen, como en el estudio principal."""
    num = sum(d * (1 - a) for d, a in pts)
    den = sum(d * d for d, _ in pts)
    return num / den


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anib", required=True)
    ap.add_argument("--mash", nargs="+", required=True, help="etiqueta=fichero")
    ap.add_argument("--skani", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    anib = {}
    cls = {}
    for r in csv.DictReader(open(a.anib), delimiter="\t"):
        k = (r["acc_a"], r["acc_b"])
        anib[k] = float(r["anib"])
        cls[k] = (r["sp_a"], r["sp_b"], r["same_species"] == "1")
    keep = set(anib)

    mash = {}
    for spec in a.mash:
        tag, path = spec.split("=", 1)
        mash[tag] = load_mash(path, keep)
    sk = load_skani(a.skani, keep)

    with open(a.out + "_pairs.tsv", "w") as fh:
        cols = ["acc_a", "acc_b", "sp_a", "sp_b", "same_species", "anib", "skani", "af_min"]
        cols += ["mash_%s" % t for t in mash]
        fh.write("\t".join(cols) + "\n")
        for k in sorted(keep):
            sa, sb, same = cls[k]
            row = [k[0], k[1], sa, sb, "1" if same else "0", "%.4f" % anib[k]]
            row.append("%.4f" % sk[k][0] if k in sk else "")
            row.append("%.2f" % sk[k][1] if k in sk else "")
            for t in mash:
                row.append("%.6f" % mash[t][k] if k in mash[t] else "")
            fh.write("\t".join(row) + "\n")

    with open(a.out + "_summary.txt", "w") as fh:
        def out(s=""):
            print(s)
            fh.write(s + "\n")

        out("Puente entre escalas sobre %d pares con ANIb" % len(anib))
        out("  con skani: %d   %s" % (len(sk), "  ".join(
            "con mash %s: %d" % (t, len(m)) for t, m in mash.items())))
        out()

        # --- desfase ANIb frente a skani, en la banda que decide
        both = [(anib[k], sk[k][0]) for k in keep if k in sk]
        out("ANIb frente a skani, %d pares con los dos:" % len(both))
        out("  r = %+.4f" % pearson([x for x, _ in both], [y for _, y in both]))
        out("  %-16s %8s %8s %9s" % ("banda de ANIb", "pares", "skani", "desfase"))
        for lo, hi in ((99, 100), (98, 99), (97, 98), (96, 97), (95, 96),
                       (94, 95), (93.5, 94), (93, 93.5), (92, 93), (90, 92), (80, 90)):
            v = [(x, y) for x, y in both if lo <= x < hi]
            if len(v) >= 5:
                out("  %-16s %8d %8.2f %+9.2f"
                    % ("%.1f - %.1f" % (lo, hi), len(v),
                       statistics.median(y for _, y in v),
                       statistics.median(y - x for x, y in v)))
        out()

        # --- la conversión de Mash, reajustada dentro de este género
        for t, m in mash.items():
            pts = [(m[k], sk[k][0] / 100.0) for k in keep if k in m and k in sk]
            if len(pts) < 50:
                continue
            near = [(d, x) for d, x in pts if 0 < d <= 0.13]
            k1 = fit_through_origin(near)
            out("Conversion sobre skani, Mash %s, %d pares con d <= 0.13:" % (t, len(near)))
            out("  ajuste local  ANI = 1 - %.4f d" % k1)
            out("  global del estudio principal: 1.1231")
            dev = max(abs((1 - 1.12 * d) - x) for d, x in near)
            out("  la forma global 1 - 1.12 d se aparta como mucho %.2f puntos aqui"
                % (100 * dev))
            ptsb = [(m[k], anib[k] / 100.0) for k in keep if k in m]
            nearb = [(d, x) for d, x in ptsb if 0 < d <= 0.13]
            out("  contra ANIb en cambio:  ANI = 1 - %.4f d" % fit_through_origin(nearb))
            out()

        # --- donde cae su umbral en el eje de Mash, medido
        out("Donde cae el umbral de Bi et al. en cada eje:")
        for edge in (93.3846, 93.9482):
            row = ["ANIb %.4f" % edge]
            v = [sk[k][0] for k in keep if k in sk and abs(anib[k] - edge) < 0.35]
            if v:
                row.append("skani %.2f" % statistics.median(v))
            for t, m in mash.items():
                w = [m[k] for k in keep if k in m and abs(anib[k] - edge) < 0.35]
                if w:
                    row.append("mash %s d = %.4f" % (t, statistics.median(w)))
            out("  " + "   ".join(row))


if __name__ == "__main__":
    main()
