#!/usr/bin/env bash
# Los tres ejes sobre los mismos genomas: Mash a dos tamaños de sketch y skani.
#
# El objeto es medir el desfase entre escalas en lugar de suponerlo. Bi et al.
# publican ANIb (pyani con BLASTN); nuestra conversión está calibrada contra
# skani; y la distancia Mash es lo que queremos poder usar en su lugar. Con los
# tres sobre los mismos pares, la traducción se mide.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
cd "$FUSO"
cut -f2 data/manifest_540.tsv | tail -n +2 > data/paths.txt
n=$(wc -l < data/paths.txt); echo "genomas: $n"

for s in 10000 100000; do
    [[ -s data/fuso_s${s}.msh ]] || "$MASH" sketch -p 8 -k 21 -s $s -o data/fuso_s${s} -l data/paths.txt >/dev/null 2>&1
    # sin corte de cribado: son pocos genomas y el rango completo interesa
    "$MASH" dist -p 8 data/fuso_s${s}.msh data/fuso_s${s}.msh > data/mash_s${s}.raw
    echo "  s=$s: $(wc -l < data/mash_s${s}.raw) filas"
done

[[ -s data/skani_ani.tsv ]] || "$SKANI" triangle -l data/paths.txt -E -t 8 -s 80 --min-af 5 -o data/skani_ani.tsv >/dev/null 2>&1
echo "  skani: $(( $(wc -l < data/skani_ani.tsv) - 1 )) pares"
