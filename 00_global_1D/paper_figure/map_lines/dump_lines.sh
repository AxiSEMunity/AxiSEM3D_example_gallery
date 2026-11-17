#!/bin/bash
# Dump shorelines (GSHHG level 1 only: no lakes, no rivers) and political borders (national, state) as ASCII polylines.
set -e
cd "$(dirname "$0")"     # needs GMT 6 (gmt coast)
R=-R-140/-50/12/62
for res in l i; do
  gmt coast $R -D$res -A0/0/1 -W -M > shore_$res.txt
  gmt coast $R -D$res -N1 -M > border_national_$res.txt
  gmt coast $R -D$res -N2 -M > border_state_$res.txt
done
gmt --version
wc -l *.txt
