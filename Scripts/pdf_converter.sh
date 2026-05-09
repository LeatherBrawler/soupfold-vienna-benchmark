#!/bin/sh
for f in ../Plots/eps/*.eps; do
    if [ -f "$f" ]; then
        filename=$(basename "$f" .eps)
        ps2pdf -dEPSCrop "$f" "../Plots/pdf/${filename}.pdf"
    fi
done