#!/bin/sh

# Use the first argument as target directory, default to "../Plots" if not provided
TARGET_DIR="${1:-../Plots}"
EPS_DIR="$TARGET_DIR/eps"
PDF_DIR="$TARGET_DIR/pdf"

# Ensure the pdf directory exists
mkdir -p "$PDF_DIR"

for f in "$EPS_DIR"/*.eps; do
    if [ -f "$f" ]; then
        filename=$(basename "$f" .eps)
        ps2pdf -dEPSCrop "$f" "$PDF_DIR/${filename}.pdf"
    fi
done