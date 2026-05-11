#!/bin/sh

# Use the first argument as target directory, default to "../Plots" if not provided
TARGET_DIR="${1:-../Plots}"
EPS_DIR="$TARGET_DIR/eps"
PNG_DIR="$TARGET_DIR/png"

# Ensure the png directory exists
mkdir -p "$PNG_DIR"

for f in "$EPS_DIR"/*.eps; do
    if [ -f "$f" ]; then
        filename=$(basename "$f" .eps)

        gs -dSAFER -dBATCH -dNOPAUSE \
           -sDEVICE=pngalpha \
           -r600 \
           -sOutputFile="$PNG_DIR/${filename}.png" \
           "$f"
    fi
done