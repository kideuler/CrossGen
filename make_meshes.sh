#!/bin/bash

source clean.sh
source build.sh --all

export NP=30

# Create output directory
OUT_DIR="$(pwd)/data/meshes"
mkdir -p "$OUT_DIR"

# Loop over geometry files and generate meshes
GEOM_DIR="$(pwd)/data/geometry"
EXEC="$(pwd)/build/Mesh2Dgmsh"

if [ ! -x "$EXEC" ]; then
	echo "Error: executable not found: $EXEC" >&2
	exit 1
fi

shopt -s nullglob
for geo in "$GEOM_DIR"/*.geo; do
	echo "Meshing $(basename "$geo") with NP=$NP"
	"$EXEC" "$geo" "$OUT_DIR" "$NP" || {
		echo "Failed to mesh $geo" >&2
	}
done
shopt -u nullglob

echo "All done. Meshes are in $OUT_DIR"
cd build