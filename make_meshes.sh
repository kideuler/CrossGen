#!/bin/bash

source clean.sh
source build.sh --all

export NP=120

# Geometry lives in per-material-kind subdirectories (singlemat/, multimat/,
# ...); each one is meshed into the matching subdirectory of data/meshes.
GEOM_ROOT="$(pwd)/data/geometry"
OUT_ROOT="$(pwd)/data/meshes"
EXEC="$(pwd)/build/Mesh2Dgmsh"

if [ ! -x "$EXEC" ]; then
	echo "Error: executable not found: $EXEC" >&2
	exit 1
fi

shopt -s nullglob
for geo in "$GEOM_ROOT"/*/*.geo; do
	kind="$(basename "$(dirname "$geo")")"
	out_dir="$OUT_ROOT/$kind"
	mkdir -p "$out_dir"
	echo "Meshing $kind/$(basename "$geo") with NP=$NP"
	"$EXEC" "$geo" "$out_dir" "$NP" || {
		echo "Failed to mesh $geo" >&2
	}
done
shopt -u nullglob

echo "All done. Meshes are in $OUT_ROOT"
cd build
cmake ..
