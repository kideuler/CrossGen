#!/bin/bash

# Options:
#   --mambo   also mesh the planar faces of every .step/.stp file under
#             MAMBO_DIR (recursively) with Mesh3Dto2Dfacesgmsh, writing
#             data/meshes/mambo/<stem>_<face>.obj -- only the nice faces, and
#             each shape once across all the files: data/meshes/mambo/
#             shapes.txt is the registry the runs share (see the header of
#             src/utils/Mesh3Dto2Dfacesgmsh.cxx).
MAMBO=false
for arg in "$@"; do
	case "$arg" in
		--mambo) MAMBO=true ;;
		*) echo "Unknown option: $arg" >&2; exit 1 ;;
	esac
done

source clean.sh
source build.sh --all

export NP=100
MAMBO_DIR=$HOME/Documents/mambo

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

if [ "$MAMBO" = true ]; then
	EXEC3D="$(pwd)/build/Mesh3Dto2Dfacesgmsh"
	MAMBO_OUT="$OUT_ROOT/mambo"
	if [ ! -x "$EXEC3D" ]; then
		echo "Error: executable not found: $EXEC3D" >&2
		exit 1
	fi
	if [ ! -d "$MAMBO_DIR" ]; then
		echo "Error: MAMBO_DIR not found: $MAMBO_DIR" >&2
		exit 1
	fi
	mkdir -p "$MAMBO_OUT"
	# A model can lose faces between runs, and the registry of shapes already
	# written decides which faces are new, so drop both first.
	rm -f "$MAMBO_OUT"/*.obj "$MAMBO_OUT/shapes.txt"
	while IFS= read -r -d '' step; do
		stem="$(basename "${step%.*}")"
		# The output name is a printf format, so escape any literal '%'. The
		# backslash matters: zsh reads an unescaped leading % in the pattern
		# as "at the end", which appended %% to every name (B0%_2.obj).
		fmt="$MAMBO_OUT/${stem//\%/%%}_%d.obj"
		echo "Meshing mambo/$(basename "$step") with NP=$NP"
		"$EXEC3D" "$step" "$fmt" 0 "$NP" --registry "$MAMBO_OUT/shapes.txt" || {
			echo "Failed to mesh $step" >&2
		}
	done < <(find "$MAMBO_DIR" -type f \( -iname '*.step' -o -iname '*.stp' \) -print0 | sort -z)
fi

echo "All done. Meshes are in $OUT_ROOT"
cd build
cmake ..
