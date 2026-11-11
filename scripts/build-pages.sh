#!/bin/bash
# Build the browser demo and copy only the public site files to a new directory (the GitHub
# Pages site, see .github/workflows/pages.yml). Adapted from egg-solver's scripts/build-pages.sh.
# usage: scripts/build-pages.sh OUTPUT_DIRECTORY
set -euo pipefail
cd "$(dirname "$0")/.."
output=${1:?usage: scripts/build-pages.sh OUTPUT_DIRECTORY}
if [ -e "$output" ]; then
    echo "output already exists: $output (choose a new directory)" >&2
    exit 2
fi
scripts/build-wasm.sh
mkdir -p "$output/dist"
cp wasm/index.html wasm/demo-worker.mjs wasm/omelettehedron.mjs wasm/omelettehedron-thread.mjs \
    wasm/pattern-input.mjs wasm/pattern-picture.mjs wasm/coi-sw.js "$output/"
# (every module the page imports must be listed above: a missing one only shows in a browser)
for m in $(grep -ho "from '\./[^']*'" "$output"/*.mjs "$output/index.html" | sed "s/from '\.\/\(.*\)'/\1/" | sort -u); do
    [ -f "$output/$m" ] || { echo "build-pages: $m is imported but not copied" >&2; exit 1; }
done
cp wasm/dist/omelettehedron.wasm wasm/dist/omelettehedron-threads.wasm \
    wasm/dist/omelettehedron-threads-medium.wasm "$output/dist/"
# (GitHub Pages: no Jekyll processing)
touch "$output/.nojekyll"
printf 'Pages site: %s\n' "$output"
