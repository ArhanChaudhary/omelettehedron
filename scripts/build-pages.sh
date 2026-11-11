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
    wasm/coi-sw.js "$output/"
cp wasm/dist/omelettehedron.wasm wasm/dist/omelettehedron-threads.wasm "$output/dist/"
# (GitHub Pages: no Jekyll processing)
touch "$output/.nojekyll"
printf 'Pages site: %s\n' "$output"
