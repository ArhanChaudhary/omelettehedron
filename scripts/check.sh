#!/usr/bin/env bash
# Pre-commit checks: formatting, clippy on all targets (warnings are errors), the test suite
# (test profile, so debug assertions are checked) and a scan for runtime assertions.
# usage: scripts/check.sh
set -euo pipefail
cd "$(dirname "$0")/.."
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
# No assert!/assert_eq!/assert_ne! outside #[cfg(test)] items: non-test code must not panic.
# debug_assert! and compile-time `const _: () = assert!(...);` are allowed.
awk '
    FNR == 1 { skip = 0 }
    /^[[:space:]]*#\[cfg\(test\)\]/ { skip = 1; depth = 0; opened = 0; next }
    skip {
        o = gsub(/\{/, "{"); c = gsub(/\}/, "}")
        depth += o - c
        if (o > 0) opened = 1
        if ((opened && depth <= 0) || (!opened && /;[[:space:]]*$/)) skip = 0
        next
    }
    /^[[:space:]]*\/\// { next }
    /const _: \(\) = assert!\(/ { next }
    /(^|[^_[:alnum:]])assert(_eq|_ne)?!/ { print FILENAME ":" FNR ": " $0; bad = 1 }
    END { if (bad) { print "runtime assertions found outside test code" > "/dev/stderr"; exit 1 } }
' $(find src -name '*.rs')

# Backend matrix (as in egg-solver's check.sh): the scalar fallback on this machine, clippy for
# the other targets, the library tests on aarch64 under qemu and on wasm32 under bun's WASI.
# A missing tool is reported as "not checked", never as passed. FTO_CHECK_FAST=1 skips it.
if [ "${FTO_CHECK_FAST:-0}" != 1 ]; then
    not_checked=()
    cargo test --lib --features scalar
    for t in aarch64-apple-darwin aarch64-unknown-linux-gnu wasm32-wasip1 wasm32-unknown-unknown; do
        if rustup target list --installed | grep -qx "$t"; then
            cargo clippy --lib --target "$t" -- -D warnings
        else
            not_checked+=("clippy $t (target not installed)")
        fi
    done
    if command -v qemu-aarch64 >/dev/null && command -v aarch64-linux-gnu-gcc >/dev/null; then
        CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER=aarch64-linux-gnu-gcc \
        CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_RUNNER="qemu-aarch64 -L /usr/aarch64-linux-gnu" \
            cargo test --release --lib --target aarch64-unknown-linux-gnu
    else
        not_checked+=("aarch64 tests (qemu-aarch64 or aarch64-linux-gnu-gcc missing)")
    fi
    bun=${BUN:-$(command -v bun || echo /data/tools/bun/bun)}
    if [ -x "$bun" ]; then
        out=$(RUSTFLAGS="-C target-feature=+simd128" cargo test --release --lib --target wasm32-wasip1 \
            --no-run --message-format=json 2>/dev/null |
            sed -n 's/.*"executable":"\([^"]*\.wasm\)".*/\1/p' | tail -1)
        timeout 600 "$bun" scripts/wasi-run.mjs "$out" --test-threads=1
    else
        not_checked+=("wasm32 tests (bun missing)")
    fi
    for n in "${not_checked[@]}"; do echo "not checked: $n"; done
fi
echo "all checks passed"
