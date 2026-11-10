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
echo "all checks passed"
