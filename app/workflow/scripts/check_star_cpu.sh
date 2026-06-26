#!/usr/bin/env bash
set -euo pipefail

if ! STAR --version >/dev/null 2>&1; then
    {
        echo "***********************************************************************"
        echo ""
        echo "ERROR: STAR failed to execute on this node ($(hostname))."
        echo "It might be an illegal-instruction crash (CPU lacks AVX2/required ISA)."
        echo "CPU flags:"
        grep -o -m1 'avx[^ ]*' /proc/cpuinfo || echo "  no avx flags found"
        echo ""
        echo "***********************************************************************"
    }
    exit 1
fi