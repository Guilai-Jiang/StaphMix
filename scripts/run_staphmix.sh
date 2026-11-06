#!/usr/bin/env bash
# Run StaphMix from its conda environment without activating it first.
#
#   ./scripts/run_staphmix.sh -q genome.fna -s -c -m -v -a -p sample -o results
#
# Override the environment name with STAPHMIX_ENV if it is not "staphmix".
set -euo pipefail
exec conda run -n "${STAPHMIX_ENV:-staphmix}" --no-capture-output staphmix "$@"
