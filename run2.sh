cd "$(dirname "$0")/.."
for q in "$@"; do python3 -I scripts/01_build.py $q > work/q/${q^^}.log 2>&1 & done; wait; tail -qn2 $(for q in "$@"; do echo work/q/${q^^}.log; done)
