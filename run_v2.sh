cd "$(dirname "$0")/.."; export OUT_W=$HOME/fw/v2/w OUT_X=$HOME/fw/v2/x SKIP_ALL=1
for q in "$@"; do python3 -I scripts/01_build.py $q > $OUT_W/${q^^}.log 2>&1 & done; wait; tail -qn2 $(for q in "$@"; do echo $OUT_W/${q^^}.log; done) | grep -v DONE
