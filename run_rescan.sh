cd "$(dirname "$0")/.."; for q in "$@"; do python3 -I scripts/07_rescan.py $q & done; wait
