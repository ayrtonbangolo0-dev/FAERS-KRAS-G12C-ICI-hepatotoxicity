#!/bin/bash
# Step 2: FDA-recommended deduplication (latest FDA_DT per CASEID, ties -> highest PRIMARYID),
# remove FDA-deleted cases, count unique report-PT pairs in the deduplicated database.
# Usage: bash 02_dedup.sh <workdir containing all_demo.tsv, all_reac.tsv, deleted_caseids.txt>
set -e; S="$(cd "$(dirname "$0")" && pwd)"; cd "${1:?workdir}"; export LC_ALL=C
sort -u deleted_caseids.txt > deleted_sorted.txt
sort -t$'\t' -k2,2 -k3,3n -k1,1n -S 1G -T . all_demo.tsv | awk -F'\t' 'prev!=$2 && NR>1{print keep} {prev=$2; keep=$2"\t"$1} END{print keep}' > case_latest.tsv
join -t$'\t' -v1 <(sort -t$'\t' -k1,1 -S 1G case_latest.tsv) deleted_sorted.txt | cut -f2 | sort -u -S 1G > kept_primaryids.txt
wc -l < all_demo.tsv > n_raw.txt; wc -l < kept_primaryids.txt > n_kept.txt
sort -t$'\t' -k1,1 -S 1500M -T . all_reac.tsv | uniq > reac_sorted.tsv
join -t$'\t' kept_primaryids.txt reac_sorted.tsv > reac_kept.tsv
cut -f2 reac_kept.tsv | sort -S 1G -T . | uniq -c | awk '{c=$1; $1=""; sub(/^ /,""); print $0"\t"c}' > pt_counts.tsv
awk -F'\t' 'NR==FNR{h[tolower($1)]=1; next} (tolower($2) in h){print $1}' "$S/hepato_pts.txt" reac_kept.tsv | sort -u | wc -l > n_hepato_reports.txt
echo "raw=$(cat n_raw.txt) kept=$(cat n_kept.txt) hepato=$(cat n_hepato_reports.txt)"
