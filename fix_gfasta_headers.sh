#!/usr/bin/env bash
# usage: bash fix_gfasta_headers_v3.sh FILE1.fa.gz [FILE2.fa.gz ...]
#
# Patches already-produced --gfasta FASTA files in place: whitelist-based
# header sanitization. Any character NOT in [A-Za-z0-9._-] (letters,
# digits, dot, underscore, hyphen) is replaced with '_', one-for-one --
# so ',', '$', '(', ')', '"', whitespace, etc. are all neutralized in a
# single pass, without merging or dropping any distinguishing content
# (unlike stripping a whole field, which can collide paralogous loci --
# see the earlier '_<digits>[,<digits>...]'-stripping attempt, which
# produced ~1700 duplicate headers on real data).
#
# Rationale: --gfasta uses the header line as gene_name and (via
# split-pipe's internal, non-CSV-aware parsing of its own results files)
# characters like ',' or '$' inside it can corrupt downstream parsing
# (e.g. molinfo.py join_csv_chunk_files). Rather than chase each new
# offending character one at a time, whitelist what's safe.
#
# Keeps a .bak of each original .gz file. Warns (does not fail) if the
# fix still produces duplicate headers.

set -euo pipefail

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 FILE1.fa.gz [FILE2.fa.gz ...]" >&2
    exit 1
fi

for f in "$@"; do
    if [[ ! -f "$f" ]]; then
        echo "SKIP: $f not found" >&2
        continue
    fi

    echo "Fixing headers in $f"
    tmp_fixed=$(mktemp "${f%.gz}.fixed.XXXXXX")

    zcat "$f" | awk '
    /^>/ {
        header = substr($0, 2)
        gsub(/[^A-Za-z0-9._-]/, "_", header)
        print ">" header
        next
    }
    { print }
    ' > "$tmp_fixed"

    n_dup=$(grep "^>" "$tmp_fixed" | sed 's/>//' | sort | uniq -d | wc -l)
    if [[ "$n_dup" -gt 0 ]]; then
        echo "WARNING: $n_dup duplicate headers after fixing $f (first 10 shown):" >&2
        grep "^>" "$tmp_fixed" | sed 's/>//' | sort | uniq -d | head -10 >&2
    else
        echo "  No duplicate headers."
    fi

    cp "$f" "${f}.bak3"
    gzip -c "$tmp_fixed" > "$f"
    rm "$tmp_fixed"

    n_before=$(zcat "${f}.bak3" | grep -c "^>")
    n_after=$(zcat "$f" | grep -c "^>")
    echo "  Records before: $n_before, after: $n_after (backup: ${f}.bak3)"
    echo ""
done

echo "Done. Remember to re-run the index-building step (SPLITPIPE_INDEX_GFASTA) with the patched FASTA."