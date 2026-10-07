#!/usr/bin/env bash
# usage: bash Edit_fna_gfasta_TOGA2.sh -g Myotis_myotis_Bat1k -mt KT901455
# usage: bash Edit_fna_gfasta_TOGA2.sh -g Rhinolophus_ferrumequinum_Bat1k -mt MK089295
# usage: bash Edit_fna_gfasta_TOGA2.sh -g Myotis_emarginatus_Bat1k -mt ?????
# usage: bash Edit_fna_gfasta_TOGA2.sh -g Carollia_perspicillata_Bat1k -mt NC_022422
#
# Prepares a TOGA2 nucleotide.fa (query_annotation-derived, spliced transcript
# sequences) + mitochondrial gene sequences for split-pipe mkref via --gfasta.
#
# Unlike the previous GTF-based pipeline, this version does NOT build or feed
# a GTF to split-pipe. --gfasta auto-generates one internal "gene"/"exon"
# record per FASTA record, spanning 1..len(seq) -- which is always correct
# by construction, since TOGA2's nucleotide.fa is already the spliced
# transcript sequence, not a genomic-coordinate locus. Feeding TOGA2's
# query_annotation.gtf (genome/scaffold coordinates) alongside nucleotide.fa
# (spliced sequence) silently produced invalid exon coordinates in STAR's
# genome index -- that's the mismatch this rewrite avoids entirely.

set -euo pipefail

GENOME=""
ACCESSION=""

usage() {
    echo "Usage: $0 -g GENOME -mt ACCESSION" >&2
    echo "   -g, --genome           Mandatory. Name of the genome to use. e.g. Myotis_myotis_bat1k" >&2
    echo "   -mt, --mtaccession     Mandatory. Accession number of the mitochondrion genome. e.g. KT901455" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -g|--genome) GENOME="$2"; shift 2 ;;
    -mt|--mtaccession) ACCESSION="$2"; shift 2 ;;
    *)
         echo "Invalid option: $1" >&2
         usage
         ;;
    esac
done

if [[ -z "$GENOME" || -z "$ACCESSION" ]]; then
    echo "Error: -g or -mt missing" >&2
    exit 1
fi

LOGFILE="${GENOME}"_pipeline.log
SPNAME=$(echo "$GENOME" | cut -d'_' -f1,2)
SCRIPT="/dlocal/home/2019013/Data/ViroCaen/scRNABats/Scripts"

FINAL_FASTA="${GENOME}_nucleotide_forGfasta.fa"

# ---------------------------------
# Clean up outputs from a previous (possibly failed) run
# ---------------------------------
echo ""
echo "Cleaning up previous outputs for ${GENOME} / ${SPNAME}_${ACCESSION}"

rm -f "${GENOME}_nucleotide_newHeader.fa" \
      "${GENOME}_nucleotide_dedup.fa" \
      "${GENOME}_fasta_headers.txt" \
      "${GENOME}_fasta_duplicate_headers.txt" \
      "${FINAL_FASTA}" \
      "${FINAL_FASTA}.gz" \
      "${SPNAME}_${ACCESSION}.gb" \
      "${SPNAME}_${ACCESSION}_genes.fa" \
      "${SPNAME}_${ACCESSION}_genes.gtf"
# ---------------------------------

# ---------------------------------
# Cleanup function + trap: remove intermediate/partial files if the script fails
# ---------------------------------
cleanup_on_failure() {
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        echo "" >&2
        echo "!! Script failed (exit code ${exit_code}) -- cleaning up partial outputs !!" >&2
        rm -f "${GENOME}_nucleotide_newHeader.fa" \
              "${GENOME}_nucleotide_dedup.fa" \
              "${GENOME}_fasta_headers.txt" \
              "${GENOME}_fasta_duplicate_headers.txt" \
              "${FINAL_FASTA}" \
              "${FINAL_FASTA}.gz" \
              "${SPNAME}_${ACCESSION}.gb" \
              "${SPNAME}_${ACCESSION}_genes.fa" \
              "${SPNAME}_${ACCESSION}_genes.gtf"
        echo "Partial outputs removed. See ${LOGFILE} for details." >&2
    fi
}
trap cleanup_on_failure EXIT
# ---------------------------------

echo ""
echo "================================"
echo "Pipeline started: $(date)"
echo "GENOME     : $GENOME"
echo "ACCESSION  : $ACCESSION"
echo "SPNAME     : $SPNAME"
echo "================================"
echo ""

# Redirect all output to screen and log file
exec > >(tee -a "$LOGFILE") 2>&1

# Step 1: Clean FASTA headers
# - Header format is 'ENS#GeneName#nombre,nombre,...' -- keep only the first
#   two '#'-delimited fields (ENScode_GeneName); the trailing comma-separated
#   field breaks split-pipe's internal CSV parsing (a comma inside gene_name
#   gets CSV-quoted, and split-pipe's naive .split(",") on results files
#   misparses it downstream, e.g. in molinfo join_csv_chunk_files)
# - collapse any remaining whitespace to '_' as a safety net: --gfasta uses
#   the header line as gene_name, and split-pipe parses only the first
#   whitespace-delimited token as the record name

echo "Step 1: Clean FASTA headers (whitelist [A-Za-z0-9._-], everything else -> _)"
awk '
/^>/ {
    header = substr($0, 2)
    gsub(/[^A-Za-z0-9._-]/, "_", header)
    print ">" header
    next
}
{ print }
' "${GENOME}"_nucleotide.fa > "${GENOME}"_nucleotide_newHeader.fa

# Step 2: Remove empty fasta sequences
echo "Step 2: Removing empty sequences from multifasta"
awk '
/^>/ {
    if (seq != "") { print header; print seq }
    header = $0; seq = ""; next
}
{ seq = seq $0 }
END {
    if (seq != "") { print header; print seq }
}
' "${GENOME}"_nucleotide_newHeader.fa > "${GENOME}"_nucleotide_dedup.fa

n_before=$(grep -c "^>" "${GENOME}"_nucleotide_newHeader.fa)
n_after=$(grep -c "^>" "${GENOME}"_nucleotide_dedup.fa)
echo "Removed $((n_before - n_after)) empty sequences ($n_before -> $n_after)"
echo ""

mv "${GENOME}"_nucleotide_dedup.fa "${GENOME}"_nucleotide_newHeader.fa

# Step 3: Check for duplicate headers
# --gfasta uses each header as a unique chromosome name; duplicates would
# silently collide in the STAR genome index (one overwrites the other)
echo "Step 3: Checking for duplicate FASTA headers"
grep "^>" "${GENOME}"_nucleotide_newHeader.fa | sed 's/>//' | sort > "${GENOME}"_fasta_headers.txt
sort "${GENOME}"_fasta_headers.txt | uniq -d > "${GENOME}"_fasta_duplicate_headers.txt

n_dup=$(wc -l < "${GENOME}"_fasta_duplicate_headers.txt)
if [[ "$n_dup" -gt 0 ]]; then
    echo "WARNING: ${n_dup} duplicate FASTA headers found (first 10 shown):" >&2
    head -n 10 "${GENOME}"_fasta_duplicate_headers.txt >&2
    echo "These will collide as chromosome names in the STAR index -- fix upstream (TOGA2 naming) before proceeding." >&2
    exit 1
else
    echo "No duplicate headers found."
fi
echo ""

# Step 4: Get mitochondrial genome (per-gene sequences only; no GTF needed)
echo ""
echo "Fetching mitochondrial genome: $ACCESSION"
python3 "${SCRIPT}"/get_mitochondrion.py "${ACCESSION}" "${SPNAME}"

# Same header-cleaning safety net applied to mito gene headers
awk '
/^>/ {
    header = substr($0, 2)
    gsub(/#/, "_", header)
    gsub(/[ \t]+/, "_", header)
    print ">" header
    next
}
{ print }
' "${SPNAME}"_"${ACCESSION}"_genes.fa > "${SPNAME}"_"${ACCESSION}"_genes_newHeader.fa
mv "${SPNAME}"_"${ACCESSION}"_genes_newHeader.fa "${SPNAME}"_"${ACCESSION}"_genes.fa

# Step 5: Combine nuclear (TOGA2 spliced transcripts) + mitochondrial gene fasta
echo "Step 5: Combining nuclear + mitochondrial FASTA for --gfasta"
cat "${GENOME}"_nucleotide_newHeader.fa "${SPNAME}"_"${ACCESSION}"_genes.fa \
    > "${FINAL_FASTA}"

# Re-check duplicates across the combined file (nuclear vs mito name clash)
combined_dup=$(grep "^>" "${FINAL_FASTA}" | sed 's/>//' | sort | uniq -d)
if [[ -n "$combined_dup" ]]; then
    echo "WARNING: duplicate headers between nuclear and mitochondrial FASTA:" >&2
    echo "$combined_dup" >&2
    exit 1
fi

echo "Final FASTA records: $(grep -c '^>' "${FINAL_FASTA}")"

# Compress and check fasta
echo ""
echo "Start compression of fasta"
gzip "${FINAL_FASTA}"
gzip -t "${FINAL_FASTA}.gz" && \
    echo "fasta gz status: OK" || echo "fasta gz status: FAILED"

# Cleaning
rm -f "${GENOME}_nucleotide_newHeader.fa" \
      "${GENOME}_fasta_headers.txt" \
      "${GENOME}_fasta_duplicate_headers.txt" \
      "${SPNAME}_${ACCESSION}.gb" \
      "${SPNAME}_${ACCESSION}_genes.fa" \
      "${SPNAME}_${ACCESSION}_genes.gtf"
	  
echo ""
echo "Finished! Use with split-pipe mkref like:"
echo "  split-pipe --mode mkref --genome_name ${GENOME%%_*} --gfasta ${GENOME%%_*} ${FINAL_FASTA}.gz --output_dir GENOME_INDEX --parfile parfile.txt"
