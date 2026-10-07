#!/usr/bin/env bash
# usage: bash Edit_fna_gtf_TOGA2.sh -g Myotis_myotis_Bat1k -mt KT901455
# usage: bash Edit_fna_gtf_TOGA2.sh -g Rhinolophus_ferrumequinum_Bat1k -mt MK089295
# usage: bash Edit_fna_gtf_TOGA2.sh -g Myotis_emarginatus_Bat1k -mt 
# usage: bash Edit_fna_gtf_TOGA2.sh -g Carollia_perspicillata_Bat1k -mt NC_022422

# safety for failing steps
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

# ---------------------------------
# Clean up outputs from a previous (possibly failed) run
# ---------------------------------
echo ""
echo "Cleaning up previous outputs for ${GENOME} / ${SPNAME}_${ACCESSION}"

rm -f "${GENOME}_replaced.gtf" \
      "${GENOME}_replaced_field1.gtf" \
      "${GENOME}_final.gtf" \
      "${GENOME}_final_filtered.gtf" \
      "${GENOME}_nucleotide_newHeader.fa" \
      "${GENOME}_fasta_headers.txt" \
      "${GENOME}_gtf_seqnames.txt" \
      "${GENOME}_fasta_not_in_gtf.txt" \
      "${GENOME}_gtf_not_in_fasta.txt" \
      "${GENOME}_nucleotide_filtered_mtDNA_forParse.fa" \
      "${GENOME}_nucleotide_filtered_mtDNA_forParse.fa.gz" \
      "${GENOME}_filtered_mtDNA_forParse.gtf" \
      "${GENOME}_filtered_mtDNA_forParse.gtf.gz" \
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
        echo "!! Script failed (exit code ${exit_code}) — cleaning up partial outputs !!" >&2
        rm -f "${GENOME}_replaced.gtf" \
              "${GENOME}_replaced_field1.gtf" \
              "${GENOME}_final.gtf" \
              "${GENOME}_final_filtered.gtf" \
              "${GENOME}_nucleotide_newHeader.fa" \
              "${GENOME}_fasta_headers.txt" \
              "${GENOME}_gtf_seqnames.txt" \
              "${GENOME}_fasta_not_in_gtf.txt" \
              "${GENOME}_gtf_not_in_fasta.txt" \
              "${GENOME}_nucleotide_filtered_mtDNA_forParse.fa" \
              "${GENOME}_nucleotide_filtered_mtDNA_forParse.fa.gz" \
              "${GENOME}_filtered_mtDNA_forParse.gtf" \
              "${GENOME}_filtered_mtDNA_forParse.gtf.gz" \
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

# Step 1: Replace # with _ ONLY in the attributes column ($9)
echo "Step 1: Replace # with _ in attributes"
awk 'BEGIN {FS=OFS="\t"} /^#/ {next} {
    gsub(/#/, "_", $9);
    print $0;
}' "${GENOME}.gtf" > "${GENOME}_replaced.gtf"

# Step 2: Replace $1 for gene lines with transcript_id
echo "Step 2: Replace \$1 for gene lines with transcript_id"
awk 'BEGIN {FS=OFS="\t"}
NR==FNR {
    # First pass: Build transcript_map
    if ($3 == "transcript") {
        match($9, /gene_id "[^"]+"/, gene_arr);
        match($9, /transcript_id "[^"]+"/, transcript_arr);
        if (gene_arr[0] != "" && transcript_arr[0] != "") {
            gene_id = gene_arr[0];
            gsub(/gene_id "/, "", gene_id);
            gsub(/"$/, "", gene_id);

            transcript_id = transcript_arr[0];
            gsub(/transcript_id "/, "", transcript_id);
            gsub(/"$/, "", transcript_id);

            if (!(gene_id in transcript_map)) {
                transcript_map[gene_id] = transcript_id;
            }
        }
    }
    next;
}
{
    # Second pass: Replace $1 for gene lines
    if ($3 == "gene") {
        match($9, /gene_id "[^"]+"/, gene_arr);
        if (gene_arr[0] != "") {
            gene_id = gene_arr[0];
            gsub(/gene_id "/, "", gene_id);
            gsub(/"$/, "", gene_id);

            if (gene_id in transcript_map) {
                $1 = transcript_map[gene_id];
            }
        }
    }
    print $0;
}' "${GENOME}_replaced.gtf" "${GENOME}_replaced.gtf" > "${GENOME}_replaced_field1.gtf"

# Step 2.5: Synthesize missing "gene" lines from "transcript" lines
echo "Step 2.5: Checking for gene lines, synthesizing if absent"

if ! awk -F'\t' '$3=="gene"{f=1} END{exit !f}' "${GENOME}_replaced_field1.gtf"; then
    echo "No 'gene' lines found — synthesizing from 'transcript' lines"

    awk 'BEGIN {FS=OFS="\t"}
    $3=="transcript" {
        match($9, /gene_id "[^"]+"/, gid_arr);
        gene_id = gid_arr[0];
        gsub(/gene_id "/, "", gene_id);
        gsub(/"$/, "", gene_id);

        match($9, /transcript_id "[^"]+"/, tid_arr);
        transcript_id = tid_arr[0];
        gsub(/transcript_id "/, "", transcript_id);
        gsub(/"$/, "", transcript_id);

        # Use the first transcript_id seen for this gene_id, same rule as Step 2
        if (!(gene_id in grep_tid)) grep_tid[gene_id] = transcript_id;

        if (!(gene_id in gmin) || $4 < gmin[gene_id]) gmin[gene_id] = $4;
        if (!(gene_id in gmax) || $5 > gmax[gene_id]) gmax[gene_id] = $5;
        gsource[gene_id] = $2;
        gstrand[gene_id] = $7;
    }
    END {
        for (g in gmin) {
            printf "%s\t%s\tgene\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_biotype \"protein_coding\";\n",
                grep_tid[g], gsource[g], gmin[g], gmax[g], gstrand[g], g;
        }
    }' "${GENOME}_replaced_field1.gtf" > "${GENOME}_synthetic_genes.gtf"

    cat "${GENOME}_replaced_field1.gtf" "${GENOME}_synthetic_genes.gtf" > "${GENOME}_replaced_field1_withgenes.gtf"
    mv "${GENOME}_replaced_field1_withgenes.gtf" "${GENOME}_replaced_field1.gtf"
    rm "${GENOME}_synthetic_genes.gtf"

    echo "Synthesized $(awk -F'\t' '$3=="gene"' "${GENOME}_replaced_field1.gtf" | wc -l) gene lines"
else
    echo "Gene lines already present, skipping synthesis"
fi

# Step 3: Replace $1 with transcript_id for transcript/exon/CDS features AND add gene_biotype
echo "Step 3: Replace \$1 for transcript features and add gene_biotype"
awk 'BEGIN {FS=OFS="\t"} {
    # Replace $1 with transcript_id for features that have it (transcript, exon, CDS, etc.)
    if ($3 != "gene") {
        if (match($9, /transcript_id "[^"]+"/, arr)) {
            transcript_id = arr[0];
            gsub(/transcript_id "/, "", transcript_id);
            gsub(/"$/, "", transcript_id);
            $1 = transcript_id;
        }
    }
    # Ensure gene features have gene_biotype if missing
    else if ($3 == "gene") {
        if ($9 !~ /gene_biotype/) {
            # Add space before gene_biotype if $9 is non-empty
            if ($9 != "") {
                $9 = $9 " gene_biotype \"protein_coding\";";
            } else {
                $9 = "gene_biotype \"protein_coding\";";
            }
        }
    }
    print $0;
}' "${GENOME}_replaced_field1.gtf" > "${GENOME}_final.gtf"

awk 'BEGIN{FS=OFS="\t"} {
    order = ($3=="gene") ? 0 : ($3=="transcript") ? 1 : 2
    print order"\t"$0
}' "${GENOME}_final.gtf" \
| sort -t $'\t' -k2,2 -k1,1n -k5,5n \
| cut -f2- \
> "${GENOME}_final_sorted.gtf"

mv "${GENOME}_final_sorted.gtf" "${GENOME}_final.gtf"

# Change fasta header format: substitute # by _
echo ""
echo "Replace # with _ in FASTA headers"
echo ""
sed 's/#/_/g' "${GENOME}"_nucleotide.fa > "${GENOME}"_nucleotide_newHeader.fa

# Remove empty fasta sequences
echo "Removing empty sequences from multifasta"
awk '
/^>/ {
    if (seq != "") { print header; print seq }
    header = $0; seq = ""; next
}
{ seq = seq $0 }
END {
    if (seq != "") { print header; print seq }
}
' "${GENOME}"_nucleotide_newHeader.fa > "${GENOME}"_nucleotide_newHeader_noempty.fa

n_before=$(grep -c "^>" "${GENOME}"_nucleotide_newHeader.fa)
n_after=$(grep -c "^>" "${GENOME}"_nucleotide_newHeader_noempty.fa)
echo "Removed $((n_before - n_after)) empty sequences ($n_before -> $n_after)"
echo ""

mv "${GENOME}"_nucleotide_newHeader_noempty.fa "${GENOME}"_nucleotide_newHeader.fa

echo "Extract FASTA headers"
echo ""
grep "^>" "${GENOME}"_nucleotide_newHeader.fa | sed 's/>//' | sort -u > "${GENOME}"_fasta_headers.txt

echo "Extract GTF seqnames"
echo ""
awk '{print $1}' "${GENOME}"_final.gtf | sort -u > "${GENOME}"_gtf_seqnames.txt

echo "Compare FASTA headers and GTF seqnames"
echo""
comm -23 "${GENOME}"_fasta_headers.txt "${GENOME}"_gtf_seqnames.txt > "${GENOME}"_fasta_not_in_gtf.txt
comm -13 "${GENOME}"_fasta_headers.txt "${GENOME}"_gtf_seqnames.txt > "${GENOME}"_gtf_not_in_fasta.txt

echo "Print counts"
echo "FASTA headers not in GTF: $(wc -l < "${GENOME}"_fasta_not_in_gtf.txt)"
echo "GTF seqnames not in FASTA: $(wc -l < "${GENOME}"_gtf_not_in_fasta.txt)"

# Remove GTF lines for which there is no corresponding FASTA sequence
echo ""
echo "Remove GTF lines without matching FASTA sequences"
echo ""
if [[ -s "${GENOME}_gtf_not_in_fasta.txt" ]]; then
    awk -F'\t' 'NR==FNR {excl[$0]; next} !($1 in excl)' \
        "${GENOME}_gtf_not_in_fasta.txt" "${GENOME}_final.gtf" \
        > "${GENOME}_final_filtered.gtf"
    echo "Removed $(wc -l < "${GENOME}_gtf_not_in_fasta.txt") GTF lines without matching FASTA sequences."
else
    echo "No GTF lines to remove (all seqnames have matching FASTA headers)."
    cp "${GENOME}_final.gtf" "${GENOME}_final_filtered.gtf"
fi

echo "Final GTF lines: $(wc -l < "${GENOME}_final_filtered.gtf")"

# Get mitochondrial genome
echo ""
echo "Fetching mitochondrial genome: $ACCESSION"
python3 "${SCRIPT}"/get_mitochondrion.py "${ACCESSION}" "${SPNAME}"

cat "${GENOME}"_nucleotide_newHeader.fa "${SPNAME}"_"${ACCESSION}"_genes.fa \
    > "${GENOME}"_nucleotide_filtered_mtDNA_forParse.fa

cat "${GENOME}_final_filtered.gtf" "${SPNAME}"_"${ACCESSION}"_genes.gtf \
    > "${GENOME}"_filtered_mtDNA_forParse.gtf

# Safety net: consolidate any attribute field accidentally split across extra tab-separated columns
awk 'BEGIN{FS=OFS="\t"} {
    if (NF > 9) {
        attr = $9
        for (i=10; i<=NF; i++) attr = attr" "$i
        $9 = attr
        NF = 9
    }
    print
}' "${GENOME}"_filtered_mtDNA_forParse.gtf > "${GENOME}"_filtered_mtDNA_forParse.gtf.tmp
mv "${GENOME}"_filtered_mtDNA_forParse.gtf.tmp "${GENOME}"_filtered_mtDNA_forParse.gtf

# check field number in gtf
echo ""
awk -F'\t' '{count[NF]++} END {
    for (n in count) printf "Number of lines with %d fields: %d\n", n, count[n]
    if (9 in count && length(count) == 1)
        print "OK: all lines have 9 fields"
    else
        print "WARNING: some lines do not have 9 fields"
}' "${GENOME}"_filtered_mtDNA_forParse.gtf | sort -n

# Compress and check fasta
echo ""
echo "Start compression of fa and gtf"
gzip "${GENOME}"_nucleotide_filtered_mtDNA_forParse.fa
gzip -t "${GENOME}"_nucleotide_filtered_mtDNA_forParse.fa.gz && \
    echo "fasta mtDNA gz status: OK" || echo "fasta mtDNA gz status: FAILED"

# Compress and check gtf
echo ""
gzip "${GENOME}"_filtered_mtDNA_forParse.gtf
gzip -t "${GENOME}"_filtered_mtDNA_forParse.gtf.gz && \
    echo "gtf mtDNA gz status: OK" || echo "gtf mtDNA gz status: FAILED"

# Cleaning
rm "${GENOME}"_fasta_not_in_gtf.txt "${GENOME}"_gtf_seqnames.txt "${GENOME}"_fasta_headers.txt "${GENOME}"_gtf_not_in_fasta.txt
rm "${SPNAME}"_"${ACCESSION}".gb
rm "${GENOME}_replaced_field1.gtf" "${GENOME}_replaced.gtf" "${GENOME}_final.gtf"

echo "Finshed!"
