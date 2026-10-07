#!/usr/bin/env bash
#SBATCH --job-name=ParseNIP                            # Job name
#SBATCH --output=ParseNIP.%A.out                       # File to which stdout will be written
#SBATCH --error=ParseNIP.%A.err                        # File to which stderr will be written
#SBATCH --partition=tcourt                             # Partition
#SBATCH --cpus-per-task=1                              # Number of cpus ask per task
#SBATCH --time=01-00:00                                # Runtime in DD-HH:MM
#SBATCH --mem=10G                                      # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                                # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr            # Email address

# ---------------------------------
# Version v1.4
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 20/08/2026
# v1.1: add --gfasta mode (-m gfasta), fix OUTPUTSPLITPIPE typo
# v1.2: add --outdir choose name of output directory (default results)
# v1.3: validate genome file globs explicitly (fail fast on 0 or >1 matches)
# v1.4: add demultiplexing support (-d/--demultiplex, -k/--kit, default WT);
#       new 'demultiplex' step to run demultiplexing on its own
# ---------------------------------

# ---------------------------------
# Submission:
# sbatch ParseNIP.sh -r $RUN -g $REFGENOME -s index|mapping|all|demultiplex [-m gtf|gfasta] [-o outputdir] [-d] [-k kit]
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# Variable
GROUP="ViroCaen"
PROJECT="scRNABats"
RUN=""
REFGENOME=""
STEP=""
MODE="gtf"   # gtf (default, --fasta/--gtf) or gfasta (--gfasta, no GTF)
OUTPUTDIR="results"
DEMULTIPLEX=false
KIT="WT"

usage() {
    echo "Usage: $0 -r RUN -g REFGENOME -s STEP [-m gtf|gfasta] [-o outputdir] [-d] [-k kit]" >&2
	echo "   -r, --run           Mandatory. correpond to the list of your fastq" >&2
	echo "   -g, --genome        Mandatory for -s index/mapping/all. Name of the reference genome to use e.g. Myotis_myotis" >&2
	echo "   -s, --step          Mandatory. Step of analyses to run: 'index', 'mapping', 'all', or 'demultiplex'" >&2
	echo "   -m, --mode          Optional. Genome index build mode: 'gtf' (default, --fasta/--gtf) or 'gfasta' (--gfasta, no GTF; TOGA2-derived pre-extracted FASTA)" >&2
    echo "   -o, --outdir        Optional. Name of the output directory (default: results)" >&2
    echo "   -d, --demultiplex   Optional. Also demultiplex raw FASTQ by sample well group (fastq_sep_groups.py); combines with -s index/mapping/all, or use -s demultiplex on its own" >&2
    echo "   -k, --kit           Optional. Kit name for demultiplexing (default: WT)" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "   # Index + mapping + demultiplexing, all together" >&2
    echo "   sbatch $0 -r Test -g Myotis_myotis -s all -m gfasta -d" >&2
    echo "" >&2
    echo "   # Demultiplexing only, no index/mapping" >&2
    echo "   sbatch $0 -r Test -s demultiplex -k WT" >&2
    echo "" >&2
    echo "   # Index + mapping only, no demultiplexing" >&2
    echo "   sbatch $0 -r Test -g Myotis_myotis -s all -m gfasta" >&2
    echo "" >&2
    echo "   # Mapping only, index already built" >&2
    echo "   sbatch $0 -r Test -g Myotis_myotis -s mapping" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -r|--run) RUN="$2"; shift 2 ;;
	-g|--genome) REFGENOME="$2"; shift 2 ;;
	-s|--step)  STEP="$2"; shift 2 ;;
	-m|--mode)  MODE="$2"; shift 2 ;;
    -o|--outdir)  OUTPUTDIR="$2"; shift 2 ;;
    -d|--demultiplex) DEMULTIPLEX=true; shift 1 ;;
    -k|--kit)  KIT="$2"; shift 2 ;;
    *) 
         echo "Invalid option: $1" >&2
         usage
         ;;
    esac
done

if [[ -z "$RUN" || -z "$STEP" ]]; then
    echo "Error missing -r or -s" >&2
    exit 1
fi

if [[ "$STEP" != "demultiplex" && -z "$REFGENOME" ]]; then
    echo "Error missing -g (required for -s index/mapping/all)" >&2
    exit 1
fi

if [[ "$MODE" != "gtf" && "$MODE" != "gfasta" ]]; then
    echo "Error: -m/--mode must be 'gtf' or 'gfasta' (got '$MODE')" >&2
    exit 1
fi

if [[ "$STEP" != "index" && "$STEP" != "mapping" && "$STEP" != "all" && "$STEP" != "demultiplex" ]]; then
    echo "Error: -s/--step must be 'index', 'mapping', 'all', or 'demultiplex' (got '$STEP')" >&2
    exit 1
fi

# -s demultiplex implies demultiplexing, even without -d explicitly
if [[ "$STEP" == "demultiplex" ]]; then
    DEMULTIPLEX=true
fi

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PARSENIP="/soft/2019013/Logiciels/ParseNIP"
GENOMEINDEX="/dlocal/home/2019013/Databases/SplitpipeIndex"
SCRIPTS="/${FOLDER}/Scripts"
INPUTILLUMINA="/${FOLDER}/InputIllumina"
REFGENOMEPATH="/${FOLDER}/ReferenceGenomes"
OUTPUTSPLITPIPE="/${FOLDER}/OutputSplitpipe"
OUTPUTDEMULTIPLEXING="/${FOLDER}/OutputDemultiplexing"

# Resolve a glob pattern to exactly one file, or fail with a clear message.
# Usage: resolve_one_file "description" "glob-pattern"
# Echoes the resolved path on success.

resolve_one_file() {
    local desc="$1"
    local pattern="$2"
    local matches=()

    shopt -s nullglob
    matches=( $pattern )
    shopt -u nullglob

    if [[ "${#matches[@]}" -eq 0 ]]; then
        echo "Error: no file found for $desc (pattern: $pattern)" >&2
        exit 1
    elif [[ "${#matches[@]}" -gt 1 ]]; then
        echo "Error: multiple files match $desc (pattern: $pattern), expected exactly one:" >&2
        printf '  %s\n' "${matches[@]}" >&2
        exit 1
    fi

    echo "${matches[0]}"
}

# Genome index build args, depending on -m/--mode
# gtf mode    : output of the old Edit_fna_gtf_TOGA2.sh    -> ${REFGENOME}_Bat1k_nucleotide_filtered_mtDNA_forParse.fa.gz / .gtf.gz
# gfasta mode : output of Edit_fna_gfasta_TOGA2.sh          -> ${REFGENOME}_Bat1k_nucleotide_forGfasta.fa.gz
# Only needed for 'index'/'all' -- 'mapping' uses an existing --genome_dir,
# and 'demultiplex' does not touch the genome at all.
if [[ "$STEP" == "index" || "$STEP" == "all" ]]; then
    if [[ "$MODE" == "gtf" ]]; then
        GENOME_FASTA=$(resolve_one_file "reference FASTA" "${REFGENOMEPATH}/${REFGENOME}_*.fa.gz")
        GENOME_GTF=$(resolve_one_file "reference GTF" "${REFGENOMEPATH}/${REFGENOME}_*.gtf.gz")
        GENOME_BUILD_ARGS=(
            --fasta "$GENOME_FASTA"
            --gtf "$GENOME_GTF"
            --genome_name "${REFGENOME}"
        )
    else
        GENOME_GFASTA=$(resolve_one_file "reference gfasta FASTA" "${REFGENOMEPATH}/TOGA2/${REFGENOME}*_forGfasta.fa.gz")
        GENOME_BUILD_ARGS=(
            --gfasta "$GENOME_GFASTA"
            --genome_name "${REFGENOME}"
        )
    fi
fi

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
if [[ "$STEP" == "index" || "$STEP" == "all" ]]; then
    echo "genome index build mode: $MODE"
fi
echo "demultiplexing: $DEMULTIPLEX (kit: $KIT)"

# Demultiplex args, appended to index/mapping/all nextflow calls when -d is set
DEMUX_ARGS=()
FASTQ_SAMPLIST_ARGS=()
if [[ "$DEMULTIPLEX" == true ]]; then
    DEMUX_ARGS=( --demultiplex --kit "$KIT" )
    FASTQ_SAMPLIST_ARGS=( --fastq "${LOCAL_WORK_DIR}/fastqs/" --samp_list "${SCRIPTS}/List_Samples_${RUN}.txt" )
fi
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load tools/nextflow/25.10.4 py_env/miniconda3/25.7.0
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${INPUTILLUMINA}"/"${RUN}"/*.fq.gz "${LOCAL_WORK_DIR}"

cd "${LOCAL_WORK_DIR}"

mkdir -p fastqs/

mv *.fq.gz fastqs/

echo "${LOCAL_WORK_DIR}"
# ---------------------------------

# ---------------------------------
# splitpipe
# ---------------------------------
if [[ "${STEP}" == "index" ]]; then

   echo "Build index of one reference genome"

   nextflow run \
   "${PARSENIP}"/main.nf \
   -config "${PARSENIP}"/nextflow_austral.config \
   "${GENOME_BUILD_ARGS[@]}" \
   "${FASTQ_SAMPLIST_ARGS[@]}" \
   "${DEMUX_ARGS[@]}" \
   -profile austral \
   --outdir "${OUTPUTDIR}"/

elif [[ "${STEP}" == "mapping" ]]; then

    echo "Starting mapping and seurat"

   nextflow run \
   "${PARSENIP}"/main.nf \
   -config "${PARSENIP}"/nextflow_austral.config \
   -profile austral \
   --samp_list "${SCRIPTS}/List_Samples_${RUN}.txt" \
   --fastq "${LOCAL_WORK_DIR}/fastqs/" \
   --genome_dir "${GENOMEINDEX}"/"${REFGENOME}"_index/ \
   "${DEMUX_ARGS[@]}" \
   --outdir "${OUTPUTDIR}"/

elif [[ "${STEP}" == "all" ]]; then

    echo "Starting full pipiline, index, mapping and seurat"

   nextflow run \
   "${PARSENIP}"/main.nf \
   -config "${PARSENIP}"/nextflow_austral.config \
   -profile austral \
   --fastq "${LOCAL_WORK_DIR}/fastqs/" \
   --samp_list "${SCRIPTS}/List_Samples_${RUN}.txt" \
   "${GENOME_BUILD_ARGS[@]}" \
   "${DEMUX_ARGS[@]}" \
   --outdir "${OUTPUTDIR}"/

elif [[ "${STEP}" == "demultiplex" ]]; then

    echo "Starting demultiplexing only"

   nextflow run \
   "${PARSENIP}"/main.nf \
   -config "${PARSENIP}"/nextflow_austral.config \
   -profile austral \
   --fastq "${LOCAL_WORK_DIR}/fastqs/" \
   --samp_list "${SCRIPTS}/List_Samples_${RUN}.txt" \
   --demultiplex --kit "$KIT" \
   --outdir "${OUTPUTDIR}"/
fi
# ---------------------------------

# ---------------------------------
# Move Outputs
# ---------------------------------
echo "Move outputs files"

# Destination for the full output directory:
# - 'demultiplex' alone (no genome/mapping involved) goes ONLY to
#   OutputDemultiplexing, there's nothing else to send to OutputSplitpipe.
# - everything else (index/mapping/all, with or without -d) keeps going
#   to OutputSplitpipe as before.
if [[ "${STEP}" == "demultiplex" ]]; then
    DEST="${OUTPUTDEMULTIPLEXING}/${RUN}/"

    if [[ ! -d "${DEST}" ]]; then
        mkdir -p "${DEST}"
    fi

    cp -r "${OUTPUTDIR}"/demultiplexing/* "${DEST}"
else
    DEST="${OUTPUTSPLITPIPE}/${RUN}/"

    # -d was combined with index/mapping/all: the demultiplexing/ subfolder
    # landed inside OUTPUTDIR alongside everything else. Pull it out to
    # OutputDemultiplexing FIRST, then remove it from the local copy so it
    # doesn't also get duplicated into OutputSplitpipe below.
    if [[ -d "${OUTPUTDIR}/demultiplexing" ]]; then
        DEMUX_DEST="${OUTPUTDEMULTIPLEXING}/${RUN}/"

        if [[ ! -d "${DEMUX_DEST}" ]]; then
            mkdir -p "${DEMUX_DEST}"
        fi

        cp -r "${OUTPUTDIR}"/demultiplexing/* "${DEMUX_DEST}"
        rm -r "${OUTPUTDIR}/demultiplexing"
    fi

    if [[ ! -d "${DEST}" ]]; then
        mkdir -p "${DEST}"
    fi

    cp -r "${OUTPUTDIR}" "${DEST}"
fi

echo "Done."
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/ParseNIP.*.err "${FOLDER}/Scripts/Sbatch_log"
# ---------------------------------