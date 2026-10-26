#!/bin/bash
# BIP version: the single place it is defined. Bump it for each release; it is written
# to publication_metrics.txt so every run records the version that produced it.
BIP_VERSION="1.1.0"
# Usage: bash BIP.sh --platform <illumina_pe|pacbio|nanopore|other> [options] <RunName>
# Run with no arguments (other than --platform/<RunName>) to use the defaults below.

export PATH="$PATH"

# ---------------------------
# Set user-specific arguments (all overridable via --flag or env var; flags win)
# ---------------------------

cores_to_leave="${cores_to_leave:-3}"
wkdir="${wkdir:-/BIP/Barcoding}"
scripts_directory="${scripts_directory:-/BIP/SCRIPTS}"
reference_library_directory="${reference_library_directory:-/BIP/REFS}"
ref_seq_corr="${ref_seq_corr:-}"

# Where the input reads and the parameters workbook live. Both default to the working
# directory (--wd). --fastq takes a directory holding the .fastq/.fastq.gz files; --params
# takes the .xlsx file itself, or a directory holding exactly one .xlsx. Relative paths are
# taken from the directory BIP is launched in. Reads in a --fastq directory other than the
# working directory are only read, never modified.
fastq_dir="${fastq_dir:-}"
param_path="${param_path:-}"

# Taxonomic rank at which the user-supplied identification is enforced during barcode
# inference. An OTU whose taxonomy mismatches the specimen at this rank, or at any rank
# above it, cannot be designated as the target barcode. Ranks BELOW this one do not
# prevent target designation; they are used to break ties between candidate OTUs.
# Ranks use the same letters as SINTAX: k, p, c, o, f, g, s. Default is o (order).
min_rank_match="${min_rank_match:-o}"

platform=""
minreads=5
runid=""

# SINTAX taxonomy assignment. FALSE (default): one thread with a fixed random seed, so the
# same input always gives the same taxonomy and probabilities. TRUE: all available cores
# and no fixed seed, which is faster but lets results vary slightly from run to run.
parallel_sintax="${parallel_sintax:-FALSE}"
# Seed used when parallel_sintax is FALSE.
sintax_randseed=1

#~#~#~#~#~#~#~#~#~#~#
# Advanced parameters
#~#~#~#~#~#~#~#~#~#~#

# This value is used as a multiplier with UMI lengths during demultiplexing. We recommend 0.75 for UMIs greater than or equal to 12 nucleotides, and 1.0 for UMIs with fewer than 12 nucleotides.
umi_overlap_min=0.75
# This value is used as a multiplier with primer lengths. We recommend 0.75.
primer_overlap_min=0.75
# Percentage of mismatches allowed in the forward/reverse UMIs. For UMIs shorter than 8 bp, we recommend 0.0.
error_umi1=0.125
error_umi2=0.125
# Percentage of mismatches allowed in the forward/reverse primers.
error_primer1=0.2
error_primer2=0.2
# Controls how "satellite OTUs" are removed. Default will result in a min read of 5.
minreadssatellite=0.00000625
# Minimum overlap needed for a BIN to match to an OTU. Default is 75% of expected amplicon size.
minoverlap=0.75
# Parameters for VSEARCH's uchime_denovo command (chimera removal).
chimera_abskew=10
chimera_mindiv=0.0005
# Parameters for VSEARCH's usearch_global command (BIN matching).
bin_maxhits=1
bin_maxaccepts=1
# In-frame sequence lengths accepted for the protein-coding COI-5P marker. An OTU whose
# length is not in this list is flagged as an indel (and typed as a NUMT), which makes it
# ineligible as a target barcode. Comma-separated integers, passed to bip3. Update this
# list as more becomes known about legitimate COI-5P barcode lengths.
coi_expected_lengths=640,643,646,649,652,655,658,661,664,667,670,673

# ---------------------------
# Parse arguments properly
# ---------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        -platform|--platform)
            platform="$2"
            if [[ "$platform" != "illumina_pe" && "$platform" != "pacbio" && "$platform" != "nanopore" && "$platform" != "other" ]]; then
                echo "ERROR: Invalid platform '$platform'. Must be one of: illumina_pe, pacbio, nanopore, other"
                exit 1
            fi
            shift 2
            ;;
        --minreads)
            minreads="$2"
            shift 2
            ;;
        --fastq)
            fastq_dir="$2"
            shift 2
            ;;
        --params)
            param_path="$2"
            shift 2
            ;;
        --min_rank_match)
            min_rank_match="$2"
            shift 2
            ;;
        --wd)
            wkdir="$2"
            shift 2
            ;;
        --scripts)
            scripts_directory="$2"
            shift 2
            ;;
        --refs)
            reference_library_directory="$2"
            shift 2
            ;;
        --ref_seq_corr)
            ref_seq_corr="$2"
            shift 2
            ;;
        --cores_to_leave)
            cores_to_leave="$2"
            shift 2
            ;;
        --umi_overlap_min)
            umi_overlap_min="$2"
            shift 2
            ;;
        --primer_overlap_min)
            primer_overlap_min="$2"
            shift 2
            ;;
        --error_umi1)
            error_umi1="$2"
            shift 2
            ;;
        --error_umi2)
            error_umi2="$2"
            shift 2
            ;;
        --error_primer1)
            error_primer1="$2"
            shift 2
            ;;
        --error_primer2)
            error_primer2="$2"
            shift 2
            ;;
        --minreadssatellite)
            minreadssatellite="$2"
            shift 2
            ;;
        --minoverlap)
            minoverlap="$2"
            shift 2
            ;;
        --chimera_abskew)
            chimera_abskew="$2"
            shift 2
            ;;
        --chimera_mindiv)
            chimera_mindiv="$2"
            shift 2
            ;;
        --bin_maxhits)
            bin_maxhits="$2"
            shift 2
            ;;
        --bin_maxaccepts)
            bin_maxaccepts="$2"
            shift 2
            ;;
        --coi_expected_lengths)
            coi_expected_lengths="$2"
            shift 2
            ;;
        --parallel_sintax)
            parallel_sintax="$2"
            shift 2
            ;;
        -*)
            echo "ERROR: Unknown option $1"
            echo "Usage: bash BIP.sh --platform <illumina_pe|pacbio|nanopore|other> [options] <RunName>"
            exit 1
            ;;
        *)
            if [[ -z "$runid" ]]; then
                runid="$1"
            else
                echo "ERROR: Multiple run names provided."
                exit 1
            fi
            shift
            ;;
    esac
done

if [[ -z "$runid" ]]; then
    echo "ERROR: No run name provided."
    echo "Usage: bash BIP.sh --platform <illumina_pe|pacbio|nanopore|other> [options] <RunName>"
    exit 1
fi

if [[ -z "$platform" ]]; then
    echo "ERROR: --platform is required. Must be one of: illumina_pe, pacbio, nanopore, other"
    exit 1
fi

# Accept either the SINTAX rank letter or the full rank name, any case; store the letter.
min_rank_match="$(printf '%s' "$min_rank_match" | tr '[:upper:]' '[:lower:]')"
case "$min_rank_match" in
    k|kingdom) min_rank_match=k ;;
    p|phylum)  min_rank_match=p ;;
    c|class)   min_rank_match=c ;;
    o|order)   min_rank_match=o ;;
    f|family)  min_rank_match=f ;;
    g|genus)   min_rank_match=g ;;
    s|species) min_rank_match=s ;;
    *)
        echo "ERROR: Invalid --min_rank_match '$min_rank_match'. Must be one of: k, p, c, o, f, g, s"
        exit 1
        ;;
esac

case "${parallel_sintax^^}" in
    TRUE|T)  parallel_sintax=TRUE ;;
    FALSE|F) parallel_sintax=FALSE ;;
    *)
        echo "ERROR: Invalid --parallel_sintax '$parallel_sintax'. Must be TRUE or FALSE"
        exit 1
        ;;
esac

if [[ ! "$coi_expected_lengths" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
    echo "ERROR: coi_expected_lengths must be a comma-separated list of integers (got '$coi_expected_lengths')"
    exit 1
fi

# Bounds of the in-frame COI-5P window, used to gate the steps that only make sense for a
# full-length barcode. Derived here so coi_expected_lengths stays the single definition.
coi_min_expected=$(printf '%s' "$coi_expected_lengths" | tr ',' '\n' | sort -n | head -1)
coi_max_expected=$(printf '%s' "$coi_expected_lengths" | tr ',' '\n' | sort -n | tail -1)

# Resolve --fastq / --params to absolute paths now, while a relative path still means the
# directory BIP was launched from (the script changes into the working directory below).
if [[ -n "$fastq_dir" ]]; then
    if [[ ! -d "$fastq_dir" ]]; then
        echo "ERROR: --fastq must be a directory containing the .fastq/.fastq.gz files (got '$fastq_dir')"
        exit 1
    fi
    fastq_dir="$(cd "$fastq_dir" && pwd)"
fi
if [[ -n "$param_path" ]]; then
    if [[ -d "$param_path" ]]; then
        param_path="$(cd "$param_path" && pwd)"
    elif [[ -f "$param_path" ]]; then
        param_path="$(cd "$(dirname "$param_path")" && pwd)/$(basename "$param_path")"
    else
        echo "ERROR: --params not found: '$param_path' (give the .xlsx file, or a directory holding one)"
        exit 1
    fi
fi

cores=$(($(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu) - $cores_to_leave))
if [[ "$parallel_sintax" == "TRUE" ]]; then
    sintax_opts=(--threads "$cores")
else
    sintax_opts=(--threads 1 --randseed "$sintax_randseed")
fi

cd "$wkdir"

# ---------------------------
# Start duration log
# ---------------------------

START_TIME=$(date +%s)
START_LABEL=$(date "+%Y-%m-%d %H:%M:%S")
printf "%-30s %-22s %-22s %-15s\n" "Stage" "Start" "End" "Duration"  > duration_log.txt
printf "%-30s %-22s %-22s %-15s\n" "-----" "-----" "---" "--------" >> duration_log.txt

# Derive behaviour flags from platform
do_pear=0
chimera_early=0
[[ "$platform" == "illumina_pe" ]] && do_pear=1 && chimera_early=1
[[ "$platform" == "pacbio" ]] && chimera_early=1

runid=$(echo "$runid" \
    | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' \
    | sed 's/ /_/g' \
    | sed 's/[^[:alnum:]_.-]//g')

counts_file="$wkdir/${runid}_readcounts.tsv"
echo -e "runid\tstage\tprimer_pair\treads" > "$counts_file"

rescue_log="$wkdir/${runid}_symmetrical_umi_rescue_log.txt"

# ---------------------------
# PARAMETERS SPREADSHEET: COLUMN NAMES
# ---------------------------
#
# THE ONLY PLACE THAT KNOWS THE LAYOUT OF THE PARAMETERS SPREADSHEET.
#
# Every column is found by looking its name up in the header row, so nothing
# below refers to a spreadsheet column by number. Columns may be added, removed
# or reordered in either tab without touching a line of code.
#
#   * RENAMED a column in the spreadsheet?  Change its value here.
#   * ADDED a column?                       Add a COLNAME_*/DICTCOL_* line here,
#                                           then use $COL_<NAME> where you need it.

# Sheet (tab) names
export SHEET_UMIS_PRIMERS="UMIs and Primers"
export SHEET_DICT_UPDATE="Dictionary Update"

# "UMIs and Primers" tab
export COLNAME_PLATE="Plate"
export COLNAME_WELL="Well"
export COLNAME_SAMPLE="Sample"
export COLNAME_FWD_UMI="Forward UMI"
export COLNAME_REV_UMI="Reverse UMI"
export COLNAME_FWD_PRIMER_NAME="Forward Primer Name"
export COLNAME_REV_PRIMER_NAME="Reverse Primer Name"
export COLNAME_NEG_CTRL="Negative Control"

# "Dictionary Update" tab (and dictionary_BIP.tsv, which shares these headings).
# The first two must be the primer-name pair used to key the dictionary.
export DICTCOL_FWD_PRIMER_NAME="Forward Primer Name"
export DICTCOL_REV_PRIMER_NAME="Reverse Primer Name"
export DICTCOL_FWD_PRIMER_SEQ="Forward Primer Sequence"
export DICTCOL_REV_PRIMER_SEQ="Reverse Primer Sequence"
export DICTCOL_MARKER="Marker"
export DICTCOL_REFLIB="Reference Library"
export DICTCOL_MINLEN="Min amplicon length"
export DICTCOL_MAXLEN="Max amplicon length"
export DICTCOL_AMPLEN="Target amplicon length"
export DICTCOL_OTU_PRIMARY="Primary OTU clustering threshold"
export DICTCOL_OTU_SECONDARY="Secondary OTU clustering threshold"
export DICTCOL_TAX_PROB="Tax assign probability threshold"
export DICTCOL_BIN_THRESH="BIN assign threshold"

# ---------------------------
# Define functions
# ---------------------------

# Find the 1-based position of a named column in a TSV's header row.
# Fails loudly rather than silently reading the wrong column.
resolve_col() {
    local wanted="$1" file="$2" idx
    idx=$(awk -F'\t' -v want="$wanted" '
        NR==1 {
            for (i = 1; i <= NF; i++) {
                h = $i
                gsub(/\r/, "", h)
                gsub(/^[ \t]+|[ \t]+$/, "", h)
                if (h == want) { print i; exit }
            }
        }' "$file")
    if [[ -z "$idx" ]]; then
        echo "ERROR: column \"$wanted\" not found in $file" >&2
        echo "       Header is: $(head -1 "$file")" >&2
        return 1
    fi
    printf '%s' "$idx"
}

log_count() {
    local stage="$1"
    local primer="$2"
    local file="$3"

    [[ ! -f "$file" ]] && return

    local count
    count=$(seqkit stats -T "$file" | awk 'NR==2 {print $4}')

    echo -e "${runid}\t${stage}\t${primer}\t${count}" >> "$counts_file"
}

get_length_bounds() {
    local primer="$1"

    awk -F'\t' -v primer="$primer" '
    function trim(x) {
        gsub(/\r/, "", x)
        gsub(/^[ \t]+|[ \t]+$/, "", x)
        return x
    }

    NR==1 {
        for (i=1; i<=NF; i++) h[$i]=i
        next
    }

    {
        f = trim($h["Forward Primer Name"])
        r = trim($h["Reverse Primer Name"])
        key = f "_" r

        if (key == primer) {
            min = $h["Min amplicon length"]
            max = $h["Max amplicon length"]

            print min "\t" max
            exit
        }
    }
    ' "$dictfile"
}

filter_by_length() {
    local primer="$1"
    local infile="$2"
    local outfile="$3"

    read min max < <(get_length_bounds "$primer")

    seqkit seq -m "$min" -M "$max" "$infile" > "$outfile"

    log_count "LengthFilter" "$primer" "$outfile"
}

sanitize_headers() {
    local infile="$1"
    local outfile="$2"

    seqkit replace \
        -p '^([^[:space:]_]+).*' \
        -r '$1' \
        "$infile" > "$outfile"
}

log_step() {
    local stage="$1"
    local start_label="$2"
    local end_label="$3"
    local elapsed_sec="$4"
    local dur="$(( elapsed_sec / 3600 ))h $(( (elapsed_sec % 3600) / 60 ))m"
    printf "%-30s %-22s %-22s %-15s\n" "$stage" "$start_label" "$end_label" "$dur" >> $wkdir/duration_log.txt
}

# Primer pairs whose BIN matching failed. The run still finishes, then reports them and
# exits with an error (see the end of the script).
bin_failures=()
# Primer pairs where almost no reads matched a UMI pair, repeated at the end of the run.
umi_warnings=()

process_all_primers() {
    for d in */; do
        name="${d%/}"
        fq="$PWD/${name}/${name}.fastq"   # absolute path
        [[ ! -f "$fq" ]] && continue
        process_primer_dir "$name" "$fq"
    done
}

# Per-sample OTU work. Runs inside the sample's own private folder (see sampletask below),
# so every file listing here sees only this sample's files. Records how it ended in
# <sampleid>.status: ok, no_otus, or error: <reason>.
sampletask_core() {

    local fasta="$1"
    local sampleid="${fasta%.fasta}"

    if ! grep -q '^>' "${fasta}"; then
        echo "No OTUs remaining for $sampleid after primary clustering"
        echo "no_otus" > "${sampleid}.status"
        return
    fi

    # ---------------------------
    # Primer trimming
    # ---------------------------

    echo -e "\n#################################################################"
    echo -e "   Trimming primers from $sampleid"
    echo -e "#################################################################\n"

    awk '
    /^>/ {
        print
        next
    }

    {
        gsub(/[a-z]/, "", $0)
        print
    }
    ' "$fasta" > "${sampleid}_trimmed.fasta"

    mv "${sampleid}_trimmed.fasta" "$fasta"


    echo -e "\n#################################################################"
    echo -e "   Dereplicating $sampleid reads prior to chimera screening"
    echo -e "#################################################################\n"

    # ---------------------------
    # Dereplicate full dataset
    # ---------------------------

    vsearch --derep_fulllength "$fasta" \
        --output "${sampleid}_derep.fasta" \
        --sizeout \
        --threads 1

    # ---------------------------
    # Chimera screening
    # ---------------------------

    echo -e "\n#################################################################"
    echo -e "   Chimera screening for $sampleid"
    echo -e "#################################################################\n"

    if [[ "$chimera_early" -eq 1 && "$ampsize" -ge 200 ]]; then

        vsearch --uchime_denovo "${sampleid}_derep.fasta" \
            --abskew "$chimera_abskew" \
            --mindiv "$chimera_mindiv" \
            --nonchimeras "${sampleid}_nochim.fasta" \
            --threads 1

    else

        cp "${sampleid}_derep.fasta" "${sampleid}_nochim.fasta"

    fi

    # replace original fasta with dereplicated nonchimera fasta
    mv "${sampleid}_nochim.fasta" "$fasta"

    # cleanup derep intermediate
    rm -f "${sampleid}_derep.fasta"

    # ---------------------------
    # Count reads after chimera removal
    # ---------------------------

    awk '
    /^>/ {
        n=1

        # POSIX match() + RSTART/RLENGTH, so this works in mawk (the Docker image) as well
        # as gawk; ";size=" is 6 characters.
        if (match($0, /;size=[0-9]+/)) {
            n=substr($0, RSTART+6, RLENGTH-6)
        }

        sum+=n
    }

    END {
        print sum
    }
    ' "$fasta" > "${sampleid}_afterchimerareadcount.tmp"

    # Nothing may survive dereplication (reads under vsearch's 32 nt minimum) or chimera
    # screening. Clustering an empty file writes an empty .uc, which the check below would
    # report as a vsearch failure, so record that the sample has no OTUs instead.
    if ! grep -q '^>' "$fasta"; then
        echo "No OTUs remaining for $sampleid after dereplication and chimera screening"
        echo "no_otus" > "${sampleid}.status"
        return
    fi

    # ---------------------------
    # Primary OTU clustering
    # ---------------------------

    echo -e "\n#################################################################"
    echo -e "   Primary OTU clustering for $sampleid"
    echo -e "#################################################################\n"

    vsearch --cluster_fast "$fasta" \
        --id "$primaryclust" \
        --clusters "${sampleid}_OTU" \
        --uc "${sampleid}_primary.uc" \
        --iddef 3 \
        --threads 1

    # ---------------------------
    # Collect the cluster files once, by exact name (<sampleid>_OTU<number>), and check the
    # count against vsearch's own cluster count. Collecting before any renaming means a
    # renamed file can never be picked up a second time.
    # ---------------------------

    local f n otu_files=() kept=()
    for f in "${sampleid}"_OTU*; do
        n="${f#"${sampleid}_OTU"}"
        [[ -f "$f" && "$n" =~ ^[0-9]+$ ]] && otu_files+=("$f")
    done
    # vsearch writes a line per read to the .uc file, so a missing or empty one means it failed.
    if [[ ! -s "${sampleid}_primary.uc" ]]; then
        echo "error: primary clustering (vsearch) produced no output" > "${sampleid}.status"
        return
    fi
    n=$(awk '$1 == "C"' "${sampleid}_primary.uc" | wc -l)
    rm -f "${sampleid}_primary.uc"
    if [[ "${#otu_files[@]}" -ne "$n" ]]; then
        echo "error: primary clustering made $n clusters but ${#otu_files[@]} cluster files were found" > "${sampleid}.status"
        return
    fi

    # ---------------------------
    # Remove low-depth OTUs. --minreads counts READS: each record is a dereplicated sequence
    # carrying ;size=N, so the reads in a cluster are the sum of its sizes, not its records.
    # ---------------------------

    local read_count
    for f in "${otu_files[@]}"; do
        read_count=$(awk '/^>/ { n = 1; if (match($0, /;size=[0-9]+/)) n = substr($0, RSTART+6, RLENGTH-6); s += n } END { print s+0 }' "$f")
        if [[ "$read_count" -lt "$minreads" ]]; then
            rm -f "$f"
        else
            kept+=("$f")
        fi
    done

    # ---------------------------
    # Stop if no OTUs remain
    # ---------------------------

    if [[ "${#kept[@]}" -eq 0 ]]; then
        echo "No OTUs remaining for $sampleid after minreads filtering"
        echo "no_otus" > "${sampleid}.status"
        return
    fi

    # ---------------------------
    # Rename OTU read files (explicit list)
    # ---------------------------

    for f in "${kept[@]}"; do
        mv "$f" "${f}_Reads.fasta"
    done

    # ---------------------------
    # Generate primary OTU consensuses
    # ---------------------------

    Rscript "$scripts_directory/bip1_otu-consensus.R" "$sampleid" "$(pwd)" primary

    # ---------------------------
    # Secondary OTU clustering
    # ---------------------------

    echo -e "\n#################################################################"
    echo -e "   Secondary OTU clustering for $sampleid"
    echo -e "#################################################################\n"

    vsearch --cluster_fast "${sampleid}_OTUs.tmp" \
        --id "$finalclust" \
        --clusters "${sampleid}_FinalOTU" \
        --uc "${sampleid}_final.uc" \
        --iddef 3 \
        --threads 1

    rm -f "${sampleid}_OTUs.tmp"

    # ---------------------------
    # Rename secondary OTUs (explicit list, checked against vsearch's cluster count)
    # ---------------------------

    local final_files=()
    for f in "${sampleid}"_FinalOTU*; do
        n="${f#"${sampleid}_FinalOTU"}"
        [[ -f "$f" && "$n" =~ ^[0-9]+$ ]] && final_files+=("$f")
    done
    # vsearch writes a line per read to the .uc file, so a missing or empty one means it failed.
    if [[ ! -s "${sampleid}_final.uc" ]]; then
        echo "error: secondary clustering (vsearch) produced no output" > "${sampleid}.status"
        return
    fi
    n=$(awk '$1 == "C"' "${sampleid}_final.uc" | wc -l)
    rm -f "${sampleid}_final.uc"
    if [[ "${#final_files[@]}" -ne "$n" ]]; then
        echo "error: secondary clustering made $n clusters but ${#final_files[@]} cluster files were found" > "${sampleid}.status"
        return
    fi
    for f in "${final_files[@]}"; do
        mv "$f" "${f}_Reads.fasta"
    done

    # ---------------------------
    # Generate secondary OTU consensuses
    # ---------------------------

    Rscript "$scripts_directory/bip1_otu-consensus.R" "$sampleid" "$(pwd)" secondary

    # ---------------------------
    # Restore original reads into final OTUs
    # ---------------------------

    for f in "${sampleid}"_FinalOTU*_Reads.fasta; do

        [[ ! -f "$f" ]] && continue

        > "${f}.out"

        grep '^>' "$f" \
            | cut -d '>' -f2 \
            | cut -d '|' -f1,2 \
            | tr '|' '_' \
            | while IFS= read -r line; do

                cat "${line}_Reads.fasta" >> "${f}.out"

            done

        mv "${f}.out" "$f"

    done

    # ---------------------------
    # Cleanup primary OTUs
    # ---------------------------

    rm -f "${sampleid}"_OTU*_Reads.fasta

    # ---------------------------
    # Final renaming
    # ---------------------------

    rename 's/FinalOTU/OTU/g' "${sampleid}"_FinalOTU*_Reads.fasta
    rename 's/.fasta/.fas/g' "${sampleid}"_OTU*_Reads.fasta

    echo "ok" > "${sampleid}.status"
}

# Run one sample's OTU work in its own folder, otu_work/<sampleid>/, then move its outputs
# back. Many samples run at once; when they all worked in the shared primer folder, a
# directory listing taken while other jobs were creating, renaming and deleting files could
# miss this sample's files on some filesystems (seen on Docker Desktop bind mounts on macOS),
# silently dropping reads. A private folder is only ever touched by its own job. A sample that
# does not finish cleanly leaves <sampleid>.otu_failed, which stops the run after the step.
sampletask() {

    local fasta="$1"
    local sampleid="${fasta%.fasta}"
    local sdir="otu_work/${sampleid}"
    local status f

    if ! mkdir -p "$sdir" || ! mv -- "$fasta" "$sdir/"; then
        echo "could not set up $sdir" > "${sampleid}.otu_failed"
        return
    fi

    ( cd "$sdir" && sampletask_core "$fasta" )

    status=$(cat "$sdir/${sampleid}.status" 2>/dev/null)
    rm -f "$sdir/${sampleid}.status"
    for f in "$sdir"/*; do
        [[ -e "$f" ]] && mv -- "$f" .
    done
    rmdir "$sdir" 2>/dev/null

    case "$status" in
        ok|no_otus) ;;
        "")         echo "OTU step did not finish" > "${sampleid}.otu_failed" ;;
        *)          echo "$status" > "${sampleid}.otu_failed" ;;
    esac
}

process_primer_dir() {
    local name="$1"
    local fq="$2"
    local ampsize
    local minlen
    local maxlen
    local otu_primary
    local otu_secondary
    local tax_prob
    local bin_thresh
    local fwd_primer_seq
    local rev_primer_seq
    local startdir="$PWD"
    local marker
    local reference_library

    cd "$name" || exit 1

    IFS=$'\t' read ampsize minlen maxlen otu_primary otu_secondary \
         tax_prob fwd_primer_seq rev_primer_seq marker reflib bin_thresh <<< "$(

    awk -F'\t' -v OFS='\t' -v p="$name" \
        -v cfp="$COL_FWD_PRIMER" -v crp="$COL_REV_PRIMER" \
        -v c_amp="$COL_AMPLEN" -v c_min="$COL_MINLEN" -v c_max="$COL_MAXLEN" \
        -v c_prim="$COL_OTU_PRIMARY" -v c_sec="$COL_OTU_SECONDARY" -v c_tax="$COL_TAX_PROB" \
        -v c_fseq="$COL_FWD_PRIMER_SEQ" -v c_rseq="$COL_REV_PRIMER_SEQ" \
        -v c_marker="$COL_MARKER" -v c_reflib="$COL_REFLIB" -v c_bin="$COL_BIN_THRESH" '
    NR>1 {

        key = $cfp "_" $crp

        if (key == p) {

            # order must match the `read` below:
            # ampsize minlen maxlen otu_primary otu_secondary tax_prob
            # fwd_primer_seq rev_primer_seq marker reflib bin_thresh
            print \
                $c_amp, \
                $c_min, \
                $c_max, \
                $c_prim, \
                $c_sec, \
                $c_tax, \
                $c_fseq, \
                $c_rseq, \
                $c_marker, \
                $c_reflib, \
                $c_bin

            exit
        }
    }
    ' $wkdir/parameters.tsv
    )"

    primaryclust="$otu_primary"
    finalclust="$otu_secondary"

    echo "Processing $name..."

    # ---------------------------
    # Filter by length
    # ---------------------------

    T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

    local len_fq="${name}.len.fastq"
    filter_by_length "$name" "$fq" "$len_fq"

    T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
    log_step "${name}:Size filtering" "$TL" "$T2L" $(( T2 - T ))

    # ---------------------------
    # Detect UMI symmetry (based on primer set)
    # ---------------------------
    T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

    if diff \
        <(
            awk -F'\t' -v p="$name" -v cfp="$COL_FWD_PRIMER" -v crp="$COL_REV_PRIMER" -v c="$COL_FWD_UMI" '
            NR>1 {
                key = $cfp "_" $crp

                if (key == p) {
                    gsub(/\r/, "", $c)
                    gsub(/^[ \t]+|[ \t]+$/, "", $c)
                    print $c
                }
            }
            ' $wkdir/parameters.tsv
        ) \
        <(
            awk -F'\t' -v p="$name" -v cfp="$COL_FWD_PRIMER" -v crp="$COL_REV_PRIMER" -v c="$COL_REV_UMI" '
            NR>1 {
                key = $cfp "_" $crp

                if (key == p) {
                    gsub(/\r/, "", $c)
                    gsub(/^[ \t]+|[ \t]+$/, "", $c)
                    print $c
                }
            }
            ' $wkdir/parameters.tsv
        ) \
        >/dev/null; then

        umi_mode="symmetrical"

    else

        umi_mode="asymmetrical"

    fi

    echo "UMI mode: $umi_mode"

    # ---------------------------
    # Sanitize read headers
    # ---------------------------
    local clean_fq="${name}.clean.fastq"
    sanitize_headers "$len_fq" "$clean_fq"
    mv "$clean_fq" "$len_fq"

    # ---------------------------
    # Linked UMI demultiplexing
    # ---------------------------

    local umi_untrim="${name}.umi_untrim.fastq"
    local umi_demux="{name}.fastq"

    local linked_umis="${name}.linked_umis.adapters"
    awk -F'\t' \
        -v ov="$umi_overlap_min" -v err1="$error_umi1" -v err2="$error_umi2" \
        -v c_fu="$COL_FWD_UMI" -v c_ru="$COL_REV_UMI" \
        -v c_fp="$COL_FWD_PRIMER" -v c_rp="$COL_REV_PRIMER" -v p="$name" '
    function revcomp(seq,   i, c, out) {
        seq = toupper(seq); out = ""
        for (i = length(seq); i > 0; i--) {
            c = substr(seq, i, 1)
            if      (c=="A") out = out "T"
            else if (c=="T") out = out "A"
            else if (c=="C") out = out "G"
            else if (c=="G") out = out "C"
            else out = out "N"
        }
        return out
    }
    NR>1 && ($c_fp "_" $c_rp) == p && $c_fu!="" && $c_fu!="NA" && $c_ru!="" && $c_ru!="NA" {
        fwd = toupper($c_fu); rev = toupper($c_ru); rev_rc = revcomp(rev)
        fov = int(length(fwd) * ov); rov = int(length(rev) * ov)
        print ">" fwd "_" rev
        print fwd ";min_overlap=" fov ";max_error_rate=" err1 "..." rev_rc ";min_overlap=" rov ";max_error_rate=" err2
    }
    ' "$wkdir/parameters.tsv" > "$linked_umis"

    cutadapt -j "$cores" \
        -g file:"$linked_umis" \
        --action=trim \
        -o "$umi_demux" \
        --untrimmed-output "$umi_untrim" \
        "$len_fq"

    linked_count=0
    for f in *.fastq; do

        [[ "$f" == "${name}.fastq" ]] && continue
        [[ "$f" == "$len_fq" ]] && continue
        [[ "$f" == "$umi_untrim" ]] && continue

        c=$(seqkit stats -T "$f" | awk 'NR==2 {print $4}')
        linked_count=$((linked_count + c))
    done
    echo -e "${runid}\tLinkedUMI\t${name}\t${linked_count}" >> "$counts_file"

    # Almost no reads matching a UMI pair usually means the parameters file's UMIs don't
    # match the reads, most often Reverse UMIs entered in the wrong orientation.
    len_count=$(seqkit stats -T "$len_fq" | awk 'NR==2 {print $4}')
    if (( len_count > 0 && linked_count * 100 < len_count )); then
        umi_msg="$name: only $linked_count of $len_count reads ($(awk -v a="$linked_count" -v b="$len_count" 'BEGIN { printf "%.2f", 100 * a / b }')%) matched a UMI pair"
        echo "WARNING: $umi_msg"
        echo "  Check the UMIs in the parameters file. Reverse UMIs must be entered 5'->3' as on the"
        echo "  reverse primer oligo, i.e. reverse-complemented relative to how they appear at the end"
        echo "  of a read in the forward-primer orientation."
        umi_warnings+=("$umi_msg")
    fi


    # ---------------------------
    # Recover single-sided UMIs (symmetrical UMIs only)
    # ---------------------------

    if [[ "$umi_mode" == "symmetrical" ]]; then

        local fwd_single="${name}.fwd_single.fastq"
        local rev_single="${name}.rev_single.fastq"
        local single_all="${name}.single_all.fastq"
        local dup_ids="${name}.dup_read_ids.txt"
        local single_unique="${name}.single_unique.fastq"

        # ---------------------------
        # Find forward-only UMIs
        # ---------------------------

        # Keep the original read ID: it is what lets the tag-switch step below tell
        # "the same read matched in both passes" apart from "two different reads that
        # happen to have the same sequence" (possible with short, low-error reads).
        # The adapter name moves into the comment, where the merge step's
        # `seqkit grep -n -r` still finds it.
        cutadapt -j "$cores" \
            -g file:"$wkdir/fwd_umis.fasta" \
            --action=none \
            --rename='{id} {adapter_name}_fwd' \
            -o "$fwd_single" \
            --untrimmed-output /dev/null \
            "$umi_untrim"

        # ---------------------------
        # Find reverse-only UMIs
        # ---------------------------

        cutadapt -j "$cores" \
            -a file:"$wkdir/rev_umis_rc.fasta" \
            --action=none \
            --rename='{id} {adapter_name}_rev' \
            -o "$rev_single" \
            --untrimmed-output /dev/null \
            "$umi_untrim"

        # ---------------------------
        # Combine singleton matches
        # ---------------------------

        cat "$fwd_single" "$rev_single" > "$single_all"

        # ---------------------------
        # Detect duplicated reads (tag-switched UMIs)
        # ---------------------------

        # A tag-switched read carries a forward UMI from one sample and a reverse UMI
        # from another, so it matches in BOTH passes and lands in single_all twice -
        # once from each pass, under the SAME original read ID. The linked pass has
        # already claimed every correctly paired read, so any leftover that matched
        # both ends is tag-switched, and both copies must go.
        #
        # Keyed on the read ID rather than the sequence: two distinct reads from the
        # same sample can legitimately share a sequence (short, high-depth, low-error
        # data), and deduplicating on sequence would discard both. The ID identifies
        # the physical read.
        seqkit seq -n -i "$single_all" \
            | sort \
            | uniq -d \
            > "$dup_ids"

        # ---------------------------
        # Remove duplicated (tag-switched) reads
        # ---------------------------

        # seqkit grep matches on the sequence ID by default (not the full header),
        # which is exactly the original read ID preserved by --rename above.
        seqkit grep -v -f "$dup_ids" "$single_all" > "$single_unique"

        # ---------------------------
        # Log counts
        # ---------------------------

        local fwd_only_count=0
        local rev_only_count=0
        local rescued_count=0

        [[ -s "$fwd_single" ]] && \
            fwd_only_count=$(seqkit stats -T "$fwd_single" | awk 'NR==2 {print $4}')

        [[ -s "$rev_single" ]] && \
            rev_only_count=$(seqkit stats -T "$rev_single" | awk 'NR==2 {print $4}')

        [[ -s "$single_unique" ]] && \
            rescued_count=$(seqkit stats -T "$single_unique" | awk 'NR==2 {print $4}')

        [[ ! -f "$rescue_log" ]] && echo -e "runid\tstage\tprimer_pair\treads" > "$rescue_log"
        echo -e "${runid}\tForwardOnlyUMI\t${name}\t${fwd_only_count}" >> "$rescue_log"
        echo -e "${runid}\tReverseOnlyUMI\t${name}\t${rev_only_count}" >> "$rescue_log"
        echo -e "${runid}\tSingletonUMIRescue\t${name}\t${rescued_count}" >> "$rescue_log"

        # ---------------------------
        # Merge linked and singleton reads into sample FASTQs
        # ---------------------------

        awk -F'\t' -v p="$name" -v cfp="$COL_FWD_PRIMER" -v crp="$COL_REV_PRIMER" \
            -v c_samp="$COL_SAMPLE" -v c_fu="$COL_FWD_UMI" -v c_ru="$COL_REV_UMI" '
        NR>1 {
            key = $cfp "_" $crp

            if (key == p) {
                print $c_samp "\t" $c_fu "\t" $c_ru
            }
        }
        ' $wkdir/parameters.tsv |

        while IFS=$'\t' read -r sample fwd_umi rev_umi; do

            fwd_umi=$(echo "$fwd_umi" | tr '[:lower:]' '[:upper:]')
            rev_umi=$(echo "$rev_umi" | tr '[:lower:]' '[:upper:]')

            local linked_file="${fwd_umi}_${rev_umi}.fastq"
            local sample_out="${sample}__${name}.fastq"

            # singleton rescue patterns
            local fwd_pattern="${fwd_umi}_fwd"
            local rev_pattern="${rev_umi}_rev"

            # start fresh
            > "$sample_out"

            # ---------------------------
            # linked 1-1 reads
            # ---------------------------

            if [[ ! -f "$linked_file" ]]; then
                echo "WARNING: missing linked UMI file for $sample ($linked_file)" >&2
            fi

            [[ -f "$linked_file" ]] && \
                cat "$linked_file" >> "$sample_out"

            # ---------------------------
            # rescued forward-only reads
            # ---------------------------

            if [[ -f "$single_unique" ]]; then

                seqkit grep \
                    -n -r \
                    -p "$fwd_pattern" \
                    "$single_unique" \
                    >> "$sample_out"
            fi

            # ---------------------------
            # rescued reverse-only reads
            # ---------------------------

            if [[ -f "$single_unique" ]]; then

                seqkit grep \
                    -n -r \
                    -p "$rev_pattern" \
                    "$single_unique" \
                    >> "$sample_out"
            fi

            # remove empty outputs
            [[ ! -s "$sample_out" ]] && rm -f "$sample_out"

        done

        rm -f "$fwd_single"
        rm -f "$rev_single"
        rm -f "$single_all"
        rm -f "$dup_ids"
        rm -f "$single_unique"
    fi


    # ---------------------------
    # Rename UMI-linked outputs (asymmetrical only)
    # ---------------------------

    if [[ "$umi_mode" == "asymmetrical" ]]; then
        # Restrict to THIS primer pair's rows. parameters.tsv holds every primer
        # pair, and new_file is keyed on $name (the current primer pair), not the
        # row's. A sample run under two markers therefore renames to the same
        # target twice, and the second mv overwrites the correctly demultiplexed
        # file with the near-empty file for the other marker's reverse UMI (which
        # only caught stray reads in this subdirectory). The symmetrical merge loop
        # already filters this way.
        awk -F'\t' -v p="$name" \
            -v c_samp="$COL_SAMPLE" -v c_fu="$COL_FWD_UMI" -v c_ru="$COL_REV_UMI" \
            -v cfp="$COL_FWD_PRIMER" -v crp="$COL_REV_PRIMER" '
        NR>1 && ($cfp "_" $crp) == p && $c_fu!="" && $c_ru!="" && $c_fu!="NA" && $c_ru!="NA" {

            fwd = toupper($c_fu)
            rev = toupper($c_ru)

            sample = $c_samp
            primer = $cfp "_" $crp

            print fwd "_" rev "\t" sample "\t" primer
        }
        ' $wkdir/parameters.tsv |
        while IFS=$'\t' read -r umi_pair sample primer; do

            old_file="${umi_pair}.fastq"
            new_file="${sample}__${name}.fastq"

            if [[ -f "$old_file" ]]; then
                mv "$old_file" "$new_file"
            fi

        done
    fi

    # ---------------------------
    # Final cleanup (shared for sym + asym)
    # ---------------------------

    rm -f "${name}.umi_untrim.fastq"
    rm -f "${name}.len.fastq"
    rm -f "${name}.fastq"
    rm -f "$linked_umis"

    # ---------------------------
    # Archive raw fastq files
    # ---------------------------
    
    cp *.fastq $wkdir/Individual_Raw_Fastq_Files/

    # ---------------------------
    # Convert to FASTA
    # ---------------------------

    for f in *.fastq; do
        [[ -e "$f" ]] || continue

        base="${f%.fastq}"
        seqkit fq2fa "$f" -o "${base}.fasta"
    done
    rm *.fastq

    T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
    log_step "${name}:Demultiplexing" "$TL" "$T2L" $(( T2 - T ))

    # ---------------------------
    # Process each FASTA file
    # ---------------------------

    T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

    # Samples run in parallel, so each writes to its own log and only this loop prints to the
    # terminal: one line per sample as it finishes. The logs are combined into one file after.
    local fastas=() sid pid otus result secs i=0 done_n=0
    local -A job_sample=() job_start=()
    for fasta in *.fasta; do
        [ -e "$fasta" ] && fastas+=("$fasta")
    done
    local n=${#fastas[@]}
    mkdir -p otu_logs
    echo -e "\n******** OTU clustering ${n} samples for ${name}, up to ${cores} at a time"

    while (( i < n || ${#job_sample[@]} )); do
        if (( i < n && ${#job_sample[@]} < cores )); then
            sid="${fastas[i]%.fasta}"
            sampletask "${fastas[i]}" > "otu_logs/${sid}.log" 2>&1 &
            job_sample[$!]="$sid"
            job_start[$!]=$SECONDS
            (( i++ ))
            continue
        fi

        pid=
        wait -n -p pid
        [[ -n "$pid" && -n "${job_sample[$pid]:-}" ]] || continue
        sid="${job_sample[$pid]}"
        secs=$(( SECONDS - job_start[$pid] ))
        unset "job_sample[$pid]" "job_start[$pid]"
        (( done_n++ ))

        if [[ -f "${sid}.otu_failed" ]]; then
            result="FAILED ($(cat "${sid}.otu_failed"))"
        else
            otus=( "${sid}"_OTU*_Reads.fas )
            [[ -e "${otus[0]}" ]] && result="${#otus[@]} OTU(s)" || result="no OTUs"
        fi
        printf '  [%*d/%d] %-30s %-10s %4ss\n' "${#n}" "$done_n" "$n" "${sid%__"$name"}" "$result" "$secs"
        if [[ -f "${sid}.otu_failed" ]]; then
            tail -n 15 "otu_logs/${sid}.log" | sed 's/^/        | /'
        fi
    done

    for fasta in "${fastas[@]}"; do
        sid="${fasta%.fasta}"
        echo "=================== ${sid} ==================="
        cat "otu_logs/${sid}.log" 2>/dev/null
        echo
    done > "${name}_OTU_clustering.log"
    rm -rf otu_logs
    echo "  Full per-sample output: ${runid}_results/${name}/${name}_OTU_clustering.log"

    # Every sample must have finished its OTU step cleanly. Stop rather than continue with
    # samples silently missing from the results.
    failed=(*.otu_failed)
    leftover=()
    for f in *_OTU* *_FinalOTU*; do
        [[ "$f" =~ _(Final)?OTU[0-9]+$ ]] && leftover+=("$f")
    done
    if (( ${#failed[@]} )) || (( ${#leftover[@]} )) || [[ -n "$(ls -A otu_work 2>/dev/null)" ]]; then
        echo "ERROR: the OTU step did not complete for every sample (primer pair $name):"
        for f in "${failed[@]}"; do echo "  ${f%.otu_failed}: $(cat "$f")"; done
        (( ${#leftover[@]} )) && echo "  unprocessed cluster files: ${leftover[*]}"
        [[ -n "$(ls -A otu_work 2>/dev/null)" ]] && echo "  sample work folders left in otu_work/: $(ls otu_work | tr '\n' ' ')"
        echo "  Results would be incomplete, so BIP has stopped. Re-run, and if it recurs, use fewer"
        echo "  cores (--cores_to_leave) or a working directory on a native filesystem."
        exit 1
    fi
    rmdir otu_work 2>/dev/null

    T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
    log_step "${name}:OTU clustering" "$TL" "$T2L" $(( T2 - T ))

    # ---------------------------
    # Sum after chimera read counts and remove temporary files
    # ---------------------------

    # Collect into an array and guard so an empty marker reports 0 instead.
    chimera_tmps=(*_afterchimerareadcount.tmp)
    if (( ${#chimera_tmps[@]} )); then
        chimerareadcount=$(awk '{s+=$1} END{print s+0}' "${chimera_tmps[@]}")
        rm -f "${chimera_tmps[@]}"
    else
        chimerareadcount=0
    fi
    echo -e "${runid}\tAfterChimRemovalOrDerep\t${name}\t${chimerareadcount}" >> "$counts_file"

    # delete any OTU component read or OTU consensus read files that are empty
    find ./ -name "*.fas" -size 0 -delete
    find ./ -name "*.tmp" -size 0 -delete

    # delete pre-clustered read files
    find ./ -name "*.fasta" -delete

    # move OTU component read files into their own folder
    mkdir -m 777 OTU_Component_Reads
    rename 's/Reads.fas/Reads.fasta/g' *Reads.fas
    for f in *Reads.fasta; do
        [[ -e "$f" ]] || continue
        mv -- "$f" OTU_Component_Reads/
    done

    # rename OTU consensus files
    rename 's/tmp/fasta/g' *finalOTUs.tmp
    rename 's/finalOTUs/OTUs/g' *finalOTUs.fasta

    # merge all OTUs into a single file, including 'empty barcode' samples
    echo -e "******** Merging all OTU consensus sequences into a single FASTA file"
    otu_files=(*OTUs.fasta)
    if (( ${#otu_files[@]} == 0 )); then
        echo "WARNING: no OTUs remaining for ${name}; skipping downstream steps for this marker"
        echo -e "${runid}\tAfterOTUClustering\t${name}\t0" >> "$counts_file"
        echo -e "${runid}\tReadsInOTUs\t${name}\t0" >> "$counts_file"
        echo -e "${runid}\tReadsInTargetBarcodes\t${name}\t0" >> "$counts_file"
        echo "Finished $name"
        cd "$startdir" || exit 1
        return 0
    fi
    cat "${otu_files[@]}" > "$runid".fasta
    rm -f "${otu_files[@]}"

    otucount=$(
        awk '
        /^>/ {
            if (match($0, /reads-[0-9]+/)) {
                s = substr($0, RSTART, RLENGTH)
                sub(/reads-/, "", s)
                sum += s
            }
        }
        END {
            print sum + 0
        }
        ' "${runid}.fasta"
    )

    echo -e "${runid}\tAfterOTUClustering\t${name}\t${otucount}" >> "$counts_file"

    # ---------------------------
    # Indel correction, NUMT detection (only if COI-5P and full-length barcode)
    # ---------------------------

    if [[ "$marker" == "COI-5P" ]] && (( ampsize >= coi_min_expected && ampsize <= coi_max_expected )); then
        T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

        echo -e "******** Auto-trim and indel correction"
        if [[ -n "$ref_seq_corr" ]]; then
            if [[ ! -f "$ref_seq_corr" ]]; then
                echo "ERROR: --ref_seq_corr file not found: $ref_seq_corr"
                exit 1
            fi
            errcorr_ref="$ref_seq_corr"
        else
            errcorr_refs=("$reference_library_directory"/reference_seqs_*.fasta)
            if [[ ${#errcorr_refs[@]} -ne 1 || ! -f "${errcorr_refs[0]}" ]]; then
                echo "ERROR: expected exactly 1 reference_seqs_*.fasta in $reference_library_directory (found ${#errcorr_refs[@]})"
                exit 1
            fi
            errcorr_ref="${errcorr_refs[0]}"
        fi
        python3.11 "$scripts_directory/bip2_homopolymer_error_fix.py" "$runid".fasta "$errcorr_ref" "$coi_expected_lengths"

        # combine trimmed/corrected and un-trimmed/un-corrected sequences into a single FASTA file
        awk '/^>/{if (seq) print seq; print; seq=""; next} {seq=seq $0} END{print seq}' problem_seqs.fasta > problem_seqs_singleline.fasta
        (cat edited_seqs.fasta; echo; cat problem_seqs_singleline.fasta) > "$runid".corrected.fasta
        rm "$runid".fasta edited_seqs.fasta problem_seqs.fasta problem_seqs_singleline.fasta

        T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
        log_step "${name}:Indel correction" "$TL" "$T2L" $(( T2 - T ))
    else
        mv "$runid".fasta "$runid".corrected.fasta
        rm "$runid".fasta
    fi

    # Chimera screening of OTU consensus sequences (non-PE / long-read mode)
    if [[ "$chimera_early" -eq 0 && "$ampsize" -ge 200 ]]; then
        T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")
        echo -e "******** Chimera screening OTU consensus sequences"
        vsearch --uchime_denovo "$runid".corrected.fasta \
            --abskew "$chimera_abskew" \
            --mindiv "$chimera_mindiv" \
            --nonchimeras "$runid".corrected.nochim.fasta \
            --threads $cores
        mv "$runid".corrected.nochim.fasta "$runid".corrected.fasta
        T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
        log_step "${name}:Chimera removal" "$TL" "$T2L" $(( T2 - T ))
    fi

    # assign taxonomy to OTUs and output summary table
    T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

    echo -e "******** Assigning taxonomy to OTU consensus sequences"
    # Exclude the search databases: BIP writes a .udb alongside the reference
    # FASTA (see BIN matching below), so without this the second run of any
    # library would match both files and abort with "Multiple files match".
    reflib_matches=$(find "$reference_library_directory" -maxdepth 1 -type f \
        \( -name "${reflib}*.fasta" -o -name "${reflib}*.fa" -o -name "${reflib}*.fna" \) 2>/dev/null)
    reflib_count=$(echo "$reflib_matches" | grep -c .)
    if [[ "$reflib_count" -eq 0 ]]; then
        echo "ERROR: No reference library file matching ${reflib}* found in $reference_library_directory"
        echo "Files available in $reference_library_directory:"
        find "$reference_library_directory" -maxdepth 1 -type f | sort | sed 's/^/  /'
        exit 1
    fi
    if [[ "$reflib_count" -gt 1 ]]; then
        echo "ERROR: Multiple files match ${reflib}* in $reference_library_directory -- use a more specific name:"
        echo "$reflib_matches" | sort | sed 's/^/  /'
        exit 1
    fi
    reflib=$(basename "$reflib_matches")
    reflib_path="$reference_library_directory/$reflib"
    if [[ "$(head -c1 "$reflib_path" 2>/dev/null)" != ">" ]]; then
        echo "ERROR: Reference library does not appear to be a FASTA file: $reflib_path"
        exit 1
    fi
    first_header=$(grep -m1 '^>' "$reflib_path")
    if [[ "$first_header" != *";tax="* ]]; then
        echo "ERROR: Reference library is not in SINTAX format (no ';tax=' in first header): $reflib_path"
        echo "  First header: $first_header"
        exit 1
    fi
    if [[ "$marker" == "COI-5P" ]] && (( ampsize >= 300 )); then
        id_part="${first_header%%;tax=*}"
        if [[ "$id_part" != *"|"* ]]; then
            echo "ERROR: BIN matching requires sequence IDs in 'ID|BIN' format, but no '|' found before ';tax=' in first header: $reflib_path"
            echo "  First header: $first_header"
            exit 1
        fi
    fi
    vsearch --sintax "$runid".corrected.fasta \
        --db "$reference_library_directory/$reflib" \
        --tabbedout "$runid".table \
        --strand plus \
        --sintax_cutoff $tax_prob \
        "${sintax_opts[@]}"

    # reformat tax id table
    echo -e "******** Reformatting SINTAX ID table"
    awk 'BEGIN{OFS="\t"; print "SeqName","TaxAssign","Strand","TaxAssignFinal";}{print}' "$runid".table > tmpfile && mv tmpfile "$runid".table

    T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
    log_step "${name}:Tax assignment" "$TL" "$T2L" $(( T2 - T ))

    # assign OTUs as target or NTS (non-target sequence) 
    T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

    echo -e "******** Inferring barcodes"
    inputreadcount=$(awk 'NR==2 {print $4}' $counts_file)
    Rscript $scripts_directory/bip3_barcode_inference.R "$runid" "$(pwd)" "$(echo "$inputreadcount * $minreadssatellite" | bc -l)" "$marker" "$name" "$min_rank_match" "$coi_expected_lengths" "$minreads"

    T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
    log_step "${name}:Barcode inference" "$TL" "$T2L" $(( T2 - T ))

    # ---------------------------
    # BIN matching (COI-5P, amplicon >= 300 bp)
    # ---------------------------
    if [[ "$marker" == "COI-5P" ]] && (( ampsize >= 300 )); then
        T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")
        echo -e "******** Performing BIN matching"

        sintax_fasta="$reference_library_directory/$reflib"
        vsearch_udb="${sintax_fasta%.fasta}.udb"

        if [[ ! -f "$vsearch_udb" ]]; then
            echo -e "******** Creating vsearch UDB from reference library..."
            sed '/^>/ s/ /_/g' "$sintax_fasta" > sintax_nospaces.fasta
            if ! vsearch --makeudb_usearch sintax_nospaces.fasta --output "$vsearch_udb" --threads $cores; then
                # Don't leave a partial UDB behind, or every later run would try to use it.
                rm -f "$vsearch_udb"
            fi
            rm sintax_nospaces.fasta
        fi

        min_overlap_bp=$(printf "%.0f" "$(echo "$ampsize * $minoverlap" | bc -l)")

        if [[ ! -f "$vsearch_udb" ]] || ! vsearch --usearch_global "$runid".corrected.fasta \
            --db "$vsearch_udb" \
            --blast6out bin_raw.txt \
            --id "$bin_thresh" \
            --maxhits "$bin_maxhits" \
            --maxaccepts "$bin_maxaccepts" \
            --threads $cores; then
            echo "ERROR: BIN matching failed for $name (vsearch could not build or search $vsearch_udb)."
            echo "  The run will continue with BIN set to NA for this primer pair."
            bin_failures+=("$name: vsearch could not build or search $vsearch_udb")
            # An empty hit list gives every OTU BIN = NA below, keeping the output format.
            : > bin_raw.txt
        fi

        awk -F'\t' -v mo="$min_overlap_bp" '
        BEGIN { OFS="\t" }
        $4 >= mo {
            n = split($2, h, "|")
            bin = (n >= 2) ? h[2] : "NA"
            sub(/;.*/, "", bin)
            print $1, bin, $3, $4
        }' bin_raw.txt | sort -k1,1 -k3nr -k4nr | awk '!seen[$1]++' > bin_results.tsv
        rm bin_raw.txt

        python3 - "${runid}__${name}__OTUDetails.tsv" bin_results.tsv <<'PYEOF'
import csv, re, sys
otu_file, bin_file = sys.argv[1], sys.argv[2]
# bip3 rewrites "reads-N" in an OTU's name after negative-control subtraction, so match
# on the name without it.
key = lambda name: re.sub(r'\|reads-\d+', '', name)
bin_map = {}
with open(bin_file) as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 4:
            bin_map[key(p[0])] = (p[1], p[2], p[3])
rows = []
with open(otu_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = reader.fieldnames + ["BIN", "Pct_Match", "Overlap_bp"]
    for row in reader:
        hit = bin_map.get(key(row["Sequence_Name"]))
        if hit:
            row["BIN"], row["Pct_Match"], row["Overlap_bp"] = hit
        else:
            row["BIN"] = row["Pct_Match"] = row["Overlap_bp"] = "NA"
        rows.append(row)
with open(otu_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
PYEOF

        rm bin_results.tsv
        T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
        log_step "${name}:BIN matching" "$TL" "$T2L" $(( T2 - T ))
    fi

    rm "$runid".table "$runid".corrected.fasta

    # count reads in target and NTS OTUs
    otureadcount=$(awk -F'\t' 'NR > 1 { sum += $4 } END {print sum}' "${runid}__${name}__OTUDetails.tsv")
    echo -e "${runid}\tReadsInOTUs\t${name}\t${otureadcount}" >> "$counts_file"

    # count reads in target OTUs
    targetreadcount=$(awk -F'\t' 'NR > 1 && $5 == "target" { sum += $4 } END { print sum }' "${runid}__${name}__OTUDetails.tsv")
    echo -e "${runid}\tReadsInTargetBarcodes\t${name}\t${targetreadcount}" >> "$counts_file"

    # clean up primer directory
    if [[ "$marker" == "COI-5P" ]] && (( ampsize >= coi_min_expected && ampsize <= coi_max_expected )); then

        mkdir -m 777 Error_Correction_Files
        mv vsearch_hit.tsv Error_Correction_Ref_Seqs.tsv
        mv uncertain_edits.fasta Error_Correction_Uncertain_Edits.fasta
        mv alignment_dir Error_Correction_Alignments
        mv Error_Correction_Ref_Seqs.tsv Error_Correction_Uncertain_Edits.fasta Error_Correction_Files
        mv ./Error_Correction_Alignments Error_Correction_Files
    fi

    # compress OTU component reads
    for f in ./OTU_Component_Reads/*Reads.fasta; do
        gzip $f
    done

    # end primer set
    echo "Finished $name"
    cd "$startdir" || exit 1
}

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
########### ACTUAL SCRIPT #####################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################

# ---------------------------
# Check working directory for presence of at least one .fastq.gz file and parameters Excel file
# ---------------------------
T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

shopt -s nullglob

cd "$wkdir" || { echo "\n\nERROR: Cannot access working directory $wkdir"; read -p "Press Enter to exit..."; exit 1; }

# Inputs come from --fastq / --params when given, otherwise from the working directory.
# Paths are absolute from here on. The read files are sorted by name, so for paired-end
# data R1 comes before R2 whatever mix of .fastq and .fastq.gz they are.
fastq_dir="${fastq_dir:-$PWD}"
fastq_files=()
while IFS= read -r f; do
    [[ -n "$f" ]] && fastq_files+=("$f")
done < <(printf '%s\n' "$fastq_dir"/*.fastq.gz "$fastq_dir"/*.fastq | LC_ALL=C sort)
# Collected into its own array: `param_file` must end up a plain scalar, because
# bash cannot export array variables and the python steps below read it from the
# environment.
if [[ -f "$param_path" ]]; then
    param_files=("$param_path")
else
    param_files=("${param_path:-$PWD}"/*.xlsx)
fi

if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo -e "\n\nERROR: Missing .fastq or .fastq.gz files in $fastq_dir"
    read -p "Press Enter to exit..."
    exit 1
fi

if [[ ${#param_files[@]} -eq 0 ]]; then
    echo "ERROR: No .xlsx parameter file found in ${param_path:-$PWD}"
    read -p "Press Enter to exit..."
    exit 1
fi

if [[ ${#param_files[@]} -gt 1 ]]; then
    echo "ERROR: Multiple .xlsx files found in ${param_path:-$PWD}. Use only one, or point --params at the file."
    read -p "Press Enter to exit..."
    exit 1
fi

param_file="${param_files[0]}"

# ---------------------------
# Convert paramaters file from Excel to TSV
# ---------------------------

tsvfile="parameters.tsv"
export param_file tsvfile

python3 - <<'EOF'
import os, sys
import pandas as pd

# Sheet and column names come from the COLUMN NAMES block at the top of this
# script; the headings are written through as-is.
param_file = os.environ["param_file"]
tsvfile    = os.environ["tsvfile"]
sheet      = os.environ["SHEET_UMIS_PRIMERS"]
neg_name   = os.environ["COLNAME_NEG_CTRL"]
samp_name  = os.environ["COLNAME_SAMPLE"]

df = pd.read_excel(param_file, sheet_name=sheet)
df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
df.columns = [str(c).strip() for c in df.columns]
# Stray spaces typed in Excel (e.g. "Diptera " or " ACGT") would otherwise break exact
# matching downstream, so trim leading/trailing whitespace from every text cell.
df = df.apply(lambda col: col.map(lambda v: v.strip() if isinstance(v, str) else v))

# The control column is found by how its header starts, so templates can add a hint such as
# "Negative Control (yes/YES/y/1)". It is renamed to the standard name, which is what
# parameters.tsv carries and what bip3 and bip4 look up.
neg_matches = [c for c in df.columns if c.lower().startswith(neg_name.lower())]
if len(neg_matches) > 1:
    sys.exit(f"ERROR: more than one column in sheet '{sheet}' starts with '{neg_name}': {neg_matches}. Keep only one.")
if neg_matches and neg_matches[0] != neg_name:
    df = df.rename(columns={neg_matches[0]: neg_name})

for required in (samp_name, neg_name):
    if required not in df.columns:
        sys.exit(f"ERROR: column '{required}' not found in sheet '{sheet}'. "
                 f"Columns present: {list(df.columns)}")

# ---------------------------------------------------------------------------
# Sample names must contain only [A-Za-z0-9._-]. 
# ---------------------------------------------------------------------------
present  = df[samp_name].notna()
original = df.loc[present, samp_name].astype(str)
cleaned  = original.str.strip().str.replace(r"[^A-Za-z0-9._-]+", "_", regex=True)

renamed = sorted({(o, c) for o, c in zip(original, cleaned) if o != c})
if renamed:
    print(f"WARNING: {len(renamed)} sample name(s) contained characters outside [A-Za-z0-9._-]; replaced with '_'.")
    for o, c in renamed[:10]:
        print(f"           {o!r} -> {c!r}")
    if len(renamed) > 10:
        print(f"           ... and {len(renamed) - 10} more")

# Two DIFFERENT names must never collapse onto the same sanitized name, or their
# reads would be silently merged into one sample. (The same name repeating across
# rows is normal - a sample occupies several wells - and is not a collision.)
sources = {}
for o, c in zip(original, cleaned):
    sources.setdefault(c, set()).add(o)
clashes = {c: sorted(v) for c, v in sources.items() if len(v) > 1}
if clashes:
    detail = "; ".join(f"{v} -> {c!r}" for c, v in sorted(clashes.items()))
    sys.exit(f"ERROR: distinct sample names collide once illegal characters are replaced with '_': {detail}. "
             f"Rename them in the '{sheet}' sheet so they stay distinct.")

df.loc[present, samp_name] = cleaned
sys.stdout.flush()

# A well is a negative control only if this column says yes, y or 1 (any case; Excel's
# 1.0 counts as 1). Anything else, blank included, means it is not a control. The column is
# rewritten as "yes" or blank so every later step (bip3, bip4) reads one consistent value.
neg_col = df[neg_name].fillna('').astype(str).str.strip().str.lower()
neg_col = neg_col.str.replace(r'^1\.0+$', '1', regex=True)
valid_yes = {'yes', 'y', '1'}
df[neg_name] = neg_col.where(neg_col.isin(valid_yes), '').map(lambda v: 'yes' if v else '')

# Negative controls should carry no taxonomy. Users sometimes fill a column down the whole
# sheet and leave identifications in control rows, so clear them here.
tax_cols = [c for c in ("Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species")
            if c in df.columns]
is_ctrl  = neg_col.isin(valid_yes)
had_tax  = is_ctrl & df[tax_cols].apply(lambda col: col.notna() & (col.astype(str).str.strip() != "")).any(axis=1)
if had_tax.any():
    df.loc[is_ctrl, tax_cols] = pd.NA
    print(f"NOTE: cleared taxonomy from {int(had_tax.sum())} negative-control row(s) "
          f"(columns: {', '.join(tax_cols)}).")

df.to_csv(tsvfile, sep="\t", index=False)
EOF

if [[ $? -ne 0 || ! -s "$tsvfile" ]]; then
    echo "ERROR: could not build $tsvfile from $param_file (see message above)."
    read -p "Press Enter to exit..."
    exit 1
fi

# ---------------------------
# Sanitize and validate the dictionary file
# ---------------------------

dictfile="$wkdir/dictionary_BIP.tsv"

# ---------------------------
# Merge Dictionary tab from Excel into local dictionary file
# ---------------------------

export dictfile

python3 - <<'EOF'
import os, sys
import pandas as pd

dictfile   = os.environ["dictfile"]
param_file = os.environ["param_file"]
sheet      = os.environ["SHEET_DICT_UPDATE"]

# Column names come from the COLUMN NAMES block at the top of this script.
dict_cols = [os.environ[k] for k in (
    'DICTCOL_FWD_PRIMER_NAME', 'DICTCOL_REV_PRIMER_NAME',
    'DICTCOL_FWD_PRIMER_SEQ',  'DICTCOL_REV_PRIMER_SEQ',
    'DICTCOL_MARKER',          'DICTCOL_REFLIB',
    'DICTCOL_MINLEN',          'DICTCOL_MAXLEN',    'DICTCOL_AMPLEN',
    'DICTCOL_OTU_PRIMARY',     'DICTCOL_OTU_SECONDARY',
    'DICTCOL_TAX_PROB',        'DICTCOL_BIN_THRESH',
)]
fwd_name, rev_name = dict_cols[0], dict_cols[1]

# A misnamed tab used to be swallowed here, leaving the dictionary silently
# empty and failing much later with a "primer pair not found". 
xl = pd.ExcelFile(param_file)
if sheet in xl.sheet_names:
    df_new = pd.read_excel(param_file, sheet_name=sheet)
    df_new = df_new.loc[:, ~df_new.columns.astype(str).str.startswith("Unnamed")]
    df_new.columns = [str(c).strip() for c in df_new.columns]

    missing = [c for c in dict_cols if c not in df_new.columns]
    if missing:
        sys.exit(f"ERROR: sheet '{sheet}' is missing column(s): {missing}. "
                 f"Columns present: {list(df_new.columns)}")

    df_new = df_new[dict_cols]
    # Trim stray spaces from every text cell, e.g. a Reference Library name typed with a
    # trailing space would otherwise fail the library lookup.
    df_new = df_new.apply(lambda col: col.map(lambda v: v.strip() if isinstance(v, str) else v))
    df_new = df_new.dropna(subset=[fwd_name, rev_name])
    df_new = df_new[df_new[fwd_name].astype(str).str.strip() != '']
else:
    print(f"NOTE: no '{sheet}' tab in {os.path.basename(param_file)}; "
          f"using the existing dictionary unchanged.")
    df_new = pd.DataFrame(columns=dict_cols)

if os.path.exists(dictfile):
    df_dict = pd.read_csv(dictfile, sep='\t')
    # A dictionary saved by an earlier run may still hold untrimmed values.
    df_dict = df_dict.apply(lambda col: col.map(lambda v: v.strip() if isinstance(v, str) else v))
else:
    df_dict = pd.DataFrame(columns=dict_cols)

for _, row in df_new.iterrows():
    fwd = str(row[fwd_name]).strip()
    rev = str(row[rev_name]).strip()
    mask = (
        df_dict[fwd_name].astype(str).str.strip() == fwd
    ) & (
        df_dict[rev_name].astype(str).str.strip() == rev
    )
    df_dict = df_dict[~mask]
    df_dict = pd.concat([df_dict, pd.DataFrame([row])], ignore_index=True)

df_dict.to_csv(dictfile, sep='\t', index=False)
EOF

[[ ! -f "$dictfile" ]] && { echo "ERROR: dictionary.tsv not found and Dictionary tab in $param_file was empty. Cannot continue."; read -p "Press Enter to exit..."; exit 1; }

# Resolve the dictionary's columns by name, so its column order is free to change
DCOL_FWD_SEQ=$(resolve_col      "$DICTCOL_FWD_PRIMER_SEQ"  "$dictfile") || exit 1
DCOL_REV_SEQ=$(resolve_col      "$DICTCOL_REV_PRIMER_SEQ"  "$dictfile") || exit 1
DCOL_MARKER=$(resolve_col       "$DICTCOL_MARKER"          "$dictfile") || exit 1
DCOL_OTU_PRIMARY=$(resolve_col  "$DICTCOL_OTU_PRIMARY"     "$dictfile") || exit 1
DCOL_OTU_SECONDARY=$(resolve_col "$DICTCOL_OTU_SECONDARY"  "$dictfile") || exit 1
DCOL_TAX_PROB=$(resolve_col     "$DICTCOL_TAX_PROB"        "$dictfile") || exit 1
DCOL_BIN_THRESH=$(resolve_col   "$DICTCOL_BIN_THRESH"      "$dictfile") || exit 1

sanitized=$(awk -F'\t' -v OFS='\t' -v refdir="$reference_library_directory" \
    -v c_fwdseq="$DCOL_FWD_SEQ" -v c_revseq="$DCOL_REV_SEQ" -v c_marker="$DCOL_MARKER" \
    -v c_prim="$DCOL_OTU_PRIMARY" -v c_sec="$DCOL_OTU_SECONDARY" \
    -v c_tax="$DCOL_TAX_PROB" -v c_bin="$DCOL_BIN_THRESH" \
    -v n_prim="$DICTCOL_OTU_PRIMARY" -v n_sec="$DICTCOL_OTU_SECONDARY" \
    -v n_tax="$DICTCOL_TAX_PROB" -v n_bin="$DICTCOL_BIN_THRESH" '

function clean(seq) {
    gsub(/[[:space:]-]/, "", seq)
    seq = toupper(seq)
    gsub(/I/, "N", seq)
    gsub(/[^ACGTRYSWKMBDHVN]/, "N", seq)
    return seq
}

function check_prob(val, colname, row) {
    if (val < 0 || val > 1) {
        printf("ERROR: %s out of range [0,1] at row %d (value=%s)\n", colname, row, val) > "/dev/stderr"
        exit 1
    }
}

NR==1 {
    print
    next
}

{
    # sanitize sequences
    $c_fwdseq = clean($c_fwdseq)
    $c_revseq = clean($c_revseq)

    # validate numeric thresholds
    check_prob($c_prim, n_prim, NR)
    check_prob($c_sec,  n_sec,  NR)
    check_prob($c_tax,  n_tax,  NR)
    if ($c_marker == "COI-5P") check_prob($c_bin, n_bin, NR)

    print
}

' "$dictfile" > "$dictfile.tmp" && mv "$dictfile.tmp" "$dictfile") || { echo "Validation failed. Exiting."; read -p "Press Enter to exit..."; rm "$dictfile.tmp" ; exit 1; }

# ---------------------------
# Integrate primer-specific parameters into parameters TSV file
# ---------------------------

# Resolve the primer-name columns on both sides by name (pre-join positions)
DCOL_FWD_NAME=$(resolve_col "$DICTCOL_FWD_PRIMER_NAME" "$dictfile")       || exit 1
DCOL_REV_NAME=$(resolve_col "$DICTCOL_REV_PRIMER_NAME" "$dictfile")       || exit 1
PCOL_FWD_NAME=$(resolve_col "$COLNAME_FWD_PRIMER_NAME" parameters.tsv)    || exit 1
PCOL_REV_NAME=$(resolve_col "$COLNAME_REV_PRIMER_NAME" parameters.tsv)    || exit 1

awk -F'\t' -v OFS='\t' \
    -v dk1="$DCOL_FWD_NAME" -v dk2="$DCOL_REV_NAME" \
    -v pk1="$PCOL_FWD_NAME" -v pk2="$PCOL_REV_NAME" '

# ---------------------------
# Load dictionary into memory
# ---------------------------
FNR==NR {
    if (NR == 1) {
        dict_header = $0
        next
    }

    key = $dk1 "|" $dk2

    if (key in dict) {
        printf("ERROR: Duplicate primer pair in dictionary: %s\n", key) > "/dev/stderr"
        exit 1
    }

    dict[key] = $0
    next
}

# ---------------------------
# Process parameters file
# ---------------------------
FNR==1 {
    # print combined header, carrying over every dictionary column except the two primer-name keys
    split(dict_header, dh, "\t")

    dict_extra = ""
    for (i=1; i<=length(dh); i++) {
        if (i == dk1 || i == dk2) continue
        dict_extra = dict_extra OFS dh[i]
    }

    print $0 dict_extra
    next
}

{
    key = $pk1 "|" $pk2

    if (!(key in dict)) {
        printf("ERROR: Primer pair not found in dictionary: %s (row %d)\n", key, FNR) > "/dev/stderr"
        exit 1
    }

    split(dict[key], vals, "\t")

    extra = ""
    for (i=1; i<=length(vals); i++) {
        if (i == dk1 || i == dk2) continue
        extra = extra OFS vals[i]
    }

    print $0 extra
}

' $dictfile parameters.tsv > parameters_joined.tsv

mv parameters_joined.tsv parameters.tsv

# ---------------------------
# Resolve every column of the JOINED parameters.tsv by name. From here on the
# script uses these variables
# ---------------------------

COL_PLATE=$(resolve_col      "$COLNAME_PLATE"            parameters.tsv) || exit 1
COL_WELL=$(resolve_col       "$COLNAME_WELL"             parameters.tsv) || exit 1
COL_SAMPLE=$(resolve_col     "$COLNAME_SAMPLE"           parameters.tsv) || exit 1
COL_FWD_UMI=$(resolve_col    "$COLNAME_FWD_UMI"          parameters.tsv) || exit 1
COL_REV_UMI=$(resolve_col    "$COLNAME_REV_UMI"          parameters.tsv) || exit 1
COL_FWD_PRIMER=$(resolve_col "$COLNAME_FWD_PRIMER_NAME"  parameters.tsv) || exit 1
COL_REV_PRIMER=$(resolve_col "$COLNAME_REV_PRIMER_NAME"  parameters.tsv) || exit 1
COL_NEG_CTRL=$(resolve_col   "$COLNAME_NEG_CTRL"         parameters.tsv) || exit 1

COL_FWD_PRIMER_SEQ=$(resolve_col "$DICTCOL_FWD_PRIMER_SEQ" parameters.tsv) || exit 1
COL_REV_PRIMER_SEQ=$(resolve_col "$DICTCOL_REV_PRIMER_SEQ" parameters.tsv) || exit 1
COL_MARKER=$(resolve_col         "$DICTCOL_MARKER"         parameters.tsv) || exit 1
COL_REFLIB=$(resolve_col         "$DICTCOL_REFLIB"         parameters.tsv) || exit 1
COL_MINLEN=$(resolve_col         "$DICTCOL_MINLEN"         parameters.tsv) || exit 1
COL_MAXLEN=$(resolve_col         "$DICTCOL_MAXLEN"         parameters.tsv) || exit 1
COL_AMPLEN=$(resolve_col         "$DICTCOL_AMPLEN"         parameters.tsv) || exit 1
COL_OTU_PRIMARY=$(resolve_col    "$DICTCOL_OTU_PRIMARY"    parameters.tsv) || exit 1
COL_OTU_SECONDARY=$(resolve_col  "$DICTCOL_OTU_SECONDARY"  parameters.tsv) || exit 1
COL_TAX_PROB=$(resolve_col       "$DICTCOL_TAX_PROB"       parameters.tsv) || exit 1
COL_BIN_THRESH=$(resolve_col     "$DICTCOL_BIN_THRESH"     parameters.tsv) || exit 1

export COL_PLATE COL_WELL COL_SAMPLE COL_FWD_UMI COL_REV_UMI COL_FWD_PRIMER COL_REV_PRIMER \
       COL_NEG_CTRL COL_FWD_PRIMER_SEQ COL_REV_PRIMER_SEQ COL_MARKER COL_REFLIB \
       COL_MINLEN COL_MAXLEN COL_AMPLEN COL_OTU_PRIMARY COL_OTU_SECONDARY \
       COL_TAX_PROB COL_BIN_THRESH

# ---------------------------
# Extract UMIs from parameters.tsv
# ---------------------------

awk -F'\t' \
-v ov="$umi_overlap_min" \
-v err1="$error_umi1" \
-v err2="$error_umi2" \
-v c_fu="$COL_FWD_UMI" \
-v c_ru="$COL_REV_UMI" '
function revcomp(seq,    i, c, out) {

    seq = toupper(seq)
    out = ""

    for (i = length(seq); i > 0; i--) {

        c = substr(seq, i, 1)

        if      (c=="A") out = out "T"
        else if (c=="T") out = out "A"
        else if (c=="C") out = out "G"
        else if (c=="G") out = out "C"
        else out = out "N"
    }

    return out
}

NR>1 && $c_fu!="" && $c_fu!="NA" && $c_ru!="" && $c_ru!="NA" {

    fwd = toupper($c_fu)
    rev = toupper($c_ru)
    rev_rc = revcomp(rev)

    fov = int(length(fwd) * ov)
    rov = int(length(rev) * ov)

    print ">" fwd "_" rev

    print fwd ";min_overlap=" fov ";max_error_rate=" err1 \
          "..." \
          rev_rc ";min_overlap=" rov ";max_error_rate=" err2
}
' parameters.tsv > linked_umis.fasta

awk -F'\t' -v ov="$umi_overlap_min" -v err1="$error_umi1" -v c_fu="$COL_FWD_UMI" '
NR>1 && $c_fu!="" && $c_fu!="NA" {
    umi = toupper($c_fu)

    if (!seen_f[umi]++) {
        f[++n] = umi
        fov[n] = int(length(umi) * ov)
    }
}
END {
    # cutadapt reads per-adapter parameters from the SEQUENCE line, not the name.
    # Putting them in the name silently falls back to the defaults
    # (max_error_rate 0.1, min_overlap 3) and leaks them into the read name.
    for (i=1;i<=n;i++)
        print ">" f[i] "\n" f[i] ";min_overlap=" fov[i] ";max_error_rate=" err1
}
' parameters.tsv > fwd_umis.fasta

# The reverse complement is done here in awk rather than by piping through seqtk.

awk -F'\t' -v ov="$umi_overlap_min" -v err2="$error_umi2" -v c_ru="$COL_REV_UMI" '
function revcomp(seq,    i, c, out) {

    seq = toupper(seq)
    out = ""

    for (i = length(seq); i > 0; i--) {

        c = substr(seq, i, 1)

        if      (c=="A") out = out "T"
        else if (c=="T") out = out "A"
        else if (c=="C") out = out "G"
        else if (c=="G") out = out "C"
        else out = out "N"
    }

    return out
}

NR>1 && $c_ru!="" && $c_ru!="NA" {
    umi = toupper($c_ru)

    if (!seen_r[umi]++) {
        r[++n] = umi
        rov[n] = int(length(umi) * ov)
    }
}
END {
    for (i=1;i<=n;i++)
        print ">" r[i] "\n" revcomp(r[i]) ";min_overlap=" rov[i] ";max_error_rate=" err2
}
' parameters.tsv > rev_umis_rc.fasta

T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
log_step "Curating input data" "$TL" "$T2L" $(( T2 - T ))

# ---------------------------
# Merge PE reads (or multiple single-end .fastq files) into a single file and archive original .fastq.gz files
# ---------------------------

T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

# Input reads may be plain .fastq or .fastq.gz. Plain .fastq inputs sitting in the working
# directory are gzipped in place, as BIP has always done, because the primer split below
# turns every *.fastq here into a primer-pair directory. Reads in a separate --fastq
# directory are read where they are and never modified; working copies go in the working
# directory.
if [[ "$(cd "$fastq_dir" && pwd -P)" == "$(pwd -P)" ]]; then
    for i in "${!fastq_files[@]}"; do
        f="${fastq_files[$i]}"
        if [[ "$f" == *.fastq ]]; then
            gzip "$f"
            fastq_files[$i]="$f.gz"
        fi
    done
fi

if [[ "$do_pear" -eq 1 ]]; then

    echo "Paired-end mode: merging reads"

    if [[ ${#fastq_files[@]} -ne 2 ]]; then
        echo "ERROR: Paired-end mode expects exactly 2 FASTQ files (R1 and R2) in $fastq_dir, found ${#fastq_files[@]}"
        exit 1
    fi

    read1="${fastq_files[0]}"
    read2="${fastq_files[1]}"

    # Uncompressed working copies for PEAR, in the working directory (not beside the inputs)
    read1_decomp="$(basename "${read1%.gz}").decomp.fastq"
    read2_decomp="$(basename "${read2%.gz}").decomp.fastq"
    for pair in "$read1:$read1_decomp" "$read2:$read2_decomp"; do
        src="${pair%%:*}"; dst="${pair#*:}"
        if [[ "$src" == *.gz ]]; then gunzip -c "$src" > "$dst"; else cp "$src" "$dst"; fi
    done

    pear -j "$cores" -f "$read1_decomp" -r "$read2_decomp" -o "$runid" > pear.log

    merged_fastq="${runid}.assembled.fastq"

    rm -f "$read1_decomp" "$read2_decomp"

    mv "$merged_fastq" all.fastq

    # PEAR also writes the reads it discarded or could not merge. Delete them now: the
    # primer split below turns every *.fastq here into a primer-pair directory, so left
    # in place they would be processed as primer pairs. They are uncompressed and can run
    # to gigabytes on a deep run; the counts that matter for QC (how many reads were
    # discarded or left unmerged) are kept in pear.log, which goes to the results.
    rm -f "${runid}.discarded.fastq" "${runid}.unassembled.forward.fastq" "${runid}.unassembled.reverse.fastq"

    # Count raw reads
    log_count "Raw" "ALL" "all.fastq"

    T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
    log_step "PE merge" "$TL" "$T2L" $(( T2 - T ))

else

    echo "Single-end mode: concatenating all FASTQ.GZ"

    python3 - "all.fastq" "${fastq_files[@]}" <<'PYEOF'
import gzip, shutil, sys

outfile = sys.argv[1]
infiles = sys.argv[2:]

with open(outfile, "wb") as out:
    for f in infiles:
        opener = gzip.open if f.endswith(".gz") else open
        with opener(f, "rb") as src:
            shutil.copyfileobj(src, out)
PYEOF

    # Count raw reads
    log_count "Raw" "ALL" "all.fastq"

    T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
    log_step "FASTQ file merge" "$TL" "$T2L" $(( T2 - T ))


fi

# ---------------------------
# Orientate reads and split based on primers
# ---------------------------

T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

primer_list="primer_pairs.tsv"

# Build unique primer pairs
awk -F'\t' -v cfp="$COL_FWD_PRIMER" -v crp="$COL_REV_PRIMER" 'NR>1 {print $cfp "\t" $crp}' parameters.tsv \
| sort -u > "$primer_list"

# Join with dictionary.tsv → primer_sequences.tsv
# (dictionary columns resolved by name; primer_list is generated above and is always fwd-name <tab> rev-name)
awk -F'\t' -v OFS='\t' \
    -v dk1="$DCOL_FWD_NAME" -v dk2="$DCOL_REV_NAME" \
    -v dfs="$DCOL_FWD_SEQ"  -v drs="$DCOL_REV_SEQ" '
function clean(x) {
    gsub(/\r/, "", x)
    gsub(/^[ \t]+|[ \t]+$/, "", x)
    return x
}

NR==FNR {
    if (FNR == 1) next
    f = clean($dk1)
    r = clean($dk2)
    fseq[f] = clean($dfs)
    rseq[r] = clean($drs)
    next
}

{
    f = clean($1)
    r = clean($2)

    if (!(f in fseq)) {
        print "ERROR: Missing forward primer:", f > "/dev/stderr"
        exit 1
    }

    if (!(r in rseq)) {
        print "ERROR: Missing reverse primer:", r > "/dev/stderr"
        exit 1
    }

    print f, r, fseq[f], rseq[r]
}
' "$dictfile" "$primer_list" > primer_sequences.tsv

awk -F'\t' \
-v ov="$primer_overlap_min" \
-v err1="$error_primer1" \
-v err2="$error_primer2" '
function revcomp(seq,    i, c, out) {
    gsub(/\r/, "", seq)
    seq = toupper(seq)
    out = ""

    for (i = length(seq); i > 0; i--) {
        c = substr(seq, i, 1)

        if      (c == "A") out = out "T"
        else if (c == "T") out = out "A"
        else if (c == "C") out = out "G"
        else if (c == "G") out = out "C"

        else if (c == "R") out = out "Y"   # A/G -> T/C
        else if (c == "Y") out = out "R"   # C/T -> G/A
        else if (c == "S") out = out "S"   # G/C stays S
        else if (c == "W") out = out "W"   # A/T stays W
        else if (c == "K") out = out "M"   # G/T -> C/A
        else if (c == "M") out = out "K"   # A/C -> T/G

        else if (c == "B") out = out "V"   # C/G/T -> G/C/A
        else if (c == "D") out = out "H"   # A/G/T -> T/C/A
        else if (c == "H") out = out "D"   # A/C/T -> T/G/A
        else if (c == "V") out = out "B"   # A/C/G -> T/G/C

        else if (c == "N") out = out "N"

        else out = out "N"
    }
    return out
}
{
    f_name = $1
    r_name = $2
    f_seq  = $3
    r_seq  = $4

    f_ov = int(length(f_seq) * ov)
    r_ov = int(length(r_seq) * ov)

    header = ">" f_name "_" r_name

    f_part = f_seq \
             ";min_overlap=" f_ov \
             ";max_error_rate=" err1 \
             ";required"

    r_seq_rc = revcomp(r_seq)

    r_part = r_seq_rc \
             ";min_overlap=" r_ov \
             ";max_error_rate=" err2 \
             ";required"

    print header
    print f_part "..." r_part
}
' primer_sequences.tsv > linked_primers.fasta

rm primer_pairs.tsv primer_sequences.tsv

cutadapt -j "$cores" \
    -g file:linked_primers.fasta \
    --revcomp \
    --discard-untrimmed \
    --action=lowercase \
    --rename="{header}" \
    -o "{name}.fastq" \
    all.fastq

rm all.fastq

for f in *.fastq; do
    name="${f%.fastq}"
    mkdir -p "$name"
    mv "$f" "$name/"
done

# Count reads per primer pair
for d in */; do
    fq="${d%/}/${d%/}.fastq"
    [[ -f "$fq" ]] || continue
    log_count "PrimerSplit" "${d%/}" "$fq"
done

T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
log_step "Primer recognition & read orientation" "$TL" "$T2L" $(( T2 - T ))

# ---------------------------
# Process each primer directory
# ---------------------------

mkdir -m 777 Individual_Raw_Fastq_Files
process_all_primers

# ---------------------------
# Generate report
# ---------------------------
T=$(date +%s); TL=$(date "+%Y-%m-%d %H:%M:%S")

echo -e "******** Generating report for $runid"
Rscript $scripts_directory/bip4_reporting.R "$runid" "$(pwd)"

T2=$(date +%s); T2L=$(date "+%Y-%m-%d %H:%M:%S")
log_step "Report generation" "$TL" "$T2L" $(( T2 - T ))

# ---------------------------
# Write publication metrics
# ---------------------------

python3 - "$wkdir/parameters.tsv" "$wkdir/publication_metrics.txt" \
          "$runid" "$platform" "$minreads" "$min_rank_match" "$BIP_VERSION" "$parallel_sintax" \
<<'PYEOF'
import sys, csv
from datetime import date

params_file, out_file, runid, platform, minreads, min_rank_match, bip_version, parallel_sintax = sys.argv[1:9]

dict_cols = [
    'Forward Primer Name', 'Reverse Primer Name',
    'Forward Primer Sequence', 'Reverse Primer Sequence',
    'Marker', 'Reference Library',
    'Min amplicon length', 'Max amplicon length', 'Target amplicon length',
    'Primary OTU clustering threshold', 'Secondary OTU clustering threshold',
    'Tax assign probability threshold', 'BIN assign threshold'
]

seen = set()
rows = []
with open(params_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        key = (row.get('Forward Primer Name', ''), row.get('Reverse Primer Name', ''))
        if key not in seen:
            seen.add(key)
            rows.append({c: row.get(c, '') for c in dict_cols})

with open(out_file, 'w') as f:
    f.write('BIP Run Parameters\n')
    f.write('==================\n')
    f.write(f'BIP version:        {bip_version}\n')
    f.write(f'Run ID:             {runid}\n')
    f.write(f'Date:               {date.today()}\n\n')
    f.write('Command line arguments\n')
    f.write('----------------------\n')
    f.write(f'Platform:           {platform}\n')
    f.write(f'Min reads:          {minreads}\n')
    f.write(f'Min rank match:     {min_rank_match}\n')
    f.write(f'Parallel SINTAX:    {parallel_sintax}\n\n')
    f.write('Primer pair parameters\n')
    f.write('----------------------\n')
    f.write('\t'.join(dict_cols) + '\n')
    for row in rows:
        f.write('\t'.join(row[c] for c in dict_cols) + '\n')

print(f'Publication metrics written to: {out_file}')
PYEOF

# Append software versions and system info
{
    echo ""
    echo "Software versions"
    echo "-----------------"
    printf "%-20s%s\n" "vsearch:"    "$(vsearch --version 2>&1 | head -1)"
    printf "%-20s%s\n" "cutadapt:"   "$(cutadapt --version 2>/dev/null)"
    printf "%-20s%s\n" "seqtk:"      "$(seqtk 2>&1 | grep -i 'version' | head -1 | sed 's/.*[Vv]ersion:*\s*//')"
    printf "%-20s%s\n" "python3:"    "$(python3 --version 2>&1)"
    printf "%-20s%s\n" "R:"          "$(R --version 2>/dev/null | head -1)"
    echo ""
    echo "System"
    echo "------"
    if [[ "$(uname -s)" == "Darwin" ]]; then
        cpu_model=$(sysctl -n machdep.cpu.brand_string 2>/dev/null)
        ram_gb="$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 )) GB"
    else
        cpu_model=$(grep 'model name' /proc/cpuinfo | head -1 | sed 's/.*:\s*//')
        ram_gb="$(awk '/MemTotal/{printf "%.0f GB", $2/1024/1024}' /proc/meminfo)"
    fi
    printf "%-20s%s\n" "OS:"         "$(uname -sr)"
    printf "%-20s%s\n" "CPU:"        "$cpu_model"
    printf "%-20s%s\n" "RAM:"        "$ram_gb"
    printf "%-20s%s\n" "CPU threads:" "$cores (of $(getconf _NPROCESSORS_ONLN) available)"
} >> "$wkdir/publication_metrics.txt"

# ---------------------------
# Clean up
# ---------------------------
tar -czf Individual_Raw_Fastq_Files.tar.gz Individual_Raw_Fastq_Files && rm -rf Individual_Raw_Fastq_Files

mkdir -m 777 Miscelleneous_Files
mv fwd_umis.fasta rev_umis_rc.fasta linked_umis.fasta linked_primers.fasta *readcounts.tsv parameters.tsv ./Miscelleneous_Files
[[ -f "$rescue_log" ]] && mv "$rescue_log" ./Miscelleneous_Files

END_TIME=$(date +%s)
END_LABEL=$(date "+%Y-%m-%d %H:%M:%S")
log_step "Total" "$START_LABEL" "$END_LABEL" $(( END_TIME - START_TIME ))

mkdir -m 777 "$runid"_results

sweep_excludes=(! -name . ! -name "${runid}_results" ! -name "$(basename "$param_file")" ! -name "$(basename "$dictfile")" ! -name "compose.yaml")
for f in "${fastq_files[@]}"; do
    sweep_excludes+=(! -name "$(basename "$f")")
done
# --fastq / --params may point into a subfolder of the working directory (e.g. /data/reads
# when Docker mounts the user's folder at /data). Keep that top-level folder out of the
# sweep, or the user's input files would be moved into the results folder.
for p in "$fastq_dir" "$(dirname "$param_file")"; do
    if [[ "$p" == "$PWD"/* ]]; then
        rel="${p#"$PWD"/}"
        sweep_excludes+=(! -name "${rel%%/*}")
    fi
done

find . -maxdepth 1 "${sweep_excludes[@]}" -exec mv {} ./"${runid}_results"/ \;

if [[ ${#umi_warnings[@]} -gt 0 ]]; then
    echo ""
    echo "WARNING: fewer than 1% of reads matched a UMI pair for:"
    for f in "${umi_warnings[@]}"; do echo "  $f"; done
    echo "  Check the UMIs in the parameters file, especially the orientation of the Reverse UMIs."
fi

if [[ ${#bin_failures[@]} -gt 0 ]]; then
    echo ""
    echo "ERROR: the run finished, but BIN matching failed, so the BIN column is NA for:"
    for f in "${bin_failures[@]}"; do echo "  $f"; done
    echo "  Check the reference library: its headers must be 'ID|BIN;tax=...', and if a .udb"
    echo "  file sits next to it, it must be a valid vsearch UDB (delete it to have BIP rebuild it)."
    exit 1
fi
