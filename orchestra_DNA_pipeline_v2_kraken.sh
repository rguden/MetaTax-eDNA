#!/bin/bash
set -euo pipefail  # Stop on error, treat unset vars as errors, fail on pipe errors

# To run this pipeline:
#   bash orchestra_DNA_pipeline_v2_kraken.sh
# To run only specific step(s):
#   bash orchestra_DNA_pipeline_v2_kraken.sh -s qc
#   bash orchestra_DNA_pipeline_v2_kraken.sh -s qc,trim,kraken,summarize_domains,filter_bacteria,generate_bacteria_otu_R,filter_archaea,generate_archaea_otu_R

echo "Setting up Conda environment"  # (Assumes your Conda env is already activated, e.g., conda activate metagenomics)

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------
# Input directory containing raw FASTQ.gz files (both R1 and R2)
INPUT_DIR="/home/rguden/Documents/ORCHESTRA/sequences/DNA_raw_data"
# Output base directory
OUTPUT_DIR="/home/rguden/Documents/ORCHESTRA/DNA_output"
# Kraken database
KRAKEN_DB="/home/rguden/Documents/databases/nematode_db/split_marker_dbs/hybrid_18S_bac_arc_db"
# Kraken outputs
KRAKEN_DIR="$OUTPUT_DIR/community_analyses/kraken_DNA_hybrid"

# Log directory and log file
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

DATESTAMP="$(date '+%Y%m%d_%H%M%S')"   # e.g., 20251010_0832
LOG_FILE="$LOG_DIR/ORCHESTRA_DNA_pipeline_kraken_${DATESTAMP}.log"
MEM_USAGE_LOG_FILE="$LOG_DIR/ORCHESTRA_DNA_pipeline_kraken_mem_${DATESTAMP}.log"
EMAIL="rodgeemae.guden@ugent.be"

# Symlinks to the most recent logs
ln -sfn "$LOG_FILE" "$LOG_DIR/latest.log"
ln -sfn "$MEM_USAGE_LOG_FILE" "$LOG_DIR/latest_mem.log"


# ------------------------------------------------------------
# Create necessary directories
# ------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/quality_control"
mkdir -p "$OUTPUT_DIR/trimming/fastp"
mkdir -p "$OUTPUT_DIR/trimming/fastqc/fastp"
mkdir -p "$LOG_DIR"
mkdir -p "$KRAKEN_DIR"

# Parameters
THREADS=16

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
echo_with_timestamp() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}
echo_log(){ echo_with_timestamp "$1"; }

send_email() {
  local status="$1"
  local message="$2"
  # Requires local MTA or sendmail alias setup
  echo -e "Subject: Metagenomics Pipeline - $status\n\n$message" | sendmail "$EMAIL" || true
}

# Error Handling
trap 'send_email "FAILED" "Metagenomics Pipeline failed at: $(date)"; exit 1' ERR

# ------------------------------------------------------------
# Start of Pipeline
# ------------------------------------------------------------
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo_with_timestamp "Pipeline started at $START_TIME"
send_email "Metagenomics Pipeline Started" "The Metagenomics pipeline started at $START_TIME"

# ------------------------------------------------------------
# Parse Command-Line Arguments for Specific Steps
# ------------------------------------------------------------
# Default order if not specified:
# qc, trim, postqc, kraken, summarize_domains, filter_bacteria, generate_bacteria_otu_R, filter_archaea, generate_archaea_otu_R, filter_nematoda,generate_nematoda_otu_R
declare -a STEPS=()
while getopts "s:" opt; do
  case ${opt} in
    s) IFS=',' read -r -a STEPS <<< "$OPTARG" ;;
    *) echo "Usage: $0 -s <step1,step2,...>" | tee -a "$LOG_FILE"; exit 1 ;;
  esac
done

if [ ${#STEPS[@]} -eq 0 ]; then
  STEPS=(qc trim postqc kraken summarize_domains filter_bacteria generate_bacteria_otu_R filter_archaea generate_archaea_otu_R filter_nematoda generate_nematoda_otu_R)
fi

echo "Steps to run: ${STEPS[*]}" | tee -a "$LOG_FILE"

# Optional: custom condition order for summarize step (space-separated), e.g.:
# export CONDITION_ORDER="Hediste_NoNoise Hediste_Noise"
CONDITION_ORDER="${CONDITION_ORDER:-}"

# ------------------------------------------------------------
# Loop through each specified step in the order provided
# ------------------------------------------------------------
for STEP in "${STEPS[@]}"; do
  echo_with_timestamp "Starting step: $STEP"
  send_email "STEP STARTED" "Metagenomics pipeline step '$STEP' started at $(date '+%Y-%m-%d %H:%M:%S')."

  # Split the step into main and sub-step if it contains a colon (not used now, but kept for extensibility)
  IFS=':' read -r MAIN_STEP SUB_STEP <<< "$STEP"

  case "$MAIN_STEP" in
    qc)
      echo_with_timestamp "Step: Quality Control (FastQC & MultiQC)"
      fastqc --noextract -t "$THREADS" -o "$OUTPUT_DIR/quality_control" "$INPUT_DIR"/*.fastq.gz &>> "$LOG_FILE"
      multiqc "$OUTPUT_DIR/quality_control" -o "$OUTPUT_DIR/quality_control/multiqc_report" &>> "$LOG_FILE"
      echo_with_timestamp "Quality Control completed!"
      ;;

    trim)
      echo_with_timestamp "Step: Trimming the adapters (fastp)"
      for R1 in "$INPUT_DIR"/*_R1.fastq.gz; do
        [[ -f "$R1" ]] || continue
        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        SAMPLE_NAME=$(basename "$R1" _R1.fastq.gz)
        if [[ ! -f "$R2" ]]; then
          echo_with_timestamp "Skipping $SAMPLE_NAME: Missing R2 file"
          continue
        fi
        FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastp"
        FASTP_R1="${FASTP_OUTPUT}/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
        FASTP_R2="${FASTP_OUTPUT}/${SAMPLE_NAME}_trimmed_R2.fastq.gz"

        echo_with_timestamp "Processing sample $SAMPLE_NAME..."
        fastp -i "$R1" -I "$R2" \
              -o "$FASTP_R1" -O "$FASTP_R2" \
              -q 20 -l 30 --trim_poly_g --detect_adapter_for_pe \
              --thread "$THREADS" \
              -j "${FASTP_OUTPUT}/${SAMPLE_NAME}_fastp.json" \
              -h "${FASTP_OUTPUT}/${SAMPLE_NAME}_fastp.html" &>> "$LOG_FILE"
        echo_with_timestamp "Fastp trimming completed for $SAMPLE_NAME!"
      done
      echo_with_timestamp "Trimming completed!"
      ;;

    postqc)
      echo_with_timestamp "Step: Post-trimming QC (FastQC & MultiQC)"
      FASTQC_FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastqc/fastp"
      mkdir -p "$FASTQC_FASTP_OUTPUT"

      if ls "$OUTPUT_DIR/trimming/fastp/"*.fastq.gz 1> /dev/null 2>&1; then
        echo_with_timestamp "Running FastQC on trimmed reads..."
        fastqc --noextract -t "$THREADS" -o "$FASTQC_FASTP_OUTPUT" "$OUTPUT_DIR/trimming/fastp/"*.fastq.gz &>> "$LOG_FILE"

        echo_with_timestamp "Running MultiQC on trimmed reads..."
        multiqc "$FASTQC_FASTP_OUTPUT" -o "$FASTQC_FASTP_OUTPUT/multiqc_report" &>> "$LOG_FILE"
      else
        echo_with_timestamp "No trimmed files found. Skipping post-trim QC."
      fi
      echo_with_timestamp "Post-trimming QC completed!"
      ;;

    kraken)
      echo_with_timestamp "Step: Taxonomic Classification with Kraken2 (paired-end)"
      FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastp"
      mkdir -p "$KRAKEN_DIR"
      
      # How to build a kraken database
            # create db folder (e.g., mkdir -p /home/rguden/Documents/databases/silva) # also create for kraken standard
            # kraken2-build --download-taxonomy --db /home/rguden/Documents/databases/silva
            # kraken2-build --special silva --db /home/rguden/Documents/databases/silva --threads 8
            # kraken2-build --db /home/rguden/Documents/databases/silva --threads 8
            
            # kraken2-build --download-taxonomy --db /home/rguden/Documents/databases/kraken2_standard
            # kraken2-build --standard --db /home/rguden/Documents/databases/kraken2_standard --threads 8
            # kraken2-build --build --db /home/rguden/Documents/databases/kraken2_standard --threads 8
            
            # mkdir -p /home/rguden/Documents/databases/Refseq_bac # create a directory
            # kraken2-build --download-taxonomy --db /home/rguden/Documents/databases/Refseq_bac --threads 6 # Download NCBI taxonomy
            # kraken2-build --download-library bacteria --db /home/rguden/Documents/databases/Refseq_bac --threads 6 # Download RefSeq bacteria library
            # kraken2-build --build --db  /home/rguden/Documents/databases/Refseq_bac # Build the Kraken2 database


            # You can also download prebuilt databases here: https://lomanlab.github.io/mockcommunity/mc_databases.html
                # kraken2-build --download-taxonomy --db /home/rguden/Documents/databases/maxikraken
                # kraken2-build --build --db "$KRAKEN_DB" --threads "$THREADS"

      for FASTP_R1 in "$FASTP_OUTPUT"/*_trimmed_R1.fastq.gz; do
        [[ -f "$FASTP_R1" ]] || continue
        SAMPLE_NAME=$(basename "$FASTP_R1" _trimmed_R1.fastq.gz)
        FASTP_R2="${FASTP_OUTPUT}/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
        if [[ ! -f "$FASTP_R2" ]]; then
          echo_with_timestamp "Skipping $SAMPLE_NAME: Missing trimmed R2 file"
          continue
        fi

        echo_with_timestamp "Running Kraken2 on $SAMPLE_NAME (paired-end)..."
        kraken2 --db "$KRAKEN_DB" \
                --threads "$THREADS" \
                --paired \
                --report "$KRAKEN_DIR/${SAMPLE_NAME}_report.txt" \
                --output "$KRAKEN_DIR/${SAMPLE_NAME}_kraken.txt" \
                "$FASTP_R1" "$FASTP_R2" &>> "$LOG_FILE"
        echo_with_timestamp "Kraken2 classification completed for $SAMPLE_NAME!"
      done

      echo_with_timestamp "Kraken2 classification completed for all samples!"
      ;;

    summarize_domains)
      echo_log "Step: Summarize Bacteria, Archaea, Nematoda, and Unclassified from Kraken2 reports (with metadata, pre-sorted)"

      # --- Inputs/paths ---
      # KRAKEN_DIR already defined above
      METADATA_FILE="/home/rguden/Documents/ORCHESTRA/metadata/DNA/sample_condition.tsv"
      OUT_TSV="$KRAKEN_DIR/domain_summary_with_metadata.tsv"

      [[ -s "$METADATA_FILE" ]] || { echo_log "ERROR: Missing metadata $METADATA_FILE"; exit 1; }

      # --- Build pre-sorted iteration list from metadata ---
      # Expected columns: sample <TAB> condition <TAB> replicate
      ORDERED_LIST="$(mktemp)"
      if [[ -n "$CONDITION_ORDER" ]]; then
        awk -F'\t' -v order="$CONDITION_ORDER" '
          BEGIN{
            n=split(order, o, /[[:space:]]+/);
            for(i=1;i<=n;i++) rank[o[i]]=i;
          }
          NR==1{next}
          NF>=3{
            cond=$2; rep=$3; samp=$1;
            r=(cond in rank?rank[cond]:9999);
            repn = (rep ~ /^[0-9]+$/ ? rep+0 : 999999);
            printf("%09d\t%09d\t%s\t%s\t%s\n", r, repn, samp, cond, rep);
          }
        ' "$METADATA_FILE" | sort -t $'\t' -k1,1n -k2,2n | cut -f3-5 > "$ORDERED_LIST"
      else
        tail -n +2 "$METADATA_FILE" | awk -F'\t' '
          NF>=3{
            cond=$2; rep=$3; samp=$1;
            repn = (rep ~ /^[0-9]+$/ ? rep+0 : 999999);
            printf("%s\t%09d\t%s\t%s\t%s\n", cond, repn, samp, cond, rep);
          }
        ' | sort -t $'\t' -k1,1 -k2,2n | cut -f3-5 > "$ORDERED_LIST"
      fi

      # --- Header ---
      echo -e "sample\tcondition\treplicate\ttotal_reads\tunclassified_reads\tunclassified_percent\tbacteria_reads\tbacteria_percent\tarchaea_reads\tarchaea_percent\tnematoda_reads\tnematoda_percent" > "$OUT_TSV"

      # --- Helper: find report file for a sample (handles _trimmed and wildcards) ---
      find_report() {
        local s="$1" d="$2"
        local c
        for c in "$d/${s}_report.txt" \
                 "$d/${s}_trimmed_report.txt" \
                 "$d/${s}"*_report.txt; do
          [[ -f "$c" ]] && { echo "$c"; return 0; }
        done
        return 1
      }

      # --- Iterate in pre-sorted order and emit rows as we go ---
      while IFS=$'\t' read -r SAMPLE CONDITION REPLICATE; do
        [[ -n "$SAMPLE" ]] || continue
        REPFILE="$(find_report "$SAMPLE" "$KRAKEN_DIR" || true)"
        if [[ -z "$REPFILE" ]]; then
          echo_log "WARN: No report for sample $SAMPLE in $KRAKEN_DIR — skipping"
          continue
        fi

        # taxids: Bacteria=2, Archaea=2157, Eukaryota=2759 (for total fallback), Nematoda=6231, Unclassified=taxid 0 / rank U
        awk -F'\t' -v sample="$SAMPLE" -v cond="$CONDITION" -v rep="$REPLICATE" -v OFS="\t" '
          BEGIN{
            total=0
            u_reads=u_pct=0
            b_reads=b_pct=0
            a_reads=a_pct=0
            e_reads=e_pct=0
            n_reads=n_pct=0
          }
          NF<6 { next }
          {
            pct=$1; cl=$2; dr=$3; rk=$4; tx=$5; nm=$6
            gsub(/^[[:space:]]+|[[:space:]]+$/,"",pct)
            gsub(/^[[:space:]]+|[[:space:]]+$/,"",cl)
            gsub(/^[[:space:]]+|[[:space:]]+$/,"",dr)
            pct+=0; cl+=0; dr+=0
            if (tx=="1" || nm ~ /^ *root$/) total=cl
            if (tx=="0" || rk=="U") { u_reads=cl; u_pct=pct }
            if (tx=="2")     { b_reads=cl; b_pct=pct }
            if (tx=="2157")  { a_reads=cl; a_pct=pct }
            if (tx=="2759")  { e_reads=cl; e_pct=pct }
            if (tx=="6231")  { n_reads=cl; n_pct=pct }
          }
          END{
            if (total==0) total = u_reads + b_reads + a_reads + e_reads
            if (total>0) {
              if (u_pct==0 && u_reads>0) u_pct = 100.0*u_reads/total
              if (b_pct==0 && b_reads>0) b_pct = 100.0*b_reads/total
              if (a_pct==0 && a_reads>0) a_pct = 100.0*a_reads/total
              if (n_pct==0 && n_reads>0) n_pct = 100.0*n_reads/total
            }
            printf "%s\t%s\t%s\t%.0f\t%.0f\t%.2f\t%.0f\t%.2f\t%.0f\t%.2f\t%.0f\t%.2f\n", \
                   sample, cond, rep, total, u_reads, u_pct, b_reads, b_pct, a_reads, a_pct, n_reads, n_pct
          }
        ' "$REPFILE" >> "$OUT_TSV"
      done < "$ORDERED_LIST"

      rm -f "$ORDERED_LIST"
      echo_log "Wrote (pre-sorted): $OUT_TSV"
      ;;

    filter_bacteria)
      echo_log "Step: Filtering Kraken2 reports to retain only Bacteria (taxid 2)"
      FILTERED_DIR="$KRAKEN_DIR/filtered_bacteria"
      mkdir -p "$FILTERED_DIR"

      shopt -s nullglob
      any=0
      for report in "$KRAKEN_DIR"/*_report.txt; do
        [[ -f "$report" ]] || continue
        any=1
        SAMPLE_NAME=$(basename "$report" _report.txt)
        OUTFILE="$FILTERED_DIR/${SAMPLE_NAME}_report_bacteria_only.txt"

        BACT_LINE=$(awk -F'\t' '
          NF>=6 { name=$6; gsub(/^ +/, "", name); if ($5=="2" || name=="Bacteria") { print NR; exit } }' "$report")

        if [[ -z "$BACT_LINE" ]]; then
          echo_log "Skipping $SAMPLE_NAME: No Bacteria entry found"
          continue
        fi

        BACT_INDENT=$(awk -v line="$BACT_LINE" -F'\t' 'NR==line { s=$6; match(s, /^ */); print RLENGTH }' "$report")

        awk -v start="$BACT_LINE" -v indent="$BACT_INDENT" -F'\t' '
          NR==start { print; next }
          NR>start  { s=$6; match(s, /^ */); cur=RLENGTH; if (cur>indent) print; else exit }' "$report" > "$OUTFILE"

        echo_log "Filtered Bacteria-only report written: $OUTFILE"
      done
      shopt -u nullglob

      if [[ $any -eq 0 ]]; then
        echo_log "WARNING: No *_report.txt files found in $KRAKEN_DIR"
      else
        echo_log "Filtering complete for all Kraken2 reports."
      fi
      ;;

    generate_bacteria_otu_R)
      echo_log "Step: Generating genus/species OTU tables from Bacteria-only reports"
      FILTERED_DIR="$KRAKEN_DIR/filtered_bacteria"
      OTU_DIR="$FILTERED_DIR/otu_table"
      R_SCRIPT="$FILTERED_DIR/combine_bacteria_reports.R"
      mkdir -p "$OTU_DIR"

      if [[ ! -f "$R_SCRIPT" ]]; then
cat > "$R_SCRIPT" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript combine_bacteria_reports.R <filtered_dir> <otu_species.tsv> <otu_genus.tsv>", call.=FALSE)
filtered_dir    <- args[1]
otu_species_out <- args[2]
otu_genus_out   <- args[3]

read_rep <- function(f){
  x <- tryCatch(read.delim(f, header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE, fill=TRUE),
                error=function(e) NULL)
  if (is.null(x) || nrow(x)==0 || ncol(x) < 6) return(NULL)
  names(x) <- c("pct","clade_reads","direct_reads","rank","taxid","name")
  x$name <- sub("^ +","", x$name)
  x$clade_reads  <- suppressWarnings(as.numeric(x$clade_reads))
  x$direct_reads <- suppressWarnings(as.numeric(x$direct_reads))
  x$pct          <- suppressWarnings(as.numeric(gsub("%","",x$pct)))
  x$taxid        <- as.character(x$taxid)
  x
}

files <- list.files(filtered_dir, pattern="_report_bacteria_only\\.txt$", full.names=TRUE)
if (length(files)==0) stop("No filtered Bacteria reports found in: ", filtered_dir)

sample_of <- function(p) sub("_report_bacteria_only\\.txt$", "", basename(p))

long_genus   <- NULL
long_species <- NULL
for (f in files) {
  df <- read_rep(f); if (is.null(df)) next
  smp <- sample_of(f)
  gdf <- df[df$rank=="G", c("taxid","name","clade_reads")]
  sdf <- df[df$rank=="S", c("taxid","name","clade_reads")]
  if (nrow(gdf)>0) { gdf$sample <- smp; long_genus   <- rbind(long_genus, gdf) }
  if (nrow(sdf)>0) { sdf$sample <- smp; long_species <- rbind(long_species, sdf) }
}

wide_cast <- function(long_df, outfile){
  if (is.null(long_df) || nrow(long_df)==0) {
    write.table(data.frame(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
    return(invisible(NULL))
  }
  long_df$key <- paste(long_df$taxid, long_df$name, sep="|")
  taxa <- sort(unique(long_df$key))
  samples <- sort(unique(long_df$sample))
  mat <- matrix(0, nrow=length(taxa), ncol=length(samples),
                dimnames=list(taxa, samples))
  for (i in seq_len(nrow(long_df))) {
    k <- long_df$key[i]; s <- long_df$sample[i]; v <- long_df$clade_reads[i]
    if (!is.na(v)) mat[k, s] <- mat[k, s] + v
  }
  out <- data.frame(taxon_key=rownames(mat), mat, check.names=FALSE)
  sp  <- strsplit(out$taxon_key, "\\|")
  out$taxid <- vapply(sp, function(z) z[1], FUN.VALUE=character(1))
  out$name  <- vapply(sp, function(z) paste(z[-1], collapse="|"), FUN.VALUE=character(1))
  out$taxon_key <- NULL
  out <- out[, c("taxid","name", samples), drop=FALSE]
  write.table(out, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
  invisible(out)
}

otuS <- wide_cast(long_species, otu_species_out)
otuG <- wide_cast(long_genus,   otu_genus_out)
cat("Wrote:\n  -", otu_species_out, "\n  -", otu_genus_out, "\n")
RSCRIPT
      fi

      OTU_SPECIES="$OTU_DIR/otu_bacteria_species.tsv"
      OTU_GENUS="$OTU_DIR/otu_bacteria_genus.tsv"

      echo_log "Running R to build Bacteria genus/species OTU tables…"
      Rscript "$R_SCRIPT" "$FILTERED_DIR" "$OTU_SPECIES" "$OTU_GENUS" \
        && echo_log "? Bacteria OTUs written to $OTU_DIR" \
        || { echo_log "? R step failed (Bacteria)"; exit 1; }
      ;;

    filter_archaea)
      echo_log "Step: Filtering Kraken2 reports to retain only Archaea (taxid 2157)"
      FILTERED_DIR="$KRAKEN_DIR/filtered_archaea"
      mkdir -p "$FILTERED_DIR"

      shopt -s nullglob
      any=0
      for report in "$KRAKEN_DIR"/*_report.txt; do
        [[ -f "$report" ]] || continue
        any=1
        SAMPLE_NAME=$(basename "$report" _report.txt)
        OUTFILE="$FILTERED_DIR/${SAMPLE_NAME}_report_archaea_only.txt"

        ARCH_LINE=$(awk -F'\t' '
          NF>=6 { name=$6; gsub(/^ +/, "", name); if ($5=="2157" || name=="Archaea") { print NR; exit } }' "$report")

        if [[ -z "$ARCH_LINE" ]]; then
          echo_log "Skipping $SAMPLE_NAME: No Archaea entry found"
          continue
        fi

        ARCH_INDENT=$(awk -v line="$ARCH_LINE" -F'\t' 'NR==line { s=$6; match(s, /^ */); print RLENGTH }' "$report")

        awk -v start="$ARCH_LINE" -v indent="$ARCH_INDENT" -F'\t' '
          NR==start { print; next }
          NR>start  { s=$6; match(s, /^ */); cur=RLENGTH; if (cur>indent) print; else exit }' "$report" > "$OUTFILE"

        echo_log "Filtered Archaea-only report written: $OUTFILE"
      done
      shopt -u nullglob

      if [[ $any -eq 0 ]]; then
        echo_log "WARNING: No *_report.txt files found in $KRAKEN_DIR"
      else
        echo_log "Filtering complete for all Kraken2 reports."
      fi
      ;;

    generate_archaea_otu_R)
      echo_log "Step: Generating genus/species OTU tables from Archaea-only reports"
      FILTERED_DIR="$KRAKEN_DIR/filtered_archaea"
      OTU_DIR="$FILTERED_DIR/otu_table"
      R_SCRIPT="$FILTERED_DIR/combine_archaea_reports.R"
      mkdir -p "$OTU_DIR"

      if [[ ! -f "$R_SCRIPT" ]]; then
cat > "$R_SCRIPT" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript combine_archaea_reports.R <filtered_dir> <otu_species.tsv> <otu_genus.tsv>", call.=FALSE)
filtered_dir    <- args[1]
otu_species_out <- args[2]
otu_genus_out   <- args[3]

read_rep <- function(f){
  x <- tryCatch(read.delim(f, header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE, fill=TRUE),
                error=function(e) NULL)
  if (is.null(x) || nrow(x)==0 || ncol(x) < 6) return(NULL)
  names(x) <- c("pct","clade_reads","direct_reads","rank","taxid","name")
  x$name <- sub("^ +","", x$name)
  x$clade_reads  <- suppressWarnings(as.numeric(x$clade_reads))
  x$direct_reads <- suppressWarnings(as.numeric(x$direct_reads))
  x$pct          <- suppressWarnings(as.numeric(gsub("%","",x$pct)))
  x$taxid        <- as.character(x$taxid)
  x
}

files <- list.files(filtered_dir, pattern="_report_archaea_only\\.txt$", full.names=TRUE)
if (length(files)==0) stop("No filtered Archaea reports found in: ", filtered_dir)

sample_of <- function(p) sub("_report_archaea_only\\.txt$", "", basename(p))

long_genus   <- NULL
long_species <- NULL
for (f in files) {
  df <- read_rep(f); if (is.null(df)) next
  smp <- sample_of(f)
  gdf <- df[df$rank=="G", c("taxid","name","clade_reads")]
  sdf <- df[df$rank=="S", c("taxid","name","clade_reads")]
  if (nrow(gdf)>0) { gdf$sample <- smp; long_genus   <- rbind(long_genus, gdf) }
  if (nrow(sdf)>0) { sdf$sample <- smp; long_species <- rbind(long_species, sdf) }
}

wide_cast <- function(long_df, outfile){
  if (is.null(long_df) || nrow(long_df)==0) {
    write.table(data.frame(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
    return(invisible(NULL))
  }
  long_df$key <- paste(long_df$taxid, long_df$name, sep="|")
  taxa <- sort(unique(long_df$key))
  samples <- sort(unique(long_df$sample))
  mat <- matrix(0, nrow=length(taxa), ncol=length(samples),
                dimnames=list(taxa, samples))
  for (i in seq_len(nrow(long_df))) {
    k <- long_df$key[i]; s <- long_df$sample[i]; v <- long_df$clade_reads[i]
    if (!is.na(v)) mat[k, s] <- mat[k, s] + v
  }
  out <- data.frame(taxon_key=rownames(mat), mat, check.names=FALSE)
  sp  <- strsplit(out$taxon_key, "\\|")
  out$taxid <- vapply(sp, function(z) z[1], FUN.VALUE=character(1))
  out$name  <- vapply(sp, function(z) paste(z[-1], collapse="|"), FUN.VALUE=character(1))
  out$taxon_key <- NULL
  out <- out[, c("taxid","name", samples), drop=FALSE]
  write.table(out, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
  invisible(out)
}

otuS <- wide_cast(long_species, otu_species_out)
otuG <- wide_cast(long_genus,   otu_genus_out)
cat("Wrote:\n  -", otu_species_out, "\n  -", otu_genus_out, "\n")
RSCRIPT
      fi

      OTU_SPECIES="$OTU_DIR/otu_archaea_species.tsv"
      OTU_GENUS="$OTU_DIR/otu_archaea_genus.tsv"

      echo_log "Running R to build Archaea genus/species OTU tables…"
      Rscript "$R_SCRIPT" "$FILTERED_DIR" "$OTU_SPECIES" "$OTU_GENUS" \
        && echo_log "? Archaea OTUs written to $OTU_DIR" \
        || { echo_log "? R step failed (Archaea)"; exit 1; }
      ;;
      
    filter_nematoda)
      echo_log "Step: Filtering Kraken2 reports to retain only Nematoda (taxid 6231)"

      # Uses the unified Kraken output dir
      FILTERED_DIR="$KRAKEN_DIR/filtered_nematoda"
      mkdir -p "$FILTERED_DIR"

      shopt -s nullglob
      any=0
      for report in "$KRAKEN_DIR"/*_report.txt; do
        [[ -f "$report" ]] || continue
        any=1

        SAMPLE_NAME=$(basename "$report" _report.txt)
        OUTFILE="$FILTERED_DIR/${SAMPLE_NAME}_report_nematoda_only.txt"

        # 1) find the Nematoda row (taxid 6231 or name "Nematoda")
        NEMA_LINE=$(awk -F'\t' '
          NF>=6 {
            name=$6; gsub(/^ +/, "", name)
            if ($5=="6231" || name=="Nematoda") { print NR; exit }
          }' "$report")

        if [[ -z "$NEMA_LINE" ]]; then
          echo_log "Skipping $SAMPLE_NAME: No Nematoda entry found"
          continue
        fi

        # 2) get indentation (leading spaces in column 6)
        NEMA_INDENT=$(awk -v line="$NEMA_LINE" -F'\t' '
          NR==line { s=$6; match(s, /^ */); print RLENGTH }' "$report")

        # 3) print Nematoda row + all deeper-indented descendants
        awk -v start="$NEMA_LINE" -v indent="$NEMA_INDENT" -F'\t' '
          NR==start { print; next }
          NR>start  {
            s=$6; match(s, /^ */); cur=RLENGTH
            if (cur>indent) print; else exit
          }' "$report" > "$OUTFILE"

        echo_log "Filtered Nematoda-only report written: $OUTFILE"
      done
      shopt -u nullglob

      if [[ $any -eq 0 ]]; then
        echo_log "WARNING: No *_report.txt files found in $KRAKEN_DIR"
      else
        echo_log "Filtering complete for all Kraken2 reports."
      fi
      ;;

    generate_nematoda_otu_R)
      echo_log "Step: Generating genus/species OTU tables from Nematoda-only reports"

      FILTERED_DIR="$KRAKEN_DIR/filtered_nematoda"
      OTU_DIR="$FILTERED_DIR/otu_table"
      R_SCRIPT="$FILTERED_DIR/combine_nematoda_reports.R"
      mkdir -p "$OTU_DIR"

      # Create the R script once if it doesn't exist
      if [[ ! -f "$R_SCRIPT" ]]; then
cat > "$R_SCRIPT" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript combine_nematoda_reports.R <filtered_dir> <otu_species.tsv> <otu_genus.tsv>", call.=FALSE)
}
filtered_dir    <- args[1]
otu_species_out <- args[2]
otu_genus_out   <- args[3]

# Pure base R to avoid extra dependencies
read_rep <- function(f){
  # Kraken2 filtered block: 6 tab columns
  # 1=%total, 2=clade_reads, 3=direct_reads, 4=rank, 5=taxid, 6=indented name
  x <- tryCatch(
    read.delim(f, header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE, fill=TRUE),
    error=function(e) NULL
  )
  if (is.null(x) || nrow(x)==0 || ncol(x) < 6) return(NULL)
  names(x) <- c("pct","clade_reads","direct_reads","rank","taxid","name")
  x$name <- sub("^ +","", x$name)
  x$clade_reads  <- suppressWarnings(as.numeric(x$clade_reads))
  x$direct_reads <- suppressWarnings(as.numeric(x$direct_reads))
  x$pct          <- suppressWarnings(as.numeric(gsub("%","",x$pct)))
  x$taxid        <- as.character(x$taxid)
  x
}

files <- list.files(filtered_dir, pattern="_report_nematoda_only\\.txt$", full.names=TRUE)
if (length(files)==0) stop("No filtered Nematoda reports found in: ", filtered_dir)

sample_of <- function(p) sub("_report_nematoda_only\\.txt$", "", basename(p))

long_genus   <- NULL
long_species <- NULL

for (f in files) {
  df <- read_rep(f)
  if (is.null(df)) next
  smp <- sample_of(f)

  # Use clade reads so higher ranks include descendants
  gdf <- df[df$rank=="G", c("taxid","name","clade_reads")]
  sdf <- df[df$rank=="S", c("taxid","name","clade_reads")]

  if (nrow(gdf)>0) { gdf$sample <- smp; long_genus   <- rbind(long_genus, gdf) }
  if (nrow(sdf)>0) { sdf$sample <- smp; long_species <- rbind(long_species, sdf) }
}

# Wide-cast using base R (rows = taxid|name, cols = samples)
wide_cast <- function(long_df, outfile){
  if (is.null(long_df) || nrow(long_df)==0) {
    write.table(data.frame(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
    return(invisible(NULL))
  }
  long_df$key <- paste(long_df$taxid, long_df$name, sep="|")
  taxa    <- sort(unique(long_df$key))
  samples <- sort(unique(long_df$sample))
  mat <- matrix(0, nrow=length(taxa), ncol=length(samples),
                dimnames=list(taxa, samples))
  for (i in seq_len(nrow(long_df))) {
    k <- long_df$key[i]; s <- long_df$sample[i]; v <- long_df$clade_reads[i]
    if (!is.na(v)) mat[k, s] <- mat[k, s] + v
  }
  out <- data.frame(taxon_key=rownames(mat), mat, check.names=FALSE)
  sp  <- strsplit(out$taxon_key, "\\|")
  out$taxid <- vapply(sp, function(z) z[1], FUN.VALUE=character(1))
  out$name  <- vapply(sp, function(z) paste(z[-1], collapse="|"), FUN.VALUE=character(1))
  out$taxon_key <- NULL
  out <- out[, c("taxid","name", samples), drop=FALSE]
  write.table(out, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
  invisible(out)
}

otuS <- wide_cast(long_species, otu_species_out)
otuG <- wide_cast(long_genus,   otu_genus_out)

# Optional: taxonomy map with rank label
rank_map <- rbind(
  if (!is.null(otuG)) data.frame(taxid=otuG$taxid, rank="G", name=otuG$name, stringsAsFactors=FALSE) else NULL,
  if (!is.null(otuS)) data.frame(taxid=otuS$taxid, rank="S", name=otuS$name, stringsAsFactors=FALSE) else NULL
)
if (!is.null(rank_map) && nrow(rank_map) > 0) {
  rank_map <- unique(rank_map)
  write.table(rank_map, file=file.path(filtered_dir,"taxonomy_map.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
}

cat("Wrote:\n  -", otu_species_out, "\n  -", otu_genus_out,
    "\n  -", file.path(filtered_dir,"taxonomy_map.tsv"), "\n")
RSCRIPT
      fi

      OTU_SPECIES="$OTU_DIR/otu_table_species.tsv"
      OTU_GENUS="$OTU_DIR/otu_table_genus.tsv"

      echo_log "Running R to build genus/species OTU tables…"
      Rscript "$R_SCRIPT" "$FILTERED_DIR" "$OTU_SPECIES" "$OTU_GENUS" \
        && echo_log "? OTU tables written to $OTU_DIR" \
        || { echo_log "? R step failed"; exit 1; }

      send_email "Nematoda OTU Table Complete" "OTU tables (genus, species) written to $OTU_DIR"
      ;;


    *)
      echo_with_timestamp "Unknown step: $MAIN_STEP"
      exit 1
      ;;
  esac

  echo_with_timestamp "Completed step: $STEP"
  send_email "STEP COMPLETED" "Metagenomics pipeline step '$STEP' completed at $(date '+%Y-%m-%d %H:%M:%S')."
done

END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo_with_timestamp "Pipeline completed at $END_TIME"
send_email "COMPLETED" "Metagenomics Pipeline completed at $END_TIME"
