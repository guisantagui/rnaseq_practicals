#!/bin/bash

# Variables
####################################################################################################

min_read_size=50
min_phred=30
strand="reverse"

# Directory stuff
####################################################################################################
cd ..
#source rnaseq_preproc/bin/activate

accessions_file="accessions_all.txt"
data_dir="data/"
ref_dir="${data_dir}reference/"
idx_dir="${ref_dir}PA14_idx/"
results_dir="results/preprocessing/"

qc_pretrim="${results_dir}qc_pretrim/"
qc_postrim="${results_dir}qc_postrim/"
trim_dir="${results_dir}trimmed/"
align_dir="${results_dir}aligned/"
counts_dir="${results_dir}counts/"

mkdir -p $ref_dir
mkdir -p $results
mkdir -p $qc_pretrim
mkdir -p $qc_postrim
mkdir -p $trim_dir
mkdir -p $idx_dir
mkdir -p $align_dir
mkdir -p $counts_dir

genome_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/625/GCA_000014625.1_ASM1462v1/GCA_000014625.1_ASM1462v1_genomic.fna.gz"
annote_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/625/GCA_000014625.1_ASM1462v1/GCA_000014625.1_ASM1462v1_genomic.gtf.gz"

genome_file=$(basename $genome_url)
annote_file=$(basename $annote_url)

genome_file="${genome_file%.gz}"
annote_file="${annote_file%.gz}"

####################################################################################################
# Pipeline                                                                                         #
####################################################################################################

# Download references
####################################################################################################

cd $ref_dir
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/014/625/GCF_000014625.1_ASM1462v1/GCF_000014625.1_ASM1462v1_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/014/625/GCF_000014625.1_ASM1462v1/GCF_000014625.1_ASM1462v1_genomic.gtf.gz

wget $genome_url
wget $annote_url

gunzip *.fna.gz
gunzip *.gtf.gz
cd ../..

# Download FASTQs
####################################################################################################
cd $data_dir

while read -r accession; do
  fastq_1="${accession}_1.fastq"
  fastq_2="${accession}_2.fastq"
  if [[ ! -f "$fastq_1" || ! -f "$fastq_2" ]]; then
    echo "Downloading accession: $accession"
    prefetch $accession
    fasterq-dump $accession
  else
    echo "$accession has already been dumped. Skipping."
  fi
done < $accessions_file

cd ..

# Preprocessing
####################################################################################################

# Create an index with the genome file
if [ ! "$(ls -A $idx_dir)" ]; then
  STAR --runThreadN 8 \
       --runMode genomeGenerate \
       --genomeDir $idx_dir \
       --genomeFastaFiles "${ref_dir}${genome_file}" \
       --sjdbGTFfile "${ref_dir}${annote_file}" \
       --sjdbOverhang 72
fi

# Read accessions.txt and for each accession run fastqc, trim reads and align with STAR
while read -r acc; do
  fastqc_out_1="${qc_pretrim}${acc}_1_fastqc.html"
  fastqc_out_2="${qc_pretrim}${acc}_2_fastqc.html"
  if [[ ! -f "$fastqc_out_1" || ! -f "$fastqc_out_2" ]]; then
    # Run FASTQC
    echo "Running FASTQC for $acc"
    fastqc "${data_dir}${acc}_1.fastq" "${data_dir}${acc}_2.fastq"  -o $qc_pretrim &
  else
    echo "FASTQC already done for: $acc. Skipping."
  fi
done < "${data_dir}${accessions_file}"
wait

while read -r acc; do
  trimmed_1="${trim_dir}/${acc}_1_val_1.fq"
  trimmed_2="${trim_dir}/${acc}_2_val_2.fq"
  if [[ ! -f "$trimmed_1" || ! -f "$trimmed_2" ]]; then
    # Trim
    echo "Trimming reads for: $acc"
    trim_galore --paired --quality $min_phred --length $min_read_size --phred33 --fastqc \
                --fastqc_args "-o ${qc_postrim}" -o $trim_dir \
                "${data_dir}${acc}_1.fastq" "${data_dir}${acc}_2.fastq" &
  else
    echo "Trimmed files already exist for: $acc. Skipping trimming."
  fi
done < "${data_dir}${accessions_file}"
wait

# Align trimmed reads to the index
while read -r acc; do
  echo "Aligning ${acc} to reference genome..."
  align_samp_dir="${align_dir}${acc}/"
  mkdir -p $align_samp_dir
  aligned_bam="${align_samp_dir}/${acc}_Aligned.sortedByCoord.out.bam"
  if [[ -f "$aligned_bam" ]]; then
    echo "Alignment already done for ${acc}. Skipping alignment."
  else
    echo "Aligning ${acc} to reference genome..."
    STAR --genomeDir $idx_dir \
         --readFilesIn "${trim_dir}${acc}_1_val_1.fq" "${trim_dir}${acc}_2_val_2.fq"\
         --runThreadN 10 \
         --outFileNamePrefix "${align_samp_dir}${acc}_" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within KeepPairs &
  fi
done < "${data_dir}${accessions_file}"
wait

# Obtain alignment stats with picard
while read -r acc; do
  echo "Obtaining alignment stats of ${acc} with Picard..."
  align_samp_dir="${align_dir}${acc}/"
  # Perform some QCs
  picard CollectAlignmentSummaryMetrics \
         R="${ref_dir}${genome_file}" \
         I="${align_samp_dir}${acc}_Aligned.sortedByCoord.out.bam" \
         O="${align_samp_dir}${acc}_picard_alignment_summary_metrics.txt"

  picard CollectInsertSizeMetrics \
         I="${align_samp_dir}${acc}_Aligned.sortedByCoord.out.bam" \
         O="${align_samp_dir}${acc}_picard_insert_size_metrics.txt" \
         H="${align_samp_dir}${acc}_insert_size_histogram.pdf" \
         M=0.5
done < "${data_dir}${accessions_file}"

# Obtain counts from the reference alignments and parse them in a CSV matrix
Rscript scripts/get_counts.R $align_dir \
        --annotFile "${ref_dir}${annote_file}" \
        --att "locus_tag" \
        --strand $strand \
        --outDir $counts_dir 

# Obtain MultiQC reports
multiqc $(find $qc_pretrim -type f -name "*.zip" -not -name "*_val_*") \
        --outdir "$qc_pretrim" -f

multiqc $(find $qc_postrim -type f -name "*.zip" -and -name "*_val_*") \
        $(find $align_ref_dir -type f -name "*_Log.final.out") \
        $(find $align_ref_dir -type f -name "*picard*") \
        --outdir "$qc_postrim" -f

# Differential expression
####################################################################################################
DE_dir="results/DE/"
mkdir -p $DE_dir

Rscript scripts/DE.R --gene_counts "${counts_dir}counts.csv" \
        --sample_info "${data_dir}sample_info.csv" \
        --outDir "$DE_dir"

# Functional enrichment
####################################################################################################
enrich_dir="results/enrichment/"
mkdir -p $enrich_dir

Rscript scripts/enrichment.R "${DE_dir}w_ercc/DESeqRes_LFCShrunk.csv" \
        --alpha 0.05 \
        --outDir "${enrich_dir}w_ercc/"

Rscript scripts/enrichment.R "${DE_dir}wo_ercc/DESeqRes_LFCShrunk.csv" \
        --alpha 0.05 \
        --outDir "${enrich_dir}wo_ercc/"