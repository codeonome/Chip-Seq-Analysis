#!/bin/bash
#SBATCH --job-name=chip-seq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Apoorva.Sharma@nyulangone.org
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH -p cpu_medium
#SBATCH -o chip-seq_%j.log
#SBATCH -e chip-seq_%j.err

module purge
module load bowtie2/2.5.3 samtools/1.20 deeptools/3.5.1 bedtools/2.30.0

# Use pre-built index
BOWTIE_INDEX="/gpfs/share/data/bowtie2/GRCh38_noalt_as"
WORK_DIR=$(pwd)

# Sample mapping based on your data
declare -A samples
samples["SRR8511937"]="CTRL_rep1"
samples["SRR8511938"]="CTRL_rep2"
samples["SRR8511943"]="DDX6_rep1"
samples["SRR8511944"]="DDX6_rep2"

echo "=== H3K27ac ChIP-seq Processing Pipeline ==="
echo "Using pre-built index: $BOWTIE_INDEX"

# Process each sample
for srr in "${!samples[@]}"; do
    sample=${samples[$srr]}
    echo "Processing $sample ($srr)..."
    
    # QC: Count input reads
    input_reads=$(zcat ${srr}.fastq.gz | wc -l | awk '{print $1/4}')
    echo "Input reads: $input_reads"
    
    # Align with bowtie2
    echo "Aligning..."
    bowtie2 -x $BOWTIE_INDEX -U ${srr}.fastq.gz -p 8 --very-sensitive \
        2> ${sample}_align.log | samtools sort -o ${sample}_raw.bam
    
    # QC: Alignment stats
    samtools flagstat ${sample}_raw.bam > ${sample}_flagstat.txt
    mapped_reads=$(samtools view -c -F 4 ${sample}_raw.bam)
    echo "Mapped reads: $mapped_reads"
    
    # Filter: Remove duplicates, low quality, and unmapped reads
    echo "Filtering..."
    samtools view -b -q 20 -F 1028 ${sample}_raw.bam > ${sample}.bam
    samtools index ${sample}.bam
    
    # QC: Final read count
    final_reads=$(samtools view -c ${sample}.bam)
    echo "Final reads: $final_reads"
    
    # Generate normalized bigWig
    echo "Creating bigWig..."
    bamCoverage -b ${sample}.bam -o ${sample}.bw \
        --normalizeUsing RPKM -p 8 --binSize 10
    
    echo "✓ Completed $sample"
    echo "  Input: $input_reads | Mapped: $mapped_reads | Final: $final_reads"
    echo "---"
done

echo "=== All samples processed! ==="
EOF
