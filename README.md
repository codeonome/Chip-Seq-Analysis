# Chip-Seq-Analysis

**H3K27ac ChIP-seq Analysis at Super-Enhancer Regions**
Objective:

Reproduce Figure 6D from Di Stefano et al.: "The RNA Helicase DDX6 Controls Cellular Plasticity by Modulating P-Body Homeostasis" Cell Stem Cell 25(5), 622-638.e13 (2019).
Analysis Focus
Examine how H3K27ac levels (measured by ChIP-sequencing in human embryonic stem cells) change at known super-enhancer regions following DDX6 knockdown.
Data Sources
• ChIP-seq data: Available as GSE dataset from the Di Stefano et al. paper
• Super-enhancer coordinates: From Hnisz et al.: "Super-Enhancers in the Control of Cell Identity and Disease" Cell 155(4), 934-947 (2013)

**Deliverables**

1. Reproduction of Figure 6D
2. Complete code and/or pipeline used for the analysis and figure generation

**FILES**
1. Aligning the raw fast files:
The full pipeline till alignment is is process-chip-seq.sh

2. Downloading the bed files and super enhancer coordinates:
downloaded from the reference papers

3. This is just an intermediate step I had to do:
The original H1.bed file contained coordinates in hg19 format (from Hnisz et al. 2013)
I aligned our FASTQ data to hg38 reference genome using the bowtie2 index
This created a coordinate system mismatch - I was looking for hg19 coordinates in hg38-aligned data
The liftOver conversion fixed this mismatch by converting coordinates from hg19 to hg38
FASTQ files → Bowtie2 (hg38 reference) → BAM files (hg38 coordinates)
                                              ↓
H1.bed (hg19) → LiftOver → H1_hg38.bed → Signal quantification → CORRECTED results

Re-quantify signal with CORRECT hg38 coordinates
echo "=== RE-RUNNING ANALYSIS WITH CORRECT COORDINATES ==="

multiBigwigSummary BED-file --BED H1_hg38_no_chr.bed \
    -b /gpfs/data/khodadadilab/home/temp/Di-Stefano-Lab-Assignment/Task-2/align/CTRL_rep1.bw \
       /gpfs/data/khodadadilab/home/temp/Di-Stefano-Lab-Assignment/Task-2/align/CTRL_rep2.bw \
       /gpfs/data/khodadadilab/home/temp/Di-Stefano-Lab-Assignment/Task-2/align/DDX6_rep1.bw \
       /gpfs/data/khodadadilab/home/temp/Di-Stefano-Lab-Assignment/Task-2/align/DDX6_rep2.bw \
    -o h3k27ac_CORRECTED_signals.npz \
    --outRawCounts h3k27ac_CORRECTED_data.tab

echo "Corrected analysis completed!"

4. plotting the results:
Results were plotted using the file plot.r
