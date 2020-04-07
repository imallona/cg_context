# Do DNMT show sequence preferences?

## Approach


### Next generation sequencing data processing

For each sample, sequencing adaptors and low quality bases were removed using cutadapt v1.16 and sickle v1.33. Read quality was checked before and after processing using FastQC v0.11.5. 

Trimmed reads were mapped against the mm9 mouse assembly using bwa-meth v0.2.0 (bwa v0.7.12-r1039). Mapping quality was evaluated with qualimap v2.2.1. Duplicates were removed with Picard MarkDuplicates (picard-tools v1.96) with `REMOVE_DUPLICATES=TRUE and REMOVE_SEQUENCING_DUPLICATES=TRUE` flags.

Alignments with MAPQ > 40 were used as input to methylation calling with methyldackel v0.3.0-3-g084d926 (using HTSlib v1.2.1) in `cytosine_report` mode for all cytosines in the genome both CpG and CpH contexts.

### Data aggregation by motif

For each sample, genome-wide reports of the per-cytosine methylation status (e.g. with a record per cytosine including its coordinate, mapping strand, number of methylated reads, number of unmethylated reads) were obtained. Next, we added the strand-aware 8-mer sequence context using the bedtools v2.27.1 suite including `slop` and `getfasta` commands. Thus, we stored the methylation status (at least a methylated read) and unmethylated (only unmethylated reads) of each cytosine and its sequence context. Script XXX.

To aggregate the data by sequence, we generated count tables with the number of methylated and unmethylated instances. The CpG count table contained 4^6 = 4096 possible sequences (since two positions, CG, are fixed), whereas the CpH count table reported 4^6*3 = 12288 possible sequences (since the `C` position is fixed and the H can only be [A,C,T]).

### Methylation enrichment and scoring

To visualize the DNA methylation per sequence 8-mer, we developed a methylation score that represents the proportion of loci with detected methylation as compared to the total, e.g. score = $\frac{M}{M+U}$ bound to 0 (no methylated locus detected) to 1 (all 8-mer instances with at least one methylated read per locus).

We next combined all 8-mer information to generate position weight matrices (PWMs) depicting DNA methylation preferences. To account for representation biases, we integrated the count tables of methylated and unmethylated 8-mers. First, we calculated the nucleotide frequency per position for the methylated and unmethylated 8-mer count tables separately. Second, we divided the proportion from the methylated frequencies by the unmethylated frequencies and log2-transformed the result. Hence, the score sign depicts the enrichment sign: positive values indicate methylation preference, and negative values trend towards unmethylation.

## Repository structure

 ![flow](./media/stranded_dnameth.png "Data flow").


1. `extract_motifs_frequency_from_bam_binary.sh`, bash script to generate the methylation count tables.
1. `extract_motifs_frequency_from_bam.sh`, bash script to generate discretized count tables
1. `accessory`, accessory scripts during the discovery phase
   * `.discarded`, discarded approaches
   * `cytosine_report`, first prototype
1.  `mapping`, bash scripts to bismark/bwa-meth retrieve, QC, map, and methylation call the single-end and paired-end reads
      * `reports_associated_to_mapping`, mainly for discovery and QC
1. `media`, SVG and PNG flowchart
1. `rmd_reports`, PWM and other plots
   * `01_motif_extract_run_nov_2019_postproc`, report (with coverage filtering)
   * `02_motif_extract_run_nov_2019_no_coverage_filtering_postproc`, report (no coverage filtering)
   * `03_stats_assessment`, attempt to evaluate significance I
   * `04_stats_assessment`, attempt to evaluate significance II (structFDR)
   * `05_ma_plots`, visualization
   * `data/counts_nested_list.RData`, used by several reports

# Contact

tuncay.baubec ta uzh tod ch
izaskun.mallona ta gmail tod com
