## Do DNMTs show sequence preferences?

### Next generation sequencing data processing

For each sample, sequencing adaptors and low quality bases were removed using cutadapt v1.16 [Martin2011](https://doi.org/10.14806/ej.17.1.200) and sickle v1.33 [Joshi and Fass 2011](https://github.com/najoshi/sickle). Read quality was checked before and after processing using FastQC v0.11.5 [Andrews et al. 2010](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 

Trimmed reads were mapped against the mm9 mouse assembly using bwa-meth v0.2.0 [Pedersen et al. 2014](https://arxiv.org/abs/1401.1129) running bwa v0.7.12-r1039 [Li 2013](https://arxiv.org/abs/1303.3997). Mapping quality was evaluated with qualimap v2.2.1 [Garcia-Alcalde et al 2012](https://academic.oup.com/bioinformatics/article/28/20/2678/206551). Duplicates were removed with Picard MarkDuplicates (picard-tools v1.96) [Broad Institute 2020](http://broadinstitute.github.io/picard/) with `REMOVE_DUPLICATES=TRUE and REMOVE_SEQUENCING_DUPLICATES=TRUE` flags.

Alignments with MAPQ > 40 were used as input to methylation calling with methyldackel v0.3.0-3-g084d926 [Ryan 2020](https://github.com/dpryan79/MethylDackel) (using HTSlib v1.2.1) in `cytosine_report` mode for all cytosines in the genome both CpG and CpH contexts.

### Data aggregation by motif

For each sample, we retrieved genome-wide reports per-cytosine. These detailed, for each cytosine, its coordinate, mapping strand, number of methylated reads and number of unmethylated reads. Next, we added the strand-aware 8-mer sequence context using the bedtools v2.27.1 [Quinlan et al 2010](https://academic.oup.com/bioinformatics/article/26/6/841/244688) suite including `slop` and `getfasta` commands. Thus, we stored the methylation status (at least a methylated read) and unmethylated (only unmethylated reads) of each cytosine and its sequence context. The procedure is available as script `extract_motifs_frequency_from_bam_binary.sh` (see below).

To aggregate the data by sequence, we generated count tables with the number of methylated and unmethylated instances. The CpG count table contained 4^6 = 4096 possible sequences (since two positions, CG, are fixed), whereas the CpH count table reported 4^6*3 = 12288 possible sequences (since the `C` position is fixed and the H can only be [A,C,T]).

### Methylation enrichment and scoring

To visualize the DNA methylation per sequence 8-mer, we developed a methylation score that represents the proportion of loci with detected methylation as compared to the total, e.g. score = $\frac{M}{M+U}$ bound to 0 (no methylated locus detected) to 1 (all 8-mer instances with at least one methylated read per locus).

For each sample, we next combined all 8-mer information to generate position weight matrices (PWMs) depicting DNA methylation preferences. To account for representation biases, we integrated both count tables of methylated and unmethylated 8-mers. First, we calculated the nucleotide frequency per position for the methylated and unmethylated 8-mer count tables separately. Second, we divided the proportion from the methylated frequencies by the unmethylated frequencies and log2-transformed the result. Hence, the score sign depicts the enrichment sign: positive values indicate methylation preference, whereas negative values suggest a trend towards unmethylation.

## Flow chart

 ![flow](./media/stranded_dnameth.png "Data flow").

## Repository structure

### Main

* `extract_motifs_frequency_from_bam_binary.sh`, bash script to generate the methylation count tables.
* `extract_motifs_frequency_from_bam.sh`, bash script to generate discretized count tables

### Mapping

*  `mapping`, bash scripts to bismark/bwa-meth retrieve, QC, map, and methylation call the single-end and paired-end reads
 * `reports_associated_to_mapping`, mainly for discovery and QC

### Accessory

* `accessory`, accessory scripts during the discovery phase
    * `.discarded`, discarded approaches
    * `cytosine_report`, first prototype

### Media

* `media`, SVG and PNG flowchart

### Exploratory reports and data (count tables)

* `rmd_reports`, PWM and other plots
    * `01_motif_extract_run_nov_2019_postproc`, report (with coverage filtering)
    * `02_motif_extract_run_nov_2019_no_coverage_filtering_postproc`, report (no coverage filtering)
    * `03_stats_assessment`, attempt to evaluate significance I
    * `04_stats_assessment`, attempt to evaluate significance II (structFDR)
    * `05_ma_plots`, visualization
    * `data/counts_nested_list.RData`, used by several reports

# Contact

* tuncay.baubec ta uzh tod ch

* izaskun.mallona ta gmail tod com
