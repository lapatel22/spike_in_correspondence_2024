spike_correspondence_figure1b
================

- <a
  href="#estimatinginterphasemitotic-h3k9ac-with-mass-spectrometry-data"
  id="toc-estimatinginterphasemitotic-h3k9ac-with-mass-spectrometry-data">EstimatingInterphase/Mitotic
  H3K9ac with Mass Spectrometry Data</a>
- <a href="#chip-seq-titration-of-h3k9ac"
  id="toc-chip-seq-titration-of-h3k9ac">ChIP-seq Titration of H3K9ac</a>
  - <a href="#trimming-.fastq-files" id="toc-trimming-.fastq-files">Trimming
    .FASTQ Files</a>
  - <a href="#alignment-.fastq---.sam"
    id="toc-alignment-.fastq---.sam">Alignment: .FASTQ -&gt; .SAM</a>
    - <a href="#genome-preparation" id="toc-genome-preparation">Genome
      Preparation:</a>
  - <a href="#remove-pcr-duplicates" id="toc-remove-pcr-duplicates">Remove
    PCR Duplicates</a>
  - <a href="#separate-alignment-file-into-species-specific-alignments"
    id="toc-separate-alignment-file-into-species-specific-alignments">Separate
    Alignment File into Species-specific Alignments</a>
  - <a href="#make-homer-tag-directories"
    id="toc-make-homer-tag-directories">Make HOMER Tag Directories</a>
  - <a href="#visualize-with-bigwigs"
    id="toc-visualize-with-bigwigs">Visualize with BigWigs</a>
  - <a href="#determine-normalization-factor"
    id="toc-determine-normalization-factor">Determine Normalization
    Factor</a>
  - <a href="#analysis-with-homer" id="toc-analysis-with-homer">Analysis
    with HOMER</a>
    - <a href="#histograms-at-tss" id="toc-histograms-at-tss">Histograms at
      TSS</a>
    - <a href="#peak-finding" id="toc-peak-finding">Peak Finding</a>
    - <a href="#quantification-at-peaks"
      id="toc-quantification-at-peaks">Quantification at Peaks</a>
  - <a href="#megapeak-get-h3k9ac-regions-to-quantify-signal"
    id="toc-megapeak-get-h3k9ac-regions-to-quantify-signal">Megapeak: Get
    H3K9ac regions to quantify signal</a>
- <a href="#import-histogram-data" id="toc-import-histogram-data">Import
  histogram data</a>
  - <a href="#calculate-signal-at-peaks-area-under-histogram-curve"
    id="toc-calculate-signal-at-peaks-area-under-histogram-curve">Calculate
    Signal at Peaks: Area under Histogram Curve</a>
  - <a href="#make-line-of-expected-signal"
    id="toc-make-line-of-expected-signal">Make line of expected signal:</a>
  - <a href="#scale-signal-set-minimum-signal-to-0-maximum-signal-to-1"
    id="toc-scale-signal-set-minimum-signal-to-0-maximum-signal-to-1">scale
    signal, set minimum signal to 0, maximum signal to 1</a>
  - <a href="#make-line-of-expected-signal-1"
    id="toc-make-line-of-expected-signal-1">Make line of expected
    signal:</a>
  - <a href="#plot-read-normalized-signal-vs-expected-signal-fig-1b"
    id="toc-plot-read-normalized-signal-vs-expected-signal-fig-1b">Plot Read
    normalized signal vs expected signal: Fig 1b</a>
- <a href="#normalize-peak-signal-to-dual-spike-ins"
  id="toc-normalize-peak-signal-to-dual-spike-ins">Normalize peak signal
  to Dual Spike-ins</a>
  - <a href="#determine-normalization-factor-from-aligned-reads"
    id="toc-determine-normalization-factor-from-aligned-reads">Determine
    normalization factor from aligned reads:</a>
  - <a href="#normalize-to-fly-ipinput"
    id="toc-normalize-to-fly-ipinput">Normalize to fly IP/input</a>
    - <a href="#calculate-signal-at-peaks-area-under-histogram-curve-1"
      id="toc-calculate-signal-at-peaks-area-under-histogram-curve-1">Calculate
      Signal at Peaks: Area under Histogram Curve</a>
    - <a
      href="#scale-signal-at-peaks-to-be-from-0-1-plot-this-normalized-signal-against-expected-signal"
      id="toc-scale-signal-at-peaks-to-be-from-0-1-plot-this-normalized-signal-against-expected-signal">Scale
      signal at peaks to be from 0-1, plot this normalized signal against
      expected signal</a>
    - <a href="#plot-read-normalized-signal-vs-expected-signal-fig-1b-1"
      id="toc-plot-read-normalized-signal-vs-expected-signal-fig-1b-1">Plot
      Read normalized signal vs expected signal: Fig 1b</a>
  - <a href="#normalize-to-yeast-spike-in"
    id="toc-normalize-to-yeast-spike-in">Normalize to yeast spike-in</a>
    - <a href="#calculate-signal-at-peaks-area-under-histogram-curve-2"
      id="toc-calculate-signal-at-peaks-area-under-histogram-curve-2">Calculate
      Signal at Peaks: Area under Histogram Curve</a>
    - <a
      href="#plot-read-normalized-signal-vs-expected-signal-supplemental-figure"
      id="toc-plot-read-normalized-signal-vs-expected-signal-supplemental-figure">Plot
      Read normalized signal vs expected signal: Supplemental Figure</a>

``` r
library(tidyverse)
```

    Warning: package 'ggplot2' was built under R version 4.3.3

    Warning: package 'lubridate' was built under R version 4.3.2

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ✔ purrr     1.0.2     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(RColorBrewer)
theme_set(theme_classic())
library(DescTools)
```

    Warning: package 'DescTools' was built under R version 4.3.3

# EstimatingInterphase/Mitotic H3K9ac with Mass Spectrometry Data

``` r
javasky_H3K9ac_massspec <- read.delim("~/Research/spike_commentary/javasky_2018_massspec_H3K9ac_only.txt")
```

``` r
javasky_H3K9ac_massspec$relative_abundance <- as.numeric(javasky_H3K9ac_massspec$relative_abundance)
```

``` r
javasky_H3K9ac_massspec$Sample <- factor(javasky_H3K9ac_massspec$Sample, levels = c("H3K9ac_unsync", "H3K9ac_sync"))
```

``` r
ggplot(javasky_H3K9ac_massspec) + 
  aes(x = Sample, y = relative_abundance) + 
  geom_col() + 
  geom_errorbar(aes(ymin = relative_abundance- sd, ymax = relative_abundance + sd), width = 0.2) + 
  labs(title = "Relative Abundance of H3K9ac in Mitotic and Interphase HeLa-S3 Cells",
       subtitle = "Quantified with Mass-Spectrometry in Javasky et al 2018", 
       x = "Sample", 
       y = "Relative Abundance (%)")
```

![](spike_correspondence_figure1b_files/figure-commonmark/javasky_H3K9ac_massspec_plot_publish-1.png)

``` r
javasky_H3K9ac_massspec_plot <- ggplot(javasky_H3K9ac_massspec) + 
  aes(x = Sample, y = relative_abundance) + 
  geom_col() + 
  geom_errorbar(aes(ymin = relative_abundance- sd, ymax = relative_abundance + sd), width = 0.2) + 
  labs(title = "Relative Abundance of H3K9ac in Mitotic and Interphase HeLa-S3 Cells",
       subtitle = "Quantified with Mass-Spectrometry in Javasky et al 2018", 
       x = "Sample", 
       y = "Relative Abundance (%)")

ggsave("javasky_H3K9ac_massspec_plot.svg", javasky_H3K9ac_massspec_plot, width = 6, height = 5)
```

# ChIP-seq Titration of H3K9ac

Sequenced on NextSeq 550, demultiplexed with Illumina bcl2fastq

## Trimming .FASTQ Files

``` bash
base=$(basename $i _R1_001.fastq.gz)
trimmomatic SE ${base}_R1_001.fastq.gz ${dir}/${base}.trim.fastq.gz \
-threads 8 \
ILLUMINACLIP:/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:2:30:7 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
```

## Alignment: .FASTQ -\> .SAM

### Genome Preparation:

Genomes Downloaded with HOMER: hg38, dm6, sacCer3

``` bash
perl /path/to/homer/configureHomer.pl -install hg38

perl /path/to/homer/configureHomer.pl -install dm6

perl /path/to/homer/configureHomer.pl -install sacCer3
```

Add chromosome suffixes to identify spike-in chromosomes

Fasta file of genome location: /homer/data/genomes/hg38/genome.fa

``` bash
## Fly dm6 Genome
# suffix added "_dm6"
sed 's/>.*/&_dm6/' genome.fa > genome_dm6.fa

# check by printing fastq headers
perl -ne 'if(/^>(\S+)/){print "$1\n"}' genome_dm6.fa

## Yeast sacCer3 Genome
# suffix added "_sac3"
sed 's/>.*/&_sac3/' genome.fa > genome_sac3.fa

# check by printing fastq headers
perl -ne 'if(/^>(\S+)/){print "$1\n"}' genome_sac3.fa
```

Combine spike-in/target genomes:

``` bash
cat ${dir}/dm6/genome_dm6.fa ${dir}/sacCer3/genome_sac3.fa ${dir}/hg38/genome.fa > genome_hg38_dm6_sac3.fa
```

Index genome:

``` bash
bwa index -p hg38_dm6_sac3 genome_hg38_dm6_sac3.fa
```

Alignment:

``` bash
bwa mem -t 4 ~/data/genome_index/genome_prefix ${file}.fastq > ${file}.sam
```

## Remove PCR Duplicates

``` bash
# single end reads: don't need to collate
# fixmate can go from sam to bam 
samtools fixmate -m ${file}.sam ${file}.fixmate.bam
samtools sort ${file}.fixmate.bam -o ${file}.sorted.bam
samtools markdup -r -s ${file}.sorted.bam ${file}.nodup.bam
```

## Separate Alignment File into Species-specific Alignments

Split alignment file into files for each chromosome

``` bash
bam splitChromosome --in ${file}.nodup.bam --out ${file}.
```

Remove *chrUn*, *random* files

Merge chromosome files to get one alignment file per species. Need to
merge spike-in species first, then remove those chromosomes

``` bash
# create dm6 spike-in alignment file
samtools merge ${file}.dm6.bam ${file}.chr*_dm6.bam

# remove dm6 spike-in chromosome files
rm ${file}.chr*_dm6.bam

# create sac3 spike-in alignment file
samtools merge ${file}.sac3.bam ${file}.chr*_sac3.bam

# remove sac3 spike-in chromosome files
rm ${file}.chr*_sac3.bam

# create target alignment file
samtools merge ${file}.hg38.bam ${file}.chr*.bam
```

Convert BAM to SAM (-h to keep header), remove suffixes from files with
sed.

``` bash
# remove first set of spike-in suffixes
samtools view -h ${file}.bam | sed -e 's/\_dm6//g' > ${file}.nosuffix.sam

# remove second set of spike-in suffixes
samtools view -h ${file}.nosuffix.sam | sed -e 's/\_sac3//g' > ${file}.nosuffix2.sam
```

## Make HOMER Tag Directories

By default, HOMER `makeTagDirectory` or `batchMakeTagDirectory` only
keep primary alignments with MAPQ \> 10. Tag directories are also
read-depth normalized to 10 million reads unless otherwise specified. To
create Tag Directories in batch mode, first create a tsv file called a
tagkey, containing the names for each input file and each desired output
tag directory.

Format of the tagkey file: TSV <br> file1-tagdir file1.nosuffix2.sam
<br> file2-tagdir file2.nosuffix2.sam <br> file3-tagdir
file3.nosuffix2.sam <br>

Then run HOMER `batchMakeTagDirectory.pl` for each species. Example with
hg38:

``` bash
batchMakeTagDirectory.pl tagkey -genome hg38 \
-cpu 8 -fragLength 150
```

Repeat for dm6, sacCer3 genomes

## Visualize with BigWigs

Can make BigWigs from bam file with Deeptools, or from Tag Directories
with HOMER `makeBigWig.pl`

``` bash
makeBigWig.pl ${file}-tagdir/ hg38 \
-webdir /path/to/webdirectory -url http://webdirectoryurl/
```

Note: you can also make BedGraphs with HOMER `makeBedGraph.pl`

Repeat for dm6, sacCer3 genomes

## Determine Normalization Factor

## Analysis with HOMER

### Histograms at TSS

Recommended parameters for H3K9ac: Make histogram of size 4kb centered
at TSS, with bin size 25bp.

``` bash
annotatePeaks.pl tss hg38 -size 4000 -hist 25 \
-d ${file}1-tagdir ${file}2-tagdir \
> histogram_tss_hg38_samples.txt
```

### Peak Finding

Recommended parameters for human histone mark peak finding: -style
histone, -size 1000, -minDist 2500. Note: yeast genome is
gene/acetylation dense, so if peak finding in yeast is desired,
paramters should be tweaked.

``` bash
findPeaks.pl ${file}-tagdir -style histone -size 1000 -minDist 2500 \
-i ${file}-input-tagdir > ${file}.regions.txt
```

### Quantification at Peaks

``` bash
annotatePeaks.pl ${file}.regions.txt -size 1000 \
-d ${file1}-tagdir ${file2}-tagdir > counts_regions_1kb_samples.txt
```

## Megapeak: Get H3K9ac regions to quantify signal

1.  Combined all samples to one “mega-ChIP-seq sample”, and all inputs
    to one “mega-input control”.
2.  Did peak finding with “mega ChIP-seq sample” with “mega-input
    control” as background. The overall increased read-depth of sample
    and input files will increase the tag threshold for peak finding,
    resulting in more significant peaks.
3.  Quantified H3K9ac signal at these peaks for every sample, making a
    histogram of size 4000 and bin size 25.

``` bash
findPeaks.pl ${file}.megasample-tagdir -style histone -size 1000 -minDist 2500 \
-i ${file}.megainput-tagdir > megasample_H3K9ac.regions.txt

annotatePeaks.pl megasample_H3K9ac.regions.txt hg38 -size 4000 -hist 25 \
-d ${file1}-tagdir ${file2}-tagdir ${file3}-tagdir .... > hist_K9ac_hg38_allsamples_LP78.txt
```

# Import histogram data

``` r
hist_K9ac_allsamples_LP78 <- read.delim("~/Research/LP_78/hist_K9ac_hg38_allsamples_LP78.txt")
```

``` r
process_my_histograms <- function(x, .x) {
    colnames(x)[1] <- "Distance_from_center"
    x <- x %>% 
    rename_with(~ gsub(".hg38.tagdir", "", .x), contains("tagdir")) %>% 
    rename_with(~ gsub(".+\\_Hela_", "HelaS3_", .x), contains("Hela")) %>%
      rename_with(~ gsub("ratio1", "100sync_0inter", .x), contains("ratio1")) %>%
      rename_with(~ gsub("ratio2", "75sync_25inter", .x), contains("ratio2")) %>%
      rename_with(~ gsub("ratio3", "50sync_50inter", .x), contains("ratio3")) %>%
      rename_with(~ gsub("ratio4", "25sync_75inter", .x), contains("ratio4")) %>%
      rename_with(~ gsub("ratio5", "5sync_95inter", .x), contains("ratio5")) %>%
      rename_with(~ gsub("ratio6", "0sync_100inter", .x), contains("ratio6")) %>%
      rename_with(~ gsub("K9ac", "H3K9ac", .x), contains("K9ac")) %>%
      rename_with(~ gsub("ac_r", "ac_rep", .x), contains("K9ac")) %>%
    rename_with(~ gsub("\\.[[:digit:]]$", "_minus", .x), contains("Tags")) %>% 
    rename_with(~ gsub("\\.\\.\\.", "_", .x), contains("Tags"))
    
    xcov <- x %>% select(contains("Coverage"))
    xcov$Distance_from_center <- x$Distance_from_center
    
    xcovlong <- 
    xcov %>% pivot_longer(
      cols = -"Distance_from_center", 
      names_to = "Sample", 
      values_to = "Coverage")
  }
```

``` r
hist_K9ac_allsamples_LP78_tidy <- 
  process_my_histograms(hist_K9ac_allsamples_LP78)
```

``` r
hist_K9ac_allsamples_LP78_sepIP <- hist_K9ac_allsamples_LP78_tidy[grep("K9ac", hist_K9ac_allsamples_LP78_tidy$Sample), ] %>% separate_wider_regex(cols = Sample, patterns = c(
  cell = ".+", 
  "\\_", 
  ratio_sync = "[:digit:]+", 
  "sync", "\\_", 
  ratio_inter = "[:digit:]+", 
  "inter", "\\_", 
  antibody = ".+", 
  "\\_", 
  replicate = ".+", 
  ".Coverage"))
```

Fix ordering of mitotic and interphase cell ratios:

``` r
hist_K9ac_allsamples_LP78_sepIP$ratio_sync <- factor(hist_K9ac_allsamples_LP78_sepIP$ratio_sync, levels = c(100, 75, 50, 25, 5, 0))

hist_K9ac_allsamples_LP78_sepIP$ratio_inter <- factor(hist_K9ac_allsamples_LP78_sepIP$ratio_inter, levels = c(0, 25, 50, 75, 95, 100))
```

``` r
ggplot(data = hist_K9ac_allsamples_LP78_sepIP) + 
  aes(x = Distance_from_center, y = Coverage, group=interaction(ratio_inter, antibody, replicate), color = ratio_inter) + 
  geom_line(alpha = 0.9, linewidth = 1.1) + 
  labs(title = "Hela H3K9ac megasample peaks", 
       x = "Distance from Peak Center") +
  scale_color_manual(
    values = colorRampPalette(brewer.pal(9, "Blues"))(8)[3:8],
    name = "Cell Ratio", 
    labels = c("0% interphase", "25% interphase", "50% interphase", "75% interphase", "95% interphase", "100% interphase")) + 
  theme_classic() + 
  theme(legend.position = c(0.84, 0.76), 
                          legend.background = element_rect(
                                  size=0.7, linetype="solid", 
                                  colour ="grey20")) 
```

    Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ℹ Please use the `linewidth` argument instead.

    Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    3.5.0.
    ℹ Please use the `legend.position.inside` argument of `theme()` instead.

![](spike_correspondence_figure1b_files/figure-commonmark/H3K9ac_titration_readnorm_histplot-1.png)

### Calculate Signal at Peaks: Area under Histogram Curve

``` r
hist_K9ac_allsamples_LP78_tidyIP <- hist_K9ac_allsamples_LP78_tidy %>% filter(grepl("H3K9ac", Sample))

samples <- unique(hist_K9ac_allsamples_LP78_tidyIP$Sample)

colnames(hist_K9ac_allsamples_LP78)[1] <- "Distance_from_center"
x <- hist_K9ac_allsamples_LP78$Distance_from_center
AUC_peaks <- matrix(data = "", nrow = length(samples), ncol = 1)
AUC_peaks <- data.frame(AUC_peaks, row.names = samples)

for (i in 1:length(samples)) {

  y <- hist_K9ac_allsamples_LP78_tidy %>% 
    filter(Sample == samples[i]) %>%
    select(Coverage)
  y <- pull(y, Coverage)

AUC_peaks[i, ] <- AUC(x, y, method = c("trapezoid"))

}
```

### Make line of expected signal:

``` r
expected_line <- c(0, 0.25, 0.5, 0.75, 0.95, 1)

percent_inter_mean <- rep(c(0, 25, 50, 75, 95, 100))

percent_mit_mean <- rep(c(100, 75, 50, 25, 5, 0))
```

### scale signal, set minimum signal to 0, maximum signal to 1

``` r
scale_AUC_minmaxnorm <- function(AUC_peaks) {
  AUC_peaks$AUC_peaks <- as.numeric(AUC_peaks$AUC_peaks)
  
  interphase <- rep(c(0, 25, 50, 75, 95, 100), each = 3)
  AUC_peaks$interphase <- as.numeric(interphase)
  
  # Minmax normalization: Scale points from 0-1
  avg_0inter <- mean(c(AUC_peaks[1,1], AUC_peaks[2,1], AUC_peaks[3,1]))
  avg_25inter <- mean(c(AUC_peaks[4,1], AUC_peaks[5,1], AUC_peaks[6,1]))
  avg_50inter <- mean(c(AUC_peaks[7,1], AUC_peaks[8,1], AUC_peaks[9,1]))
  avg_75inter <- mean(c(AUC_peaks[10,1], AUC_peaks[11,1], AUC_peaks[12,1]))
  avg_95inter <- mean(c(AUC_peaks[13,1], AUC_peaks[14,1], AUC_peaks[15,1]))
  avg_100inter <- mean(c(AUC_peaks[16,1], AUC_peaks[17,1], AUC_peaks[18,1]))
  
  AUC_peaks_minmaxnorm <- AUC_peaks %>%
    mutate(minmaxnorm = (AUC_peaks-avg_0inter)/(avg_100inter-avg_0inter) )
}
```

``` r
AUC_peaks_readnorm_minmaxnorm <- scale_AUC_minmaxnorm(AUC_peaks)
```

### Make line of expected signal:

``` r
# make variable for % mitotic cells
AUC_peaks_readnorm_minmaxnorm$mitotic <- rep(percent_mit_mean, each = 3)

# dataframe of expected line
observed_vs_expected_line <- data.frame(cbind(percent_mit_mean, expected_line))
```

### Plot Read normalized signal vs expected signal: Fig 1b

``` r
ggplot() +
  geom_point(data = AUC_peaks_readnorm_minmaxnorm, 
             aes(x = as.numeric(mitotic), y = minmaxnorm), 
             size = 2, alpha = 0.7, shape = 3, color = "grey30", stroke = 3) +
  scale_color_manual(name = "Cell Ratio", 
    labels = c("0% interphase", "25% interphase", "50% interphase", "75% interphase", "95% interphase", "100% interphase")) +
  geom_line(data = observed_vs_expected_line, 
            aes(x = as.numeric(percent_mit_mean), y = expected_line), linewidth = 1.1, color = "grey30") +
 labs(title = "Relative H3K9ac Signal", 
       subtitle = "Read Normalized, scaled 0-1",
       x = "% Mitotic Cells", 
       y = "Relative H3K9ac Signal") + 
  theme(legend.position = "none")
```

![](spike_correspondence_figure1b_files/figure-commonmark/H3K9ac_titration_readnorm_publish_dotplot-1.png)

Calculated Rsquared

``` r
get_Rsquared <- function(AUC_peaks) {
  AUC_peaks$expected <- rep(c(0, 0.25, 0.5, 0.75, 0.95, 1), each = 3)
  
  AUC_peaks$error <- AUC_peaks$minmaxnorm - AUC_peaks$expected
  
  mse = sum(((AUC_peaks$error)^2)/3)
  
  SSres = sum((AUC_peaks$error)^2)
  meanobserved <- mean(AUC_peaks$minmaxnorm)
  SStotal = sum((AUC_peaks$minmaxnorm-meanobserved)^2)
  rsquared = 1 - (SSres/SStotal)
  rsquared
}
```

``` r
get_Rsquared(AUC_peaks_readnorm_minmaxnorm)
```

    [1] 0.06828471

# Normalize peak signal to Dual Spike-ins

### Determine normalization factor from aligned reads:

``` r
H3K9ac_mitotic_titration_seqtats <- read.delim("~/Research/spike_commentary/mydata_sequencing_statistics_mitotic_titration.txt")
```

``` r
H3K9ac_mitotic_titration_seqtatsIP <- H3K9ac_mitotic_titration_seqtats %>%
  filter(type == "ip")
```

``` r
knitr::kable(H3K9ac_mitotic_titration_seqtatsIP)
```

| Sample.ID                         | type | X..mitotic.cells | X..interphase.cells |    reads | dm6.aligned | sac3.aligned | hg38.aligned | dm6.sac3 | dm6.hg38 | sac3.hg38 | hg38.peaks | hg38.FRIP | fly.ip.input.norm | yeast.ip.input.norm | dual.ip.input.norm | X   |
|:----------------------------------|:-----|-----------------:|--------------------:|---------:|------------:|-------------:|-------------:|---------:|---------:|----------:|-----------:|----------:|------------------:|--------------------:|-------------------:|:----|
| HelaS3_100sync_0inter_H3K9ac_rep1 | ip   |              100 |                   0 | 19163813 |     1062178 |      5045695 |     11688694 |   0.2105 |  0.09087 |   0.43167 |      26415 |     21.30 |            7.8269 |             55.2714 |            55.0298 | NA  |
| HelaS3_100sync_0inter_H3K9ac_rep2 | ip   |              100 |                   0 | 17838260 |      958103 |      4624313 |     10951006 |   0.2072 |  0.08749 |   0.42227 |      24495 |     20.82 |            7.5357 |             54.0679 |            53.4090 | NA  |
| HelaS3_100sync_0inter_H3K9ac_rep3 | ip   |              100 |                   0 | 17185695 |     1085502 |      5355943 |      9706357 |   0.2027 |  0.11183 |   0.55180 |      29846 |     28.75 |            9.6322 |             70.6530 |            69.0393 | NA  |
| HelaS3_75sync_25inter_H3K9ac_rep1 | ip   |               75 |                  25 | 17095802 |      637243 |      3312364 |     12182780 |   0.1924 |  0.05231 |   0.27189 |      26658 |     24.47 |            4.7339 |             35.6811 |            34.4093 | NA  |
| HelaS3_75sync_25inter_H3K9ac_rep2 | ip   |               75 |                  25 | 20167562 |      670064 |      3201780 |     15143202 |   0.2093 |  0.04425 |   0.21143 |      23593 |     20.24 |            4.0045 |             27.7467 |            27.8892 | NA  |
| HelaS3_75sync_25inter_H3K9ac_rep3 | ip   |               75 |                  25 | 15777520 |      587962 |      3102874 |     11242290 |   0.1895 |  0.05230 |   0.27600 |      27700 |     25.97 |            4.7330 |             36.2205 |            34.6758 | NA  |
| HelaS3_50sync_50inter_H3K9ac_rep1 | ip   |               50 |                  50 | 20417995 |      564595 |      3329126 |     15334517 |   0.1696 |  0.03682 |   0.21710 |      27960 |     26.82 |            3.3112 |             29.4973 |            26.3377 | NA  |
| HelaS3_50sync_50inter_H3K9ac_rep2 | ip   |               50 |                  50 | 24102172 |      680504 |      3627376 |     18461307 |   0.1876 |  0.03686 |   0.19649 |      28549 |     26.55 |            3.3147 |             26.6970 |            24.9501 | NA  |
| HelaS3_50sync_50inter_H3K9ac_rep3 | ip   |               50 |                  50 | 17272978 |      472306 |      2262079 |     13658282 |   0.2088 |  0.03458 |   0.16562 |      24498 |     22.25 |            3.1097 |             22.5027 |            22.1354 | NA  |
| HelaS3_25sync_75inter_H3K9ac_rep1 | ip   |               25 |                  75 | 26192648 |      568915 |      3115984 |     20997543 |   0.1826 |  0.02709 |   0.14840 |      27308 |     24.98 |            2.6611 |             19.1484 |            18.8880 | NA  |
| HelaS3_25sync_75inter_H3K9ac_rep2 | ip   |               25 |                  75 | 23095083 |      525579 |      2815467 |     18448518 |   0.1867 |  0.02849 |   0.15261 |      27884 |     27.24 |            2.7986 |             19.6916 |            19.6410 | NA  |
| HelaS3_25sync_75inter_H3K9ac_rep3 | ip   |               25 |                  75 | 21971781 |      429002 |      2096816 |     18136505 |   0.2046 |  0.02365 |   0.11561 |      21852 |     18.95 |            2.3232 |             14.9174 |            15.5898 | NA  |
| HelaS3_5sync_95inter_H3K9ac_rep1  | ip   |                5 |                  95 | 13961566 |      291234 |      1569811 |     11304907 |   0.1855 |  0.02576 |   0.13886 |      27084 |     25.93 |            2.5230 |             19.4755 |            18.5683 | NA  |
| HelaS3_5sync_95inter_H3K9ac_rep2  | ip   |                5 |                  95 | 23002199 |      413781 |      2036320 |     19071919 |   0.2032 |  0.02170 |   0.10677 |      24418 |     22.35 |            2.1254 |             14.9748 |            14.9262 | NA  |
| HelaS3_5sync_95inter_H3K9ac_rep3  | ip   |                5 |                  95 | 20406280 |      397349 |      2202504 |     16560318 |   0.1804 |  0.02399 |   0.13300 |      27421 |     27.75 |            2.3497 |             18.6536 |            17.5506 | NA  |
| HelaS3_0sync_100inter_H3K9ac_rep1 | ip   |                0 |                 100 | 24461245 |      388596 |      1852990 |     20952226 |   0.2097 |  0.01855 |   0.08844 |      21921 |     18.22 |            1.8550 |             11.8235 |            12.4043 | NA  |
| HelaS3_0sync_100inter_H3K9ac_rep2 | ip   |                0 |                 100 | 21247480 |      414971 |      2080661 |     17513474 |   0.1994 |  0.02369 |   0.11880 |      27893 |     28.69 |            2.3690 |             15.8824 |            16.2327 | NA  |
| HelaS3_0sync_100inter_H3K9ac_rep3 | ip   |                0 |                 100 | 20656642 |      396011 |      2143990 |     17006383 |   0.1847 |  0.02329 |   0.12607 |      27679 |     26.81 |            2.3290 |             16.8543 |            16.5786 | NA  |

## Normalize to fly IP/input

``` r
hist_K9ac_allsamples_LP78_tidyIP <- hist_K9ac_allsamples_LP78_tidy %>% 
  filter(grepl("H3K9ac", Sample))

# copy the hist_tss_hg38_LH58_cov dataframe 
peakcov_fly_ip_input_norm <- hist_K9ac_allsamples_LP78_tidyIP

sampleID <- unique(hist_K9ac_allsamples_LP78_tidyIP$Sample)

seqstatID <- sub('.Coverage', "", sampleID )

# dataframe 1: hist_K9ac_allsamples_LP78_tidy (original read-normalized data) 
# dataframe 2: H3K9ac_mitotic_titration_seqtatsIP
# dataframe 3: peakcov_avg_ip_input_norm (output df)

# When Sample rows of df1 match Sample.ID in df2, multiply Coverage column in df1 by factor in df2, assign to df3

for (i in 1:nrow(hist_K9ac_allsamples_LP78_tidyIP)) {
  if (!hist_K9ac_allsamples_LP78_tidyIP[i, 2] %in% paste0(H3K9ac_mitotic_titration_seqtatsIP$Sample.ID, ".Coverage")) {
    next()
  }
  
  # make get current sampleID, remove .Coverage 
 seqstatIDi <- sub('.Coverage', "", hist_K9ac_allsamples_LP78_tidyIP[i, 2] )
 
 # get normalization factor from sequencing stats (df3)
 # Col 14 in seqstats dataframe contains fly normalization factor
 
 normfactori <- H3K9ac_mitotic_titration_seqtatsIP[grep(seqstatIDi, H3K9ac_mitotic_titration_seqtatsIP$Sample.ID), 14]
 
 # multiply read_norm coverage by norm factor, assign to new df
peakcov_fly_ip_input_norm[i, 3] <- 
  hist_K9ac_allsamples_LP78_tidyIP[i, 3]/(normfactori)
  
}
```

``` r
peakcov_fly_ip_input_norm_sepIP <- peakcov_fly_ip_input_norm[grep("K9ac", peakcov_fly_ip_input_norm$Sample), ] %>% separate_wider_regex(cols = Sample, patterns = c(
  cell = ".+", 
  "\\_", 
  ratio_sync = "[:digit:]+", 
  "sync", "\\_", 
  ratio_inter = "[:digit:]+", 
  "inter", "\\_", 
  antibody = ".+", 
  "\\_", 
  replicate = ".+", 
  ".Coverage"))
```

Fix ordering of mitotic and interphase cell ratios:

``` r
peakcov_fly_ip_input_norm_sepIP$ratio_sync <- factor(peakcov_fly_ip_input_norm_sepIP$ratio_sync, levels = c(100, 75, 50, 25, 5, 0))

peakcov_fly_ip_input_norm_sepIP$ratio_inter <- factor(peakcov_fly_ip_input_norm_sepIP$ratio_inter, levels = c(0, 25, 50, 75, 95, 100))
```

``` r
ggplot(data = peakcov_fly_ip_input_norm_sepIP) + 
  aes(x = Distance_from_center, y = Coverage, group=interaction(ratio_inter, antibody, replicate), color = ratio_inter) + 
  geom_line(alpha = 0.9, linewidth = 1.1) + 
  labs(title = "Hela H3K9ac megasample peaks",
       subtitle = "Normalized to Fly Spike-in",
       x = "Distance from Peak Center") +
  scale_color_manual(
    values = colorRampPalette(brewer.pal(9, "Blues"))(8)[3:8],
    name = "Cell Ratio", 
    labels = c("0% interphase", "25% interphase", "50% interphase", "75% interphase", "95% interphase", "100% interphase")) + 
  theme_classic() + 
  theme(legend.position = c(0.84, 0.76), 
                          legend.background = element_rect(
                                  size=0.7, linetype="solid", 
                                  colour ="grey20")) 
```

![](spike_correspondence_figure1b_files/figure-commonmark/H3K9ac_titration_spikenorm_histplot-1.png)

### Calculate Signal at Peaks: Area under Histogram Curve

``` r
samples <- unique(peakcov_fly_ip_input_norm$Sample)

colnames(hist_K9ac_allsamples_LP78)[1] <- "Distance_from_center"
x <- hist_K9ac_allsamples_LP78$Distance_from_center

AUC_peaks <- matrix(data = "", nrow = length(samples), ncol = 1)
AUC_peaks <- data.frame(AUC_peaks, row.names = samples)

for (i in 1:length(samples)) {

  y <- peakcov_fly_ip_input_norm %>% 
    filter(Sample == samples[i]) %>%
    select(Coverage)
  y <- pull(y, Coverage)

AUC_peaks[i, ] <- AUC(x, y, method = c("trapezoid"))

}
```

### Scale signal at peaks to be from 0-1, plot this normalized signal against expected signal

``` r
AUC_peaks_fly_minmaxnorm <- scale_AUC_minmaxnorm(AUC_peaks)

AUC_peaks_fly_minmaxnorm$mitotic <- rep(percent_mit_mean, each = 3)
```

``` r
#observed_vs_expected_LP78 <- data.frame(cbind(percent_inter_mean, percent_mit_mean, expected_line))
```

### Plot Read normalized signal vs expected signal: Fig 1b

``` r
ggplot() +
  geom_point(data = AUC_peaks_fly_minmaxnorm, 
             aes(x = as.numeric(mitotic), y = minmaxnorm), 
             size = 2, alpha = 0.7, shape = 3, color = "grey30", stroke = 3) +
  scale_color_manual(name = "Cell Ratio", 
    labels = c("0% interphase", "25% interphase", "50% interphase", "75% interphase", "95% interphase", "100% interphase")) +
  geom_line(data = observed_vs_expected_line, 
            aes(x = as.numeric(percent_mit_mean), y = expected_line), linewidth = 1.1, color = "grey30") +
 labs(title = "Relative H3K9ac Signal", 
       subtitle = "Normalized to Drosophila spike-in, scaled 0-1",
       x = "% Mitotic Cells", 
       y = "Relative H3K9ac Signal") + 
  theme(legend.position = "none")
```

![](spike_correspondence_figure1b_files/figure-commonmark/H3K9ac_titration_flynorm_publish_dotplot-1.png)

Rsquared value:

``` r
get_Rsquared(AUC_peaks_fly_minmaxnorm)
```

    [1] 0.977569

Save plot:

``` r
H3K9ac_titration_flynorm_publish_dotplot <- ggplot() +
  geom_point(data = AUC_peaks_fly_minmaxnorm, 
             aes(x = as.numeric(mitotic), y = minmaxnorm), 
             size = 2, alpha = 0.7, shape = 3, color = "grey30", stroke = 3) +
  scale_color_manual(name = "Cell Ratio", 
    labels = c("0% interphase", "25% interphase", "50% interphase", "75% interphase", "95% interphase", "100% interphase")) +
  geom_line(data = observed_vs_expected_line, 
            aes(x = as.numeric(percent_mit_mean), y = expected_line), linewidth = 1.1, color = "grey30") +
 labs(title = "Relative H3K9ac Signal", 
       subtitle = "Normalized to Drosophila spike-in, scaled 0-1",
       x = "% Mitotic Cells", 
       y = "Relative H3K9ac Signal") + 
  theme(legend.position = "none")

ggsave("H3K9ac_titration_flynorm_publish_dotplot.svg", H3K9ac_titration_flynorm_publish_dotplot, width = 4.5, height = 4.5)
```

## Normalize to yeast spike-in

``` r
hist_K9ac_allsamples_LP78_tidyIP <- hist_K9ac_allsamples_LP78_tidy %>% 
  filter(grepl("H3K9ac", Sample))

# copy the hist_tss_hg38_LH58_cov dataframe 
peakcov_yeast_ip_input_norm <- hist_K9ac_allsamples_LP78_tidyIP

sampleID <- unique(hist_K9ac_allsamples_LP78_tidyIP$Sample)

seqstatID <- sub('.Coverage', "", sampleID )

# dataframe 1: hist_K9ac_allsamples_LP78_tidy (original read-normalized data) 
# dataframe 2: H3K9ac_mitotic_titration_seqtatsIP
# dataframe 3: peakcov_avg_ip_input_norm (output df)

# When Sample rows of df1 match Sample.ID in df2, multiply Coverage column in df1 by factor in df2, assign to df3

for (i in 1:nrow(hist_K9ac_allsamples_LP78_tidyIP)) {
  if (!hist_K9ac_allsamples_LP78_tidyIP[i, 2] %in% paste0(H3K9ac_mitotic_titration_seqtatsIP$Sample.ID, ".Coverage")) {
    next()
  }
  
  # make get current sampleID, remove .Coverage 
 seqstatIDi <- sub('.Coverage', "", hist_K9ac_allsamples_LP78_tidyIP[i, 2] )
 
 # get normalization factor from sequencing stats (df3)
 # Col 15 in seqstats dataframe contains yeast normalization factor

 normfactori <- H3K9ac_mitotic_titration_seqtatsIP[grep(seqstatIDi, H3K9ac_mitotic_titration_seqtatsIP$Sample.ID), 15]
 
 # multiply read_norm coverage by norm factor, assign to new df
peakcov_yeast_ip_input_norm[i, 3] <- 
  hist_K9ac_allsamples_LP78_tidyIP[i, 3]/(normfactori)
  
}
```

``` r
peakcov_yeast_ip_input_norm_sepIP <- peakcov_yeast_ip_input_norm[grep("K9ac", peakcov_yeast_ip_input_norm$Sample), ] %>% separate_wider_regex(cols = Sample, patterns = c(
  cell = ".+", 
  "\\_", 
  ratio_sync = "[:digit:]+", 
  "sync", "\\_", 
  ratio_inter = "[:digit:]+", 
  "inter", "\\_", 
  antibody = ".+", 
  "\\_", 
  replicate = ".+", 
  ".Coverage"))
```

Fix ordering of mitotic and interphase cell ratios:

``` r
peakcov_yeast_ip_input_norm_sepIP$ratio_sync <- factor(peakcov_yeast_ip_input_norm_sepIP$ratio_sync, levels = c(100, 75, 50, 25, 5, 0))

peakcov_yeast_ip_input_norm_sepIP$ratio_inter <- factor(peakcov_yeast_ip_input_norm_sepIP$ratio_inter, levels = c(0, 25, 50, 75, 95, 100))
```

``` r
ggplot(data = peakcov_yeast_ip_input_norm_sepIP) + 
  aes(x = Distance_from_center, y = Coverage, group=interaction(ratio_inter, antibody, replicate), color = ratio_inter) + 
  geom_line(alpha = 0.9, linewidth = 1.1) + 
  labs(title = "Hela H3K9ac megasample peaks",
       subtitle = "Normalized to Yeast Spike-in",
       x = "Distance from Peak Center") +
  scale_color_manual(
    values = colorRampPalette(brewer.pal(9, "Blues"))(8)[3:8],
    name = "Cell Ratio", 
    labels = c("0% interphase", "25% interphase", "50% interphase", "75% interphase", "95% interphase", "100% interphase")) + 
  theme_classic() + 
  theme(legend.position = c(0.84, 0.76), 
                          legend.background = element_rect(
                                  size=0.7, linetype="solid", 
                                  colour ="grey20")) 
```

![](spike_correspondence_figure1b_files/figure-commonmark/H3K9ac_titration_yeastnorm_histplot-1.png)

### Calculate Signal at Peaks: Area under Histogram Curve

``` r
samples <- unique(peakcov_yeast_ip_input_norm$Sample)

colnames(hist_K9ac_allsamples_LP78)[1] <- "Distance_from_center"
x <- hist_K9ac_allsamples_LP78$Distance_from_center

AUC_peaks <- matrix(data = "", nrow = length(samples), ncol = 1)
AUC_peaks <- data.frame(AUC_peaks, row.names = samples)

for (i in 1:length(samples)) {

  y <- peakcov_yeast_ip_input_norm %>% 
    filter(Sample == samples[i]) %>%
    select(Coverage)
  y <- pull(y, Coverage)

AUC_peaks[i, ] <- AUC(x, y, method = c("trapezoid"))

}
```

Processing Signal at peaks, add expected line

``` r
AUC_peaks_yeast_minmaxnorm <- scale_AUC_minmaxnorm(AUC_peaks)

AUC_peaks_yeast_minmaxnorm$mitotic <- rep(percent_mit_mean, each = 3)
```

### Plot Read normalized signal vs expected signal: Supplemental Figure

``` r
ggplot() +
  geom_point(data = AUC_peaks_yeast_minmaxnorm, 
             aes(x = as.numeric(mitotic), y = minmaxnorm), 
             size = 2, alpha = 0.7, shape = 3, color = "grey30", stroke = 3) +
  scale_color_manual(name = "Cell Ratio", 
    labels = c("0% interphase", "25% interphase", "50% interphase", "75% interphase", "95% interphase", "100% interphase")) +
  geom_line(data = observed_vs_expected_line, 
            aes(x = as.numeric(percent_mit_mean), y = expected_line), linewidth = 1.1, color = "grey30") +
 labs(title = "Relative H3K9ac Signal", 
       subtitle = "Normalized to Yeast spike-in, scaled 0-1",
       x = "% Mitotic Cells", 
       y = "Relative H3K9ac Signal") + 
  theme(legend.position = "none")
```

![](spike_correspondence_figure1b_files/figure-commonmark/H3K9ac_titration_yeastnorm_dotplot-1.png)

``` r
get_Rsquared(AUC_peaks_yeast_minmaxnorm)
```

    [1] 0.954441
