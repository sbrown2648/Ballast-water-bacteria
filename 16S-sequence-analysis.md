Ballast water bacteria 16S rRNA analysis
================
Sarah Brown
July 19, 2023

-   <a href="#1-load-packages-and-import-data"
    id="toc-1-load-packages-and-import-data">1 Load Packages and Import
    Data</a>
-   <a href="#2-preparing-the-data" id="toc-2-preparing-the-data">2
    Preparing the Data</a>
    -   <a href="#21-examine-the-total-number-of-reads-and-asvs"
        id="toc-21-examine-the-total-number-of-reads-and-asvs">2.1 Examine the
        total number of reads and ASV’s</a>
    -   <a href="#22-create-table-showing-the-number-of-reads-per-sample"
        id="toc-22-create-table-showing-the-number-of-reads-per-sample">2.2
        Create table showing the number of reads per sample</a>
    -   <a href="#23-create-a-rarefaction-curve"
        id="toc-23-create-a-rarefaction-curve">2.3 Create a rarefaction
        curve</a>
    -   <a href="#24-filter-samples" id="toc-24-filter-samples">2.4 Filter
        samples</a>
    -   <a href="#25-label-and-order-variables"
        id="toc-25-label-and-order-variables">2.5 Label and order variables</a>
    -   <a href="#26-abundance-transformation"
        id="toc-26-abundance-transformation">2.6 Abundance transformation</a>
-   <a href="#3-alpha-diversity" id="toc-3-alpha-diversity">3 Alpha
    Diversity</a>
    -   <a href="#31-calculate-shannon-diversity-per-sample"
        id="toc-31-calculate-shannon-diversity-per-sample">3.1 Calculate Shannon
        Diversity per Sample</a>
    -   <a href="#32-create-alpha-diversity-boxplots"
        id="toc-32-create-alpha-diversity-boxplots">3.2 Create Alpha Diversity
        Boxplots</a>
        -   <a href="#321-comparison-sample-type"
            id="toc-321-comparison-sample-type">3.2.1 Comparison: sample type</a>
        -   <a href="#322-comparison-sample-type-and-voyage"
            id="toc-322-comparison-sample-type-and-voyage">3.2.2 Comparison: sample
            type and voyage</a>
    -   <a href="#33-alpha-diversity-statistics"
        id="toc-33-alpha-diversity-statistics">3.3 Alpha Diversity
        statistics</a>
-   <a href="#4-ordination-plots" id="toc-4-ordination-plots">4 Ordination
    plots</a>
    -   <a href="#41-find-plot-color-hex-codes"
        id="toc-41-find-plot-color-hex-codes">4.1 Find plot color hex codes</a>
    -   <a href="#42-create-nmds-with-the-bray-curtis-dissimilarity-matrix"
        id="toc-42-create-nmds-with-the-bray-curtis-dissimilarity-matrix">4.2
        Create NMDS with the Bray-Curtis dissimilarity matrix</a>
    -   <a href="#43-create-pcoa-with-the-bray-curtis-dissimilarity-matrix"
        id="toc-43-create-pcoa-with-the-bray-curtis-dissimilarity-matrix">4.3
        Create PCoA with the Bray-Curtis dissimilarity matrix</a>
        -   <a href="#431" id="toc-431">4.3.1</a>
-   <a href="#5-preliminary-taxonomy-fix"
    id="toc-5-preliminary-taxonomy-fix">5 Preliminary taxonomy fix</a>
-   <a href="#6-bacterial-community-composition-bar-charts"
    id="toc-6-bacterial-community-composition-bar-charts">6 Bacterial
    community composition bar charts</a>
    -   <a href="#61-find-color-hex-codes" id="toc-61-find-color-hex-codes">6.1
        Find color hex codes</a>
    -   <a href="#62-family-level-individual-samples"
        id="toc-62-family-level-individual-samples">6.2 Family level, individual
        samples</a>
    -   <a href="#63-family-level-grouped-by-sample-type"
        id="toc-63-family-level-grouped-by-sample-type">6.3 Family level,
        grouped by sample type</a>
    -   <a href="#64-family-level-grouped-by-sample-type-and-voyage"
        id="toc-64-family-level-grouped-by-sample-type-and-voyage">6.4 Family
        level, grouped by sample type and voyage</a>
-   <a href="#7-deseq2-analysis-on-metacyc-data"
    id="toc-7-deseq2-analysis-on-metacyc-data">7 DESeq2 analysis on MetaCyc
    data</a>
    -   <a href="#71-create-the-deseq-object"
        id="toc-71-create-the-deseq-object">7.1 Create the DESeq object</a>
    -   <a href="#72-load-pathway-descriptions"
        id="toc-72-load-pathway-descriptions">7.2 Load pathway descriptions</a>
    -   <a href="#73-create-pca-of-all-metacyc-results-in-deseq2"
        id="toc-73-create-pca-of-all-metacyc-results-in-deseq2">7.3 Create PCA
        of all MetaCyc results in DESeq2</a>
    -   <a href="#74-create-pcoa-of-all-metacyc-results"
        id="toc-74-create-pcoa-of-all-metacyc-results">7.4 Create PCoA of all
        MetaCyc results</a>
    -   <a href="#75-bwt-vs-port-uptake" id="toc-75-bwt-vs-port-uptake">7.5 BWT
        vs port uptake</a>
        -   <a href="#751-create-dot-plot-of-results---log2fold-change"
            id="toc-751-create-dot-plot-of-results---log2fold-change">7.5.1 Create
            dot plot of results - log2fold change</a>
-   <a href="#8-session-info" id="toc-8-session-info">8 Session Info</a>

# 1 Load Packages and Import Data

Here, we load the required packages and then import the .qza files
exported from QIIME2 using the qiime2R package. The final object is a
phyloseq object called Ballast_physeq.

``` r
library(ggplot2)
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

``` r
library(plyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(picante)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:reshape2':
    ## 
    ##     smiths

``` r
library(viridis)
```

    ## Loading required package: viridisLite

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:scales':
    ## 
    ##     viridis_pal

``` r
library(qiime2R)
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:plyr':
    ## 
    ##     rename

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:nlme':
    ## 
    ##     collapse

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:plyr':
    ## 
    ##     desc

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## The following object is masked from 'package:plyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

``` r
library(patchwork)
library(RColorBrewer)
library(RColorBrewer)
library(microViz)
```

    ## microViz version 0.10.8 - Copyright (C) 2022 David Barnett
    ## ! Website: https://david-barnett.github.io/microViz
    ## ✔ Useful?  For citation details, run: `citation("microViz")`
    ## ✖ Silence? `suppressPackageStartupMessages(library(microViz))`

``` r
library(speedyseq)
```

    ## 
    ## Attaching package: 'speedyseq'

    ## The following objects are masked from 'package:phyloseq':
    ## 
    ##     filter_taxa, plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom,
    ##     tip_glom, transform_sample_counts

``` r
library(ComplexHeatmap)
```

    ## ========================================
    ## ComplexHeatmap version 2.14.0
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite either one:
    ## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    ## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##     genomic data. Bioinformatics 2016.
    ## 
    ## 
    ## The new InteractiveComplexHeatmap package can directly export static 
    ## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

``` r
library(plotly)
```

    ## 
    ## Attaching package: 'plotly'

    ## The following object is masked from 'package:ComplexHeatmap':
    ## 
    ##     add_heatmap

    ## The following object is masked from 'package:microViz':
    ## 
    ##     add_paths

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, mutate, rename, summarise

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     last_plot

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

``` r
library(seecolor)
library(tibble)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ readr     2.1.4
    ## ✔ lubridate 1.9.2     ✔ stringr   1.5.0
    ## ✔ purrr     1.0.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%()    masks IRanges::%within%()
    ## ✖ plotly::arrange()        masks dplyr::arrange(), plyr::arrange()
    ## ✖ readr::col_factor()      masks scales::col_factor()
    ## ✖ IRanges::collapse()      masks nlme::collapse(), dplyr::collapse()
    ## ✖ Biobase::combine()       masks BiocGenerics::combine(), dplyr::combine()
    ## ✖ purrr::compact()         masks plyr::compact()
    ## ✖ matrixStats::count()     masks dplyr::count(), plyr::count()
    ## ✖ IRanges::desc()          masks dplyr::desc(), plyr::desc()
    ## ✖ purrr::discard()         masks scales::discard()
    ## ✖ S4Vectors::expand()      masks tidyr::expand()
    ## ✖ dplyr::failwith()        masks plyr::failwith()
    ## ✖ plotly::filter()         masks dplyr::filter(), stats::filter()
    ## ✖ S4Vectors::first()       masks dplyr::first()
    ## ✖ dplyr::id()              masks plyr::id()
    ## ✖ dplyr::lag()             masks stats::lag()
    ## ✖ plotly::mutate()         masks dplyr::mutate(), plyr::mutate()
    ## ✖ BiocGenerics::Position() masks ggplot2::Position(), base::Position()
    ## ✖ purrr::reduce()          masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ plotly::rename()         masks S4Vectors::rename(), dplyr::rename(), plyr::rename()
    ## ✖ lubridate::second()      masks S4Vectors::second()
    ## ✖ lubridate::second<-()    masks S4Vectors::second<-()
    ## ✖ plotly::slice()          masks IRanges::slice(), dplyr::slice()
    ## ✖ plotly::summarise()      masks dplyr::summarise(), plyr::summarise()
    ## ✖ dplyr::summarize()       masks plyr::summarize()
    ## ✖ ape::where()             masks dplyr::where()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
set.seed(13745)

#Create a phyloseq object from the .qza files exported from qiime2 using 
#the qiime2R package
Ballast_physeq <- qza_to_phyloseq(
  features="Qiime output/ballast-dada2-table.qza",
  tree="Qiime output/ballast-rooted-tree.qza",
  taxonomy="Qiime output/ballast-taxonomy.qza",
  metadata = "Qiime output/Metadata.tsv"
)
```

# 2 Preparing the Data

Before performing any analyses, we need to examine the phyloseq object
and remove chloroplast and mitochondrial sequences:

``` r
#Check the rank names to make sure they are accurate
rank_names(Ballast_physeq)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

``` r
#Correct output:[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  
#                   "Genus"   "Species"

#Check sample variables
sample_variables(Ballast_physeq)
```

    ## [1] "Sample_type"   "Voyage"        "Tank"          "Sample_number"
    ## [5] "Sample_site"   "Sample_date"   "notes"

``` r
#Remove chloroplast sequences and any contaminant sequences
Ballast_physeq <- subset_taxa(Ballast_physeq, Kingdom != "d__Archaea")
Ballast_physeq <- subset_taxa(Ballast_physeq, Order != "Chloroplast")
Ballast_physeq <- subset_taxa(Ballast_physeq, Family != "Mitochondria")

#Check that contaminant sequences are removed (easiest to save as data frame and search)
taxtabl <- as.data.frame(tax_table(Ballast_physeq))

#At this point, blanks and positive controls should also be removed
Ballast_physeq = subset_samples(Ballast_physeq, Sample_type != "control")
Ballast_physeq = subset_samples(Ballast_physeq, Sample_type != "negative")
```

## 2.1 Examine the total number of reads and ASV’s

This part of the code examines the number of times each ASV occurs in
the data set and the number of reads per sample; this will help
determine which ASV’s and samples need to be filtered out of the data
set.

``` r
#Create graph showing number of reads per ASV and number of reads per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(Ballast_physeq), TRUE), 
                        sorted = 1:ntaxa(Ballast_physeq), type = "ASVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(Ballast_physeq), TRUE), sorted = 1:nsamples(Ballast_physeq), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
```

![](16S-sequence-analysis_files/figure-gfm/asv-1.png)<!-- -->

``` r
#Create a dataframe that includes only the total number of times that each ASV occurs
readsumsdf_ASVs <- readsumsdf %>% 
    filter(type == "ASVs")  
  
#How many singletons are present?
length(which(readsumsdf_ASVs$nreads <= 0))
```

    ## [1] 11

``` r
#Answer: 11
length(which(readsumsdf_ASVs$nreads == 1))
```

    ## [1] 0

``` r
#Answer: 0

#How many doubletons are present?
length(which(readsumsdf_ASVs$nreads == 2))
```

    ## [1] 29

``` r
#Answer: 29

#Create a bar plot showing the number of times each ASV occurs
ASV_occurance <- ggplot(data = readsumsdf_ASVs) +
    geom_bar(mapping = aes(x=nreads)) + 
    ylim(0, 20) +
    labs(x = "Number of reads", y = "Number of ASV's") +
    theme(text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank())
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(ASV_occurance)
```

    ## Warning: Removed 2 rows containing missing values (`geom_bar()`).

![](16S-sequence-analysis_files/figure-gfm/ASV%20abundance-1.png)<!-- -->

## 2.2 Create table showing the number of reads per sample

``` r
#Sum reads and sort by samples with least to greatest number of reads
samples_reads <- sort(sample_sums(Ballast_physeq))

#Create table 
knitr::kable(samples_reads) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

    ## Warning in !is.null(rmarkdown::metadata$output) && rmarkdown::metadata$output
    ## %in% : 'length(x) = 3 > 1' in coercion to 'logical(1)'

<div
style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">

x

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

Ship-Ballast-104-02-05-17-SA-2-Dis2-1

</td>
<td style="text-align:right;">

0

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-105-02-05-17-SA-2-Dis2-2

</td>
<td style="text-align:right;">

0

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-107-02-05-17-SA-2-Dis2B-1

</td>
<td style="text-align:right;">

0

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-108-02-05-17-SA-2-Dis2B-2

</td>
<td style="text-align:right;">

0

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-109-02-05-17-SA-2-Dis2B-3

</td>
<td style="text-align:right;">

41

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-106-02-05-17-SA-2-Dis2-3

</td>
<td style="text-align:right;">

179

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-99-02-01-17-SA-2-tank-6-uptake-2

</td>
<td style="text-align:right;">

1640

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-103-02-01-17-SA-2-tank-2-uptake-3

</td>
<td style="text-align:right;">

1808

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-102-02-01-17-SA-2-tank-2-uptake-2

</td>
<td style="text-align:right;">

1957

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-98-02-01-17-SA-2-tank-6-uptake-1

</td>
<td style="text-align:right;">

2114

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-101-02-01-17-SA-2-tank-2-uptake-1

</td>
<td style="text-align:right;">

2169

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-97-10-23-16-SA-1-Dis-5-3

</td>
<td style="text-align:right;">

2703

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-82-10-18-16-SA-1-UP-5-3

</td>
<td style="text-align:right;">

3542

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-95-10-23-16-SA-1-Dis-5-1

</td>
<td style="text-align:right;">

3569

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-100-02-01-17-SA-2-tank-6-uptake-3

</td>
<td style="text-align:right;">

3676

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-81-10-18-16-SA-1-UP-5-2

</td>
<td style="text-align:right;">

3726

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-86-10-20-16-EX-Dis-5-1

</td>
<td style="text-align:right;">

3855

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-115-02-05-17-SA-2-Dis2B-3

</td>
<td style="text-align:right;">

4006

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-87-10-20-16-EX-Dis-5-2

</td>
<td style="text-align:right;">

4194

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-96-10-23-16-SA-1-Dis-5-2

</td>
<td style="text-align:right;">

4263

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-80-10-18-16-SA-1-UP-5-1

</td>
<td style="text-align:right;">

4416

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-83-10-17-16-SA-1-UP-6-1

</td>
<td style="text-align:right;">

4701

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-111-02-05-17-SA-2-Dis6-2

</td>
<td style="text-align:right;">

4937

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-88-10-20-16-EX-Dis-5-3

</td>
<td style="text-align:right;">

5198

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-85-10-17-16-SA-1-UP-6-3

</td>
<td style="text-align:right;">

5375

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-84-10-17-16-SA-1-UP-6-2

</td>
<td style="text-align:right;">

5588

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-113-02-05-17-SA-2-Dis2B-1

</td>
<td style="text-align:right;">

5964

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-112-02-05-17-SA-2-Dis6-3

</td>
<td style="text-align:right;">

6036

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-110-02-05-17-SA-2-Dis6-1

</td>
<td style="text-align:right;">

6570

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-114-02-05-17-SA-2-Dis2B-2

</td>
<td style="text-align:right;">

7009

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-94-10-23-16-SA-1-Dis-6-3

</td>
<td style="text-align:right;">

7257

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-93-10-23-16-SA-1-Dis-6-2

</td>
<td style="text-align:right;">

8050

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-91-10-20-16-EX-UP-5-3

</td>
<td style="text-align:right;">

8904

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-92-10-23-16-SA-1-Dis-6-1

</td>
<td style="text-align:right;">

9436

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-89-10-20-16-EX-UP-5-1

</td>
<td style="text-align:right;">

10467

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-90-10-20-16-EX-UP-5-2

</td>
<td style="text-align:right;">

10776

</td>
</tr>
</tbody>
</table>

</div>

## 2.3 Create a rarefaction curve

``` r
#Rarefaction curve using vegan
#From: https://micca.readthedocs.io/en/latest/phyloseq.html
taxa_are_rows(Ballast_physeq)
```

    ## [1] TRUE

``` r
mat <- t(otu_table(Ballast_physeq))
class(mat) <- "matrix"
class(mat)
```

    ## [1] "matrix" "array"

``` r
mat <- as(t(otu_table(Ballast_physeq)), "matrix")
class(mat)
```

    ## [1] "matrix" "array"

``` r
raremax <- min(rowSums(mat))

system.time(rarecurve(mat, step = 100, sample = raremax, col = "blue", label = FALSE))
```

    ## empty rows removed

![](16S-sequence-analysis_files/figure-gfm/rarefaction-1.png)<!-- -->

    ##    user  system elapsed 
    ##    0.11    0.00    0.13

## 2.4 Filter samples

Based on the table above, there are some samples with very few to no
reads, so we’ll filter these out here. In this case, we keep all samples
with more than 1,000 reads.

Removing these samples with low sequencing depths is important because
samples with low reads typically are lower quality and have a greater
probability of containing contaminant sequences ([Weiss et al.,
2017](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y)).

``` r
#Remove samples with few reads
Ballast_physeq = subset_samples(Ballast_physeq, 
                            Sample_number != "104")
Ballast_physeq = subset_samples(Ballast_physeq, 
                            Sample_number != "105")
Ballast_physeq = subset_samples(Ballast_physeq, 
                            Sample_number != "107")
Ballast_physeq = subset_samples(Ballast_physeq, 
                            Sample_number != "108")
Ballast_physeq = subset_samples(Ballast_physeq, 
                            Sample_number != "109")
Ballast_physeq = subset_samples(Ballast_physeq, 
                            Sample_number != "106")

#Filter out low abundance ASVs
Ballast_physeq = prune_taxa(taxa_sums(Ballast_physeq) > 2, Ballast_physeq)
```

## 2.5 Label and order variables

Now we can relabel the variables; if you type in Ballast_physeq you will
see that sample_data is the matrix that holds the information that we
want to change, so we need to include sample_data in our code here.

``` r
#Order factors
sample_data(Ballast_physeq)$Sample_type <- factor(sample_data(Ballast_physeq)$Sample_type, 
               levels = c("port_uptake", "ocean_uptake", "BWT", "BWT_BWE"),
               labels = c("Port uptake", "Ocean uptake", "BWT", "BWT+BWE"))

sample_data(Ballast_physeq)$Sample_number <- factor(sample_data(Ballast_physeq)$Sample_number, 
               levels = c("80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "100", "101", "102", "103", "110", "111", "112", "113", "114", "115"),
               labels = c("80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "100", "101", "102", "103", "110", "111", "112", "113", "114", "115"))

sample_data(Ballast_physeq)$Voyage <- factor(sample_data(Ballast_physeq)$Voyage, 
               levels = c("1", "2"),
               labels = c("1", "2"))
```

## 2.6 Abundance transformation

We have to transform the count data to account for differences in
library size between samples. This is often done using rarefaction, but
rarefaction results in the elimination of valid data and can also
increase the rate of false positives when testing for ASV’s that are
differentially abundant in different sample categories.

Additionally, here are some helpful references regarding library
normalization strategies:

[Waste Not, Want Not: Why Rarefying Microbiome Data is
Inadmissible](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)

[Normalization and microbial differential abundance strategies depend on
data
characteristics](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y)

Below, we’ll use relative abundance as our variance-stabilizing method.

``` r
#Transform to relative abundance
Ballast_physeqRA = transform_sample_counts(Ballast_physeq, function(x){x / sum(x)})
```

Alternatively, you could also choose to log-transform the data:

``` r
Ballast_physeqlog = transform_sample_counts(Ballast_physeq, function(x) log(1 + x))
```

# 3 Alpha Diversity

Note that while these data are not rarefied, they are normalized. We’ll
use the Shannon diversity index (based on richness AND evenness;
examines how many different taxa are present and how evenly they’re
distributed within a sample) to analyze alpha diversity between variable
types. This means it considers both the number of species and the
inequality between species abundances.

## 3.1 Calculate Shannon Diversity per Sample

Using the estimate_richness function in the phyloseq package to
calculate Shannon diversity. The estimate_richness function can also
take measures “Chao1” “ACE” “Simpson” and “Fisher”.

Note that running this code results in an error message stating that the
provided data does not have any singletons. See here for more
information:

<https://github.com/benjjneb/dada2/issues/214>

<https://forum.qiime2.org/t/singletons-and-diversity-richness-indices/2971/9>

``` r
#Calculating richness - Shannon diversity in a new dataframe
richness <- data.frame(estimate_richness(Ballast_physeqRA, measures = c("Shannon")))
richness <- setNames(cbind(rownames(richness), richness, row.names = NULL), 
                     c("sample-id", "Shannon"))

#Add the sample metadata to the dataframe
s <- data.frame(sample_data(Ballast_physeqRA))
s <- setNames(cbind(rownames(s), s, row.names = NULL), 
              c("sample-id", "Sample_type", "Voyage", "Tank", 
                "Sample_number", "Sample_site", "Sample_date"))

alphadiv <- merge(s, richness, by = "sample-id")

#Order factors
alphadiv$Sample_type <- factor(alphadiv$Sample_type, 
                      levels = c("Port uptake", "Ocean uptake", "BWT", "BWT+BWE"),
                      labels = c("Port uptake", "Ocean uptake", "BWT", "BWT+BWE")) 

#Shows the calculated indices
knitr::kable(head(alphadiv)) %>% 
  kableExtra::kable_styling("striped", latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div
style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">

sample-id

</th>
<th style="text-align:left;">

Sample_type

</th>
<th style="text-align:left;">

Voyage

</th>
<th style="text-align:right;">

Tank

</th>
<th style="text-align:left;">

Sample_number

</th>
<th style="text-align:left;">

Sample_site

</th>
<th style="text-align:left;">

Sample_date

</th>
<th style="text-align:left;">

NA

</th>
<th style="text-align:right;">

Shannon

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

Ship-Ballast-100-02-01-17-SA-2-tank-6-uptake-3

</td>
<td style="text-align:left;">

Port uptake

</td>
<td style="text-align:left;">

2

</td>
<td style="text-align:right;">

6

</td>
<td style="text-align:left;">

100

</td>
<td style="text-align:left;">

Antioch

</td>
<td style="text-align:left;">

1-Feb-17

</td>
<td style="text-align:left;">
</td>
<td style="text-align:right;">

4.053663

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-101-02-01-17-SA-2-tank-2-uptake-1

</td>
<td style="text-align:left;">

Port uptake

</td>
<td style="text-align:left;">

2

</td>
<td style="text-align:right;">

2

</td>
<td style="text-align:left;">

101

</td>
<td style="text-align:left;">

Antioch

</td>
<td style="text-align:left;">

1-Feb-17

</td>
<td style="text-align:left;">
</td>
<td style="text-align:right;">

3.984224

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-102-02-01-17-SA-2-tank-2-uptake-2

</td>
<td style="text-align:left;">

Port uptake

</td>
<td style="text-align:left;">

2

</td>
<td style="text-align:right;">

2

</td>
<td style="text-align:left;">

102

</td>
<td style="text-align:left;">

Antioch

</td>
<td style="text-align:left;">

1-Feb-17

</td>
<td style="text-align:left;">
</td>
<td style="text-align:right;">

3.999169

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-103-02-01-17-SA-2-tank-2-uptake-3

</td>
<td style="text-align:left;">

Port uptake

</td>
<td style="text-align:left;">

2

</td>
<td style="text-align:right;">

2

</td>
<td style="text-align:left;">

103

</td>
<td style="text-align:left;">

Antioch

</td>
<td style="text-align:left;">

1-Feb-17

</td>
<td style="text-align:left;">
</td>
<td style="text-align:right;">

3.755243

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-110-02-05-17-SA-2-Dis6-1

</td>
<td style="text-align:left;">

BWT

</td>
<td style="text-align:left;">

2

</td>
<td style="text-align:right;">

6

</td>
<td style="text-align:left;">

110

</td>
<td style="text-align:left;">

Port McNeill

</td>
<td style="text-align:left;">

7-Feb-17

</td>
<td style="text-align:left;">

treatment only

</td>
<td style="text-align:right;">

1.770741

</td>
</tr>
<tr>
<td style="text-align:left;">

Ship-Ballast-111-02-05-17-SA-2-Dis6-2

</td>
<td style="text-align:left;">

BWT

</td>
<td style="text-align:left;">

2

</td>
<td style="text-align:right;">

6

</td>
<td style="text-align:left;">

111

</td>
<td style="text-align:left;">

Port McNeill

</td>
<td style="text-align:left;">

7-Feb-17

</td>
<td style="text-align:left;">

treatment only

</td>
<td style="text-align:right;">

1.681808

</td>
</tr>
</tbody>
</table>

</div>

## 3.2 Create Alpha Diversity Boxplots

### 3.2.1 Comparison: sample type

``` r
#Plot
ggplot(data=alphadiv, aes(x=Sample_type, y=Shannon), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  geom_point(position=position_dodge(width=0.75),aes(group=Sample_type)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank())
```

![](16S-sequence-analysis_files/figure-gfm/boxplot-1.png)<!-- -->

### 3.2.2 Comparison: sample type and voyage

``` r
#Plot
ggplot(data=alphadiv, aes(x=Voyage, y=Shannon), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  geom_point(position=position_dodge(width=0.75),aes(group=Sample_type)) +
  facet_wrap(~Sample_type, scale="free") +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank())
```

![](16S-sequence-analysis_files/figure-gfm/boxplot2-1.png)<!-- -->

## 3.3 Alpha Diversity statistics

Perform Wilcoxon test.

# 4 Ordination plots

## 4.1 Find plot color hex codes

``` r
#For points
brewer.pal(4, "Dark2")
```

    ## [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A"

``` r
#For lines
brewer.pal(3, "Set1")
```

    ## [1] "#E41A1C" "#377EB8" "#4DAF4A"

## 4.2 Create NMDS with the Bray-Curtis dissimilarity matrix

Here, we use the relative abundance transformed data to create an NMDS
plot.

Also, pay attention to the printout of the stress here - less than 0.05
means that the two dimensional plot is an excellent representation of
the data (i.e., the NMDS ordination is a ‘good fit’ for the data), \<0.1
and 0.2 are also okay, but \> or = to 0.3 means that the NMDS ordination
does not fit the data well.

``` r
all.nmds.source.ord <- ordinate(
  physeq = Ballast_physeqRA, 
  method = "NMDS", 
  distance = "bray"
)
```

    ## Run 0 stress 0.1420328 
    ## Run 1 stress 0.07650127 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1432305  max resid 0.3388684 
    ## Run 2 stress 0.07650665 
    ## ... Procrustes: rmse 0.002126616  max resid 0.008201207 
    ## ... Similar to previous best
    ## Run 3 stress 0.07650152 
    ## ... Procrustes: rmse 9.187896e-05  max resid 0.0002036735 
    ## ... Similar to previous best
    ## Run 4 stress 0.07650038 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001382065  max resid 0.005338761 
    ## ... Similar to previous best
    ## Run 5 stress 0.07650611 
    ## ... Procrustes: rmse 0.002310297  max resid 0.0078514 
    ## ... Similar to previous best
    ## Run 6 stress 0.07650583 
    ## ... Procrustes: rmse 0.001930166  max resid 0.007550624 
    ## ... Similar to previous best
    ## Run 7 stress 0.07650116 
    ## ... Procrustes: rmse 0.0008765777  max resid 0.003586621 
    ## ... Similar to previous best
    ## Run 8 stress 0.07650568 
    ## ... Procrustes: rmse 0.001941535  max resid 0.007508334 
    ## ... Similar to previous best
    ## Run 9 stress 0.07650702 
    ## ... Procrustes: rmse 0.002520907  max resid 0.007910959 
    ## ... Similar to previous best
    ## Run 10 stress 0.07650597 
    ## ... Procrustes: rmse 0.001922496  max resid 0.007588458 
    ## ... Similar to previous best
    ## Run 11 stress 0.07650634 
    ## ... Procrustes: rmse 0.002377276  max resid 0.007741155 
    ## ... Similar to previous best
    ## Run 12 stress 0.07650127 
    ## ... Procrustes: rmse 0.001379364  max resid 0.005331976 
    ## ... Similar to previous best
    ## Run 13 stress 0.07650126 
    ## ... Procrustes: rmse 0.0008768257  max resid 0.003737948 
    ## ... Similar to previous best
    ## Run 14 stress 0.07650569 
    ## ... Procrustes: rmse 0.00219345  max resid 0.007711747 
    ## ... Similar to previous best
    ## Run 15 stress 0.1299754 
    ## Run 16 stress 0.07650107 
    ## ... Procrustes: rmse 0.0008753033  max resid 0.003622702 
    ## ... Similar to previous best
    ## Run 17 stress 0.1299782 
    ## Run 18 stress 0.1299754 
    ## Run 19 stress 0.1339471 
    ## Run 20 stress 0.1483835 
    ## *** Best solution repeated 12 times

``` r
#Plot, color coding by Sample_type
all.nmds.sample.type <- plot_ordination(
  physeq = Ballast_physeqRA,
  ordination = all.nmds.source.ord) + 
  scale_colour_manual(values = c("Port uptake" = "#E7298A", "Ocean uptake" = "#7570B3", "BWT" = "#1B9E77", "BWT+BWE" = "#D95F02"), "Sample Type") +
  scale_shape_manual(values = c("1" = 16, "2" = 17), name = "Voyage") +
  geom_point(mapping = aes(colour = factor(Sample_type), shape = factor(Voyage), size = 5)) +
  guides(size=FALSE) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank(), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.justification = c("right", "top")) 
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(all.nmds.sample.type)
```

![](16S-sequence-analysis_files/figure-gfm/NMDS%20Bray-1.png)<!-- -->

## 4.3 Create PCoA with the Bray-Curtis dissimilarity matrix

``` r
all.pcoa.source.ord <- ordinate(
  physeq = Ballast_physeqRA, 
  method = "PCoA", 
  distance = "bray"
)

#Plot, color coding by Sample_type
all.pcoa.sample.type <- plot_ordination(
  physeq = Ballast_physeqRA,
  ordination = all.pcoa.source.ord) + 
  scale_colour_manual(values = c("Port uptake" = "#E7298A", "Ocean uptake" = "#7570B3", "BWT" = "#1B9E77", "BWT+BWE" = "#D95F02"), "Sample Type") +
  scale_shape_manual(values = c("1" = 16, "2" = 17), name = "Voyage") +
  geom_point(mapping = aes(colour = factor(Sample_type), shape = factor(Voyage), size = 5)) +
  geom_text(aes(label = Sample_number, colour = factor(Sample_type)), nudge_x=0.05, nudge_y=0.05, check_overlap = TRUE) +
  guides(size=FALSE) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank(), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.justification = c("right", "top")) +
  geom_segment(
  aes(x = -0.45, y = -0.30, xend = 0.30, yend = -0.30),
  arrow = arrow(length = unit(0.03, "npc"))) +
   geom_segment(
  aes(x = -0.08, y = 0.03, xend = 0.29, yend = -0.08),
  arrow = arrow(length = unit(0.03, "npc")), colour = "#E41A1C") +
   geom_segment(
  aes(x = -0.14, y = 0.35, xend = 0.04, yend = 0.34),
  arrow = arrow(length = unit(0.03, "npc"))) +
  geom_curve(
  aes(x = -0.14, y = 0.43, xend = 0.25, yend = 0.21),
  arrow = arrow(length = unit(0.03, "npc")), curvature = -0.3) +
  geom_segment(
  aes(x = 0.06, y = 0.34, xend = 0.29, yend = -0.08),
  arrow = arrow(length = unit(0.03, "npc")), colour = "#E41A1C") 
  

print(all.pcoa.sample.type)
```

    ## Warning: Removed 1 rows containing missing values (`geom_text()`).

![](16S-sequence-analysis_files/figure-gfm/PCoA%20Bray-1.png)<!-- -->

### 4.3.1

# 5 Preliminary taxonomy fix

Some of the taxonomic levels in the phyloseq object have NA’s or short
values such as “g\_” which are uninformative when graphing bacterial
taxonomy at those levels. The tax_fix function from the microViz package
fills in these unknown values with useful information from higher
taxonomic levels.

This should be run prior to running the following code for generating
community composition barcharts, heatmaps, or any other analysis
involving the visualization of bacterial community composition.

``` r
#First, view the tax_table
knitr::kable(head(tax_table(Ballast_physeq))) %>% 
  kableExtra::kable_styling("striped") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div
style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">

Kingdom

</th>
<th style="text-align:left;">

Phylum

</th>
<th style="text-align:left;">

Class

</th>
<th style="text-align:left;">

Order

</th>
<th style="text-align:left;">

Family

</th>
<th style="text-align:left;">

Genus

</th>
<th style="text-align:left;">

Species

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

632d2cbb9c3417e0801b120443ee6a6d

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Gammaproteobacteria

</td>
<td style="text-align:left;">

Burkholderiales

</td>
<td style="text-align:left;">

T34

</td>
<td style="text-align:left;">

T34

</td>
<td style="text-align:left;">

uncultured_beta

</td>
</tr>
<tr>
<td style="text-align:left;">

8b66f81cf7f4e30dca76fb4e533fa0e2

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Gammaproteobacteria

</td>
<td style="text-align:left;">

Burkholderiales

</td>
<td style="text-align:left;">

T34

</td>
<td style="text-align:left;">

T34

</td>
<td style="text-align:left;">

uncultured_beta

</td>
</tr>
<tr>
<td style="text-align:left;">

2fb07775d5e453e95f73c4c347db7785

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Gammaproteobacteria

</td>
<td style="text-align:left;">

Burkholderiales

</td>
<td style="text-align:left;">

T34

</td>
<td style="text-align:left;">

T34

</td>
<td style="text-align:left;">

uncultured_beta

</td>
</tr>
<tr>
<td style="text-align:left;">

1fd035b1d5f122ae0ff5bed5979d964d

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Alphaproteobacteria

</td>
<td style="text-align:left;">

SAR11_clade

</td>
<td style="text-align:left;">

Clade_I

</td>
<td style="text-align:left;">

Clade_Ia

</td>
<td style="text-align:left;">

NA

</td>
</tr>
<tr>
<td style="text-align:left;">

1a67cc1cdac43c656cc5051897569de1

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Alphaproteobacteria

</td>
<td style="text-align:left;">

SAR11_clade

</td>
<td style="text-align:left;">

Clade_I

</td>
<td style="text-align:left;">

Clade_Ia

</td>
<td style="text-align:left;">

NA

</td>
</tr>
<tr>
<td style="text-align:left;">

c1c5b88084f21a29027a786c5bbd9167

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Alphaproteobacteria

</td>
<td style="text-align:left;">

SAR11_clade

</td>
<td style="text-align:left;">

Clade_I

</td>
<td style="text-align:left;">

Clade_Ia

</td>
<td style="text-align:left;">

NA

</td>
</tr>
</tbody>
</table>

</div>

``` r
#Now fix the labels
Ballast_physeq <- tax_fix(Ballast_physeq)

#View the relabeled tax_table
knitr::kable(head(tax_table(Ballast_physeq))) %>% 
  kableExtra::kable_styling("striped") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div
style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">

Kingdom

</th>
<th style="text-align:left;">

Phylum

</th>
<th style="text-align:left;">

Class

</th>
<th style="text-align:left;">

Order

</th>
<th style="text-align:left;">

Family

</th>
<th style="text-align:left;">

Genus

</th>
<th style="text-align:left;">

Species

</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

632d2cbb9c3417e0801b120443ee6a6d

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Gammaproteobacteria

</td>
<td style="text-align:left;">

Burkholderiales

</td>
<td style="text-align:left;">

Burkholderiales Order

</td>
<td style="text-align:left;">

Burkholderiales Order

</td>
<td style="text-align:left;">

uncultured_beta

</td>
</tr>
<tr>
<td style="text-align:left;">

8b66f81cf7f4e30dca76fb4e533fa0e2

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Gammaproteobacteria

</td>
<td style="text-align:left;">

Burkholderiales

</td>
<td style="text-align:left;">

Burkholderiales Order

</td>
<td style="text-align:left;">

Burkholderiales Order

</td>
<td style="text-align:left;">

uncultured_beta

</td>
</tr>
<tr>
<td style="text-align:left;">

2fb07775d5e453e95f73c4c347db7785

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Gammaproteobacteria

</td>
<td style="text-align:left;">

Burkholderiales

</td>
<td style="text-align:left;">

Burkholderiales Order

</td>
<td style="text-align:left;">

Burkholderiales Order

</td>
<td style="text-align:left;">

uncultured_beta

</td>
</tr>
<tr>
<td style="text-align:left;">

1fd035b1d5f122ae0ff5bed5979d964d

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Alphaproteobacteria

</td>
<td style="text-align:left;">

SAR11_clade

</td>
<td style="text-align:left;">

Clade_I

</td>
<td style="text-align:left;">

Clade_Ia

</td>
<td style="text-align:left;">

Clade_Ia Genus

</td>
</tr>
<tr>
<td style="text-align:left;">

1a67cc1cdac43c656cc5051897569de1

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Alphaproteobacteria

</td>
<td style="text-align:left;">

SAR11_clade

</td>
<td style="text-align:left;">

Clade_I

</td>
<td style="text-align:left;">

Clade_Ia

</td>
<td style="text-align:left;">

Clade_Ia Genus

</td>
</tr>
<tr>
<td style="text-align:left;">

c1c5b88084f21a29027a786c5bbd9167

</td>
<td style="text-align:left;">

d\_\_Bacteria

</td>
<td style="text-align:left;">

Proteobacteria

</td>
<td style="text-align:left;">

Alphaproteobacteria

</td>
<td style="text-align:left;">

SAR11_clade

</td>
<td style="text-align:left;">

Clade_I

</td>
<td style="text-align:left;">

Clade_Ia

</td>
<td style="text-align:left;">

Clade_Ia Genus

</td>
</tr>
</tbody>
</table>

</div>

# 6 Bacterial community composition bar charts

## 6.1 Find color hex codes

I’ve included here a list of hex codes for colors that can be used to
identify taxa in plots of bacterial community composition; these hex
codes are visualized with the seecolor package.

``` r
bcc_hex <- print_color(c("black", "#440154FF", "#450659FF","#460B5EFF","#472D7AFF","#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF", "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF","#2E6E8EFF", "#2D718EFF", "#2B748EFF","#2A778EFF", "#297B8EFF","#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF","#1FA287FF", "#21A585FF", "#23A983FF","#25AC82FF", "#29AF7FFF", "#2DB27DFF","#32B67AFF", "#37B878FF", "#3CBC74FF","#57C766FF", "#5EC962FF","#A2DA37FF","#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF","darkorchid1", "darkorchid2", "darkorchid3","#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF","#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#FF6699","#F68D45FF", "#FCA537FF","#F6E726FF", "#F4ED27FF", "#00489C", "#CCCCCC", "#999999", "#A1C299","#300018"), type = "r")
```

    ## 
    ##  ------ c black #440154FF #450659FF #460B5EFF #472D7AFF #3B518BFF #3A548CFF #38598CFF #365C8DFF #34608DFF #33638DFF #31678EFF #306A8EFF #2E6E8EFF #2D718EFF #2B748EFF #2A778EFF #297B8EFF #1F958BFF #1F988BFF #1E9C89FF #1F9F88FF #1FA287FF #21A585FF #23A983FF #25AC82FF #29AF7FFF #2DB27DFF #32B67AFF #37B878FF #3CBC74FF #57C766FF #5EC962FF #A2DA37FF #DAE319FF #E4E419FF #ECE51BFF #F5E61FFF darkorchid1 darkorchid2 darkorchid3 #A21C9AFF #A62098FF #AB2394FF #AE2892FF #B22B8FFF #B6308BFF #BA3388FF #BE3885FF #C13B82FF #C53F7EFF #C8437BFF #CC4678FF #CE4B75FF #D14E72FF #D5536FFF #D7566CFF #DA5B69FF #DD5E66FF #E06363FF #FF6699 #F68D45FF #FCA537FF #F6E726FF #F4ED27FF #00489C #CCCCCC #999999 #A1C299 #300018 ------
    ## black               
    ## #440154FF           
    ## #450659FF           
    ## #460B5EFF           
    ## #472D7AFF           
    ## #3B518BFF           
    ## #3A548CFF           
    ## #38598CFF           
    ## #365C8DFF           
    ## #34608DFF           
    ## #33638DFF           
    ## #31678EFF           
    ## #306A8EFF           
    ## #2E6E8EFF           
    ## #2D718EFF           
    ## #2B748EFF           
    ## #2A778EFF           
    ## #297B8EFF           
    ## #1F958BFF           
    ## #1F988BFF           
    ## #1E9C89FF           
    ## #1F9F88FF           
    ## #1FA287FF           
    ## #21A585FF           
    ## #23A983FF           
    ## #25AC82FF           
    ## #29AF7FFF           
    ## #2DB27DFF           
    ## #32B67AFF           
    ## #37B878FF           
    ## #3CBC74FF           
    ## #57C766FF           
    ## #5EC962FF           
    ## #A2DA37FF           
    ## #DAE319FF           
    ## #E4E419FF           
    ## #ECE51BFF           
    ## #F5E61FFF           
    ## darkorchid1         
    ## darkorchid2         
    ## darkorchid3         
    ## #A21C9AFF           
    ## #A62098FF           
    ## #AB2394FF           
    ## #AE2892FF           
    ## #B22B8FFF           
    ## #B6308BFF           
    ## #BA3388FF           
    ## #BE3885FF           
    ## #C13B82FF           
    ## #C53F7EFF           
    ## #C8437BFF           
    ## #CC4678FF           
    ## #CE4B75FF           
    ## #D14E72FF           
    ## #D5536FFF           
    ## #D7566CFF           
    ## #DA5B69FF           
    ## #DD5E66FF           
    ## #E06363FF           
    ## #FF6699             
    ## #F68D45FF           
    ## #FCA537FF           
    ## #F6E726FF           
    ## #F4ED27FF           
    ## #00489C             
    ## #CCCCCC             
    ## #999999             
    ## #A1C299             
    ## #300018

``` r
print(bcc_hex)
```

    ## NULL

## 6.2 Family level, individual samples

``` r
#Examining BCC in each individual sample first
BallastSeqR_family <- Ballast_physeq %>%
  tax_glom(taxrank = "Family") %>%                    #Agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)}) %>% #Transform to rel. abundance
  psmelt() %>%                                        # Melt to long format
  group_by(Sample, 
           Kingdom, Phylum, Class, Order, Family) %>%
  filter(Abundance > 0.01) %>%                        #Filter 
  arrange(Class)

#Creating a new column that combines both the Class and Family level information
#so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed
BallastSeqR_family$Class_family <- paste(BallastSeqR_family$Class, BallastSeqR_family$Family, sep="_")

#Creating plot  
BallastSeqR_family_bar <- ggplot(BallastSeqR_family, aes(x = Sample, y = Abundance, fill = Class_family)) + 
  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +
  scale_fill_manual(values = c("black", 
                               
                               "#440154FF", "#450659FF","#460B5EFF",
                               
                               "#472D7AFF", 
                               
                               "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF", 
                               "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",
                               "#2E6E8EFF", "#2D718EFF", "#2B748EFF",
                               "#2A778EFF", "#297B8EFF",
                               
                               "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",
                               "#1FA287FF", "#21A585FF", "#23A983FF",
                               "#25AC82FF", "#29AF7FFF", "#2DB27DFF",
                               "#32B67AFF", "#37B878FF", "#3CBC74FF",
                               
                               "#57C766FF", "#5EC962FF",
                               
                               "#A2DA37FF",
                               
                               "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",
                               
                               "darkorchid1", "darkorchid2", "darkorchid3", 
                               
                               "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",
                               "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF",
                               
                               "#FF6699",
                               
                               "#F68D45FF", "#FCA537FF",
                               
                               "#F6E726FF", "#F4ED27FF", 
                               
                               
                               "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",
                               
                               "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",
                               "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",
                               "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018")) +
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 1%) \n") + xlab("Sample ID") +
  theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(), 
                        panel.grid.major = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8, colour="black"), 
                        axis.text.y = element_text(size=8, colour="black"), plot.title = element_text(hjust = 0.5), 
                        axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),
                        legend.text = element_text(size = 7))  
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(BallastSeqR_family_bar)
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20sample%20id-1.png)<!-- -->

``` r
#Make interactive plot
BallastSeqR_family_bar_interactive <- ggplotly(BallastSeqR_family_bar)

BallastSeqR_family_bar_interactive
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20sample%20id-2.png)<!-- -->

``` r
#Save html of interactive plot on knit
htmlwidgets::saveWidget(BallastSeqR_family_bar_interactive, "BallastSeqR_family_bar_interactive.html")

htmltools::tags$iframe(
  src=file.path(getwd(), "BallastSeqR_family_bar_interactive.html"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0"
)
```

<iframe src="C:/Users/SBROWN13/OneDrive - Environmental Protection Agency (EPA)/Ballast water/Ballast-water-bacteria/BallastSeqR_family_bar_interactive.html" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>

## 6.3 Family level, grouped by sample type

``` r
#Examining BCC by sample type
BallastSeqR_family_sample <- Ballast_physeq %>%
  tax_glom(taxrank = "Family") %>%                    #Agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)}) %>% #Transform to rel. abundance
  psmelt() %>%                                        # Melt to long format
  group_by(Sample_type, 
           Kingdom, Phylum, Class, Order, Family) %>%
  dplyr::summarize(Mean = 
                     mean(Abundance, na.rm=TRUE)) %>% #Calculate average
  filter(Mean > 0.01) %>%                             #Filter 
  arrange(Class)
```

    ## `summarise()` has grouped output by 'Sample_type', 'Kingdom', 'Phylum',
    ## 'Class', 'Order'. You can override using the `.groups` argument.

``` r
#Creating a new column that combines both the Class and Family level information
#so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed
BallastSeqR_family_sample$Class_family <- paste(BallastSeqR_family_sample$Class, BallastSeqR_family_sample$Family, sep="_")

#Creating plot  
BallastSeqR_sample_family_bar <- ggplot(BallastSeqR_family_sample, aes(x = Sample_type, y = Mean, fill = Class_family)) + 
  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +
  scale_fill_manual(values = c("black", "#440154FF", "#450659FF","#460B5EFF","#472D7AFF", "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF","#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF","#2E6E8EFF", "#2D718EFF", "#2B748EFF","#2A778EFF", "#297B8EFF","#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",
"#1FA287FF", "#21A585FF", "#23A983FF","#25AC82FF", "#29AF7FFF", "#2DB27DFF","#32B67AFF", "#37B878FF", "#3CBC74FF","#57C766FF", "#5EC962FF","#A2DA37FF",
"#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF","darkorchid1", "darkorchid2", "darkorchid3", "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF","#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#FF6699","#F68D45FF", "#FCA537FF","#F6E726FF", "#F4ED27FF", "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018")) +
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 1%) \n") + xlab("Sample ID") +
  theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(), 
                        panel.grid.major = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8, colour="black"), 
                        axis.text.y = element_text(size=8, colour="black"), plot.title = element_text(hjust = 0.5), 
                        axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),
                        legend.text = element_text(size = 7))  

print(BallastSeqR_sample_family_bar)
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20sample%20type-1.png)<!-- -->

``` r
#Make interactive plot
BallastSeqR_sample_family_bar_interactive <- ggplotly(BallastSeqR_sample_family_bar)

BallastSeqR_sample_family_bar_interactive
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20sample%20type-2.png)<!-- -->

``` r
#Save html of interactive plot on knit
htmlwidgets::saveWidget(BallastSeqR_sample_family_bar_interactive, "BallastSeqR_sample_family_bar_interactive.html")

htmltools::tags$iframe(
  src=file.path(getwd(), "BallastSeqR_sample_family_bar_interactive.html"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0"
)
```

<iframe src="C:/Users/SBROWN13/OneDrive - Environmental Protection Agency (EPA)/Ballast water/Ballast-water-bacteria/BallastSeqR_sample_family_bar_interactive.html" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>

## 6.4 Family level, grouped by sample type and voyage

Creating the same bar plot as above, but with subsets within sample type
showing the different voyages, since similar sample types that were
collected during different voyages cluster seperately on the NMDS and
PCoA plots.

``` r
#Examining BCC by sample type & voyage
BallastSeqR_family_voyage <- Ballast_physeq %>%
  tax_glom(taxrank = "Family") %>%                    #Agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)}) %>% #Transform to rel. abundance
  psmelt() %>%                                        # Melt to long format
  group_by(Sample_type, Voyage, 
           Kingdom, Phylum, Class, Order, Family) %>%
  dplyr::summarize(Mean = 
                     mean(Abundance, na.rm=TRUE)) %>% #Calculate average
  filter(Mean > 0.01) %>%                             #Filter 
  arrange(Class)
```

    ## `summarise()` has grouped output by 'Sample_type', 'Voyage', 'Kingdom',
    ## 'Phylum', 'Class', 'Order'. You can override using the `.groups` argument.

``` r
#Creating a new column that combines both the Class and Family level information
#so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed
BallastSeqR_family_voyage$Class_family <- paste(BallastSeqR_family_voyage$Class, BallastSeqR_family_voyage$Family, sep="_")

#Creating plot  
BallastSeqR_voyage_family_bar <- ggplot(BallastSeqR_family_voyage, aes(x = Voyage, y = Mean, fill = Class_family)) + 
  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +
  scale_fill_manual(values = c("black", "#440154FF", "#450659FF","#460B5EFF","#472D7AFF", "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF","#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF","#2E6E8EFF", "#2D718EFF", "#2B748EFF","#2A778EFF", "#297B8EFF","#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",
"#1FA287FF", "#21A585FF", "#23A983FF","#25AC82FF", "#29AF7FFF", "#2DB27DFF","#32B67AFF", "#37B878FF", "#3CBC74FF","#57C766FF", "#5EC962FF","#A2DA37FF",
"#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF","darkorchid1", "darkorchid2", "darkorchid3", "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF","#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#FF6699","#F68D45FF", "#FCA537FF","#F6E726FF", "#F4ED27FF", "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018")) +
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 1%) \n") + xlab("Voyage") +
  theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"), 
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(), 
  axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8, 
                              colour="black"), 
  axis.text.y = element_text(size=8, colour="black"), 
  plot.title = element_text(hjust = 0.5), 
  axis.ticks.x = element_line(colour="#000000", size=0.1), 
  axis.ticks.y = element_line(colour="#000000", size=0.1),
  legend.text = element_text(size = 7)) +
  facet_grid(. ~ Sample_type, margins = TRUE, scale="free")

print(BallastSeqR_voyage_family_bar)
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20sample%20type%20and%20voyage-1.png)<!-- -->

``` r
#Make interactive plot
p <- ggplotly(BallastSeqR_voyage_family_bar)

p
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20sample%20type%20and%20voyage-2.png)<!-- -->

``` r
htmlwidgets::saveWidget(p, "BallastSeqR_voyage_family_bar_interactive.html")

htmltools::tags$iframe(
  src=file.path(getwd(), "BallastSeqR_voyage_family_bar_interactive.html"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0"
)
```

<iframe src="C:/Users/SBROWN13/OneDrive - Environmental Protection Agency (EPA)/Ballast water/Ballast-water-bacteria/BallastSeqR_voyage_family_bar_interactive.html" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>

# 7 DESeq2 analysis on MetaCyc data

## 7.1 Create the DESeq object

``` r
#Load metadata as a tibble
picrust_metadata <- read_delim("PICRUSt2-analysis/PICRUSt-Metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
```

    ## Rows: 30 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): sample-id, Sample_type, Sample_site, Sample_date, notes
    ## dbl (3): Voyage, Tank, Sample_number
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Load metacyc pathway abundance as a matrix
ballast_metacyc_abundance <- read_delim("PICRUSt2-analysis/MetaCyc_version/pathways_out/path_abun_unstrat_edited.tsv", delim = "\t")
```

    ## Rows: 378 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): pathway
    ## dbl (30): Ship-Ballast-100-02-01-17-SA-2-tank-6-uptake-3, Ship-Ballast-101-0...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Check to see if metadata and pathway abundance samples are in the same order
head(ballast_metacyc_abundance)
```

    ## # A tibble: 6 × 31
    ##   pathway   Ship-Ballast-100-02-…¹ Ship-Ballast-101-02-…² Ship-Ballast-102-02-…³
    ##   <chr>                      <dbl>                  <dbl>                  <dbl>
    ## 1 1CMET2-P…                  3278.                1996.                   1749. 
    ## 2 3-HYDROX…                   193.                  67.8                    88.8
    ## 3 AEROBACT…                     0                    2.63                    0  
    ## 4 ALL-CHOR…                   244.                 392.                    226. 
    ## 5 ANAEROFR…                  2706.                1473.                   1367. 
    ## 6 ANAGLYCO…                  3422.                2091.                   1892. 
    ## # ℹ abbreviated names: ¹​`Ship-Ballast-100-02-01-17-SA-2-tank-6-uptake-3`,
    ## #   ²​`Ship-Ballast-101-02-01-17-SA-2-tank-2-uptake-1`,
    ## #   ³​`Ship-Ballast-102-02-01-17-SA-2-tank-2-uptake-2`
    ## # ℹ 27 more variables: `Ship-Ballast-103-02-01-17-SA-2-tank-2-uptake-3` <dbl>,
    ## #   `Ship-Ballast-110-02-05-17-SA-2-Dis6-1` <dbl>,
    ## #   `Ship-Ballast-111-02-05-17-SA-2-Dis6-2` <dbl>,
    ## #   `Ship-Ballast-112-02-05-17-SA-2-Dis6-3` <dbl>, …

``` r
head(picrust_metadata)
```

    ## # A tibble: 6 × 8
    ##   `sample-id`     Sample_type Voyage  Tank Sample_number Sample_site Sample_date
    ##   <chr>           <chr>        <dbl> <dbl>         <dbl> <chr>       <chr>      
    ## 1 Ship-Ballast-8… port_uptake      1     5            80 Antioch     17-Oct-16  
    ## 2 Ship-Ballast-8… port_uptake      1     5            81 Antioch     17-Oct-16  
    ## 3 Ship-Ballast-8… port_uptake      1     5            82 Antioch     17-Oct-16  
    ## 4 Ship-Ballast-8… port_uptake      1     6            83 Antioch     17-Oct-16  
    ## 5 Ship-Ballast-8… port_uptake      1     6            84 Antioch     17-Oct-16  
    ## 6 Ship-Ballast-8… port_uptake      1     6            85 Antioch     17-Oct-16  
    ## # ℹ 1 more variable: notes <chr>

``` r
#In the picrust_metadata, make the sample-id column the rownames
picrust_metadata <- column_to_rownames(picrust_metadata, var = "sample-id")

#In the pathway abundance matrix, make the pathway column the rownames
ballast_metacyc_abundance <- column_to_rownames(ballast_metacyc_abundance, var = "pathway")

#The order of the rownames in the metadata file and the column names in the pathway abundance file don't match; they need to be in the same order before we can make this a DESeq object, because DESeq will assume that the sample order is consistent between the two
all(rownames(picrust_metadata) == colnames(ballast_metacyc_abundance))
```

    ## [1] FALSE

``` r
ballast_metacyc_abundance <- ballast_metacyc_abundance[, rownames(picrust_metadata)]
all(rownames(picrust_metadata) == colnames(ballast_metacyc_abundance))
```

    ## [1] TRUE

``` r
#Make Voyage & Sample_type factors
picrust_metadata$Voyage <- factor(picrust_metadata$Voyage, 
               levels = c("1", "2"),
               labels = c("1", "2"))

picrust_metadata$Sample_type <- factor(picrust_metadata$Sample_type, 
               levels = c("port_uptake", "ocean_uptake", "BWT", "BWT_BWE"),
               labels = c("Port_uptake", "Ocean_uptake", "BWT", "BWT_BWE"))

#Create a DESeq object from the MetaCyc count data; using both Voyage and Sample_type in the design formula in order to get a clearer picture of the differences in function that can be attributed to the different sample types
ballast_deseq <- DESeqDataSetFromMatrix(countData = round(ballast_metacyc_abundance),
                              colData = picrust_metadata, 
                              design = ~ Voyage + Sample_type)
```

    ## converting counts to integer mode

``` r
ballast_deseq
```

    ## class: DESeqDataSet 
    ## dim: 378 30 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(378): 1CMET2-PWY 3-HYDROXYPHENYLACETATE-DEGRADATION-PWY ...
    ##   VALDEG-PWY VALSYN-PWY
    ## rowData names(0):
    ## colnames(30): Ship-Ballast-80-10-18-16-SA-1-UP-5-1
    ##   Ship-Ballast-81-10-18-16-SA-1-UP-5-2 ...
    ##   Ship-Ballast-114-02-05-17-SA-2-Dis2B-2
    ##   Ship-Ballast-115-02-05-17-SA-2-Dis2B-3
    ## colData names(7): Sample_type Voyage ... Sample_date notes

``` r
#Make sure that the control is the reference level so that the default log2fold changes are calculated as treatment vs. control instead of vice versa
ballast_deseq$Sample_type <- relevel(ballast_deseq$Sample_type, ref = "Port_uptake")

#May add this part in later if I find that I have technical replicates: Collapse technical replicates (i.e., multiple sequencing runs of the same library) using the collapseReplicates function?
```

## 7.2 Load pathway descriptions

Load dataframe containing pathway descriptions; this will be used to add
descriptions to the pathways identified by DESeq2 later.

``` r
#Load pathways descriptions as a matrix
ballast_metacyc_descriptions <- read_delim("PICRUSt2-analysis/MetaCyc_version/pathways_out/path_abun_unstrat_descrip.tsv", delim = "\t")
```

    ## Rows: 378 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (2): pathway, description
    ## dbl (30): Ship-Ballast-100-02-01-17-SA-2-tank-6-uptake-3, Ship-Ballast-101-0...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Keep just the pathway and descritpions columns in this dataframe
ballast_metacyc_descriptions <- select(ballast_metacyc_descriptions, c(pathway, description))
```

## 7.3 Create PCA of all MetaCyc results in DESeq2

``` r
#Transforming the DESeq data using a regularized log transformation
rld <- rlog(ballast_deseq, blind=TRUE)
```

    ## rlog() may take a few minutes with 30 or more samples,
    ## vst() is a much faster transformation

    ## -- note: fitType='parametric', but the dispersion trend was not well captured by the
    ##    function: y = a/x + b, and a local regression fit was automatically substituted.
    ##    specify fitType='local' or 'mean' to avoid this message next time.

``` r
#Plot the PCA (Principle Component Analysis, not to be confused with PCoA, Principle Coordinate Analysis)
plotPCA(rld, intgroup=c("Voyage", "Sample_type"))
```

![](16S-sequence-analysis_files/figure-gfm/MetaCyc%20PCA-1.png)<!-- -->

``` r
#Customize the PCA plot using the ggplot function
pcaData <- plotPCA(rld, intgroup=c("Voyage", "Sample_type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Sample_type, shape=Voyage)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  geom_point(mapping = aes(colour = factor(Sample_type), shape = factor(Voyage), size = 3)) +
  scale_colour_manual(values = c("Port_uptake" = "#E7298A", "Ocean_uptake" = "#7570B3", "BWT" = "#1B9E77", "BWT_BWE" = "#D95F02"), "Sample Type") +
  scale_shape_manual(values = c("1" = 16, "2" = 17), name = "Voyage") +
  guides(size=FALSE) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank(), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.justification = c("right", "top"))
```

![](16S-sequence-analysis_files/figure-gfm/MetaCyc%20PCA-2.png)<!-- -->

## 7.4 Create PCoA of all MetaCyc results

## 7.5 BWT vs port uptake

``` r
#Run DESeq
ballast_deseq <- DESeq(ballast_deseq, test="Wald")
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## -- note: fitType='parametric', but the dispersion trend was not well captured by the
    ##    function: y = a/x + b, and a local regression fit was automatically substituted.
    ##    specify fitType='local' or 'mean' to avoid this message next time.

    ## final dispersion estimates

    ## fitting model and testing

``` r
#Show results
res_ballast_deseq <- results(ballast_deseq)
res_ballast_deseq
```

    ## log2 fold change (MLE): Sample type BWT BWE vs Port uptake 
    ## Wald test p-value: Sample type BWT BWE vs Port uptake 
    ## DataFrame with 378 rows and 6 columns
    ##                                          baseMean log2FoldChange     lfcSE
    ##                                         <numeric>      <numeric> <numeric>
    ## 1CMET2-PWY                             4607.13014      0.0706633 0.0192123
    ## 3-HYDROXYPHENYLACETATE-DEGRADATION-PWY  151.84829     -9.9913875 1.8507659
    ## AEROBACTINSYN-PWY                         1.43575     -0.4873761 5.0440621
    ## ALL-CHORISMATE-PWY                      608.89594    -25.5438177 2.3107275
    ## ANAEROFRUCAT-PWY                       1649.02292     -1.5990363 0.5839131
    ## ...                                           ...            ...       ...
    ## TYRFUMCAT-PWY                          3450.10632      0.0656368 0.1326428
    ## UBISYN-PWY                             3855.84916      0.1652058 0.0620179
    ## UDPNAGSYN-PWY                          3612.18634      0.2046273 0.0569165
    ## VALDEG-PWY                                6.69169     -4.8453789 4.4467187
    ## VALSYN-PWY                             5864.68071     -0.0897082 0.0808705
    ##                                               stat      pvalue        padj
    ##                                          <numeric>   <numeric>   <numeric>
    ## 1CMET2-PWY                               3.6780148 2.35056e-04 6.84484e-04
    ## 3-HYDROXYPHENYLACETATE-DEGRADATION-PWY  -5.3985149 6.71948e-08 4.21705e-07
    ## AEROBACTINSYN-PWY                       -0.0966237 9.23025e-01 9.54105e-01
    ## ALL-CHORISMATE-PWY                     -11.0544483 2.08615e-28 1.51872e-26
    ## ANAEROFRUCAT-PWY                        -2.7384836 6.17232e-03 1.28384e-02
    ## ...                                            ...         ...         ...
    ## TYRFUMCAT-PWY                             0.494839 0.620713886   0.7383655
    ## UBISYN-PWY                                2.663841 0.007725418   0.0155362
    ## UDPNAGSYN-PWY                             3.595223 0.000324114   0.0009217
    ## VALDEG-PWY                               -1.089653 0.275866170   0.3703789
    ## VALSYN-PWY                               -1.109283 0.267308332   0.3671707

``` r
#Check which comparisons are contained in the deseq object:
resultsNames(ballast_deseq)
```

    ## [1] "Intercept"                              
    ## [2] "Voyage_2_vs_1"                          
    ## [3] "Sample_type_Ocean_uptake_vs_Port_uptake"
    ## [4] "Sample_type_BWT_vs_Port_uptake"         
    ## [5] "Sample_type_BWT_BWE_vs_Port_uptake"

``` r
#Just want to compare Port_uptake with BWT; creating a results table with only this information
res_ballast_deseq_BWT <- results(ballast_deseq, contrast=c("Sample_type","BWT","Port_uptake"), alpha=0.05)

#Show results summary
summary(res_ballast_deseq_BWT)
```

    ## 
    ## out of 378 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 92, 24%
    ## LFC < 0 (down)     : 160, 42%
    ## outliers [1]       : 14, 3.7%
    ## low counts [2]     : 15, 4%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
#Log fold change shrinkage on the effluent_TWW_vs_CON comparison
resLFC <- lfcShrink(ballast_deseq, coef="Sample_type_BWT_vs_Port_uptake", type="apeglm")
```

    ## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    ##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    ##     sequence count data: removing the noise and preserving large differences.
    ##     Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
resLFC
```

    ## log2 fold change (MAP): Sample type BWT vs Port uptake 
    ## Wald test p-value: Sample type BWT vs Port uptake 
    ## DataFrame with 378 rows and 5 columns
    ##                                          baseMean log2FoldChange     lfcSE
    ##                                         <numeric>      <numeric> <numeric>
    ## 1CMET2-PWY                             4607.13014      0.0990272 0.0112923
    ## 3-HYDROXYPHENYLACETATE-DEGRADATION-PWY  151.84829     -0.7140748 0.7503179
    ## AEROBACTINSYN-PWY                         1.43575     -0.1990689 1.0423191
    ## ALL-CHORISMATE-PWY                      608.89594     -7.9981750 2.3141556
    ## ANAEROFRUCAT-PWY                       1649.02292     -0.5337356 0.3379172
    ## ...                                           ...            ...       ...
    ## TYRFUMCAT-PWY                          3450.10632       0.593383 0.0805587
    ## UBISYN-PWY                             3855.84916       0.264767 0.0375263
    ## UDPNAGSYN-PWY                          3612.18634       0.145036 0.0344510
    ## VALDEG-PWY                                6.69169      -0.488473 1.2035081
    ## VALSYN-PWY                             5864.68071      -0.302048 0.0487033
    ##                                             pvalue        padj
    ##                                          <numeric>   <numeric>
    ## 1CMET2-PWY                             2.58843e-19 2.24331e-18
    ## 3-HYDROXYPHENYLACETATE-DEGRADATION-PWY 1.83162e-01 2.33933e-01
    ## AEROBACTINSYN-PWY                      4.68474e-01 5.19892e-01
    ## ALL-CHORISMATE-PWY                     1.56946e-14 1.02015e-13
    ## ANAEROFRUCAT-PWY                       9.10885e-02 1.22801e-01
    ## ...                                            ...         ...
    ## TYRFUMCAT-PWY                          9.29827e-14 5.73656e-13
    ## UBISYN-PWY                             1.73965e-12 8.79491e-12
    ## UDPNAGSYN-PWY                          2.78790e-05 5.79883e-05
    ## VALDEG-PWY                             4.05068e-02 5.75956e-02
    ## VALSYN-PWY                             5.48814e-10 1.88461e-09

``` r
#"It is more useful visualize the MA-plot for the shrunken log2 fold changes, 
#which remove the noise associated with log2 fold changes from low count genes 
#without requiring arbitrary filtering thresholds."
#https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts
plotMA(resLFC)
```

![](16S-sequence-analysis_files/figure-gfm/DESeq2%20MetaCyc-1.png)<!-- -->

``` r
#Convert results to dataframe
df_ballast_deseq_BWT <- as.data.frame(res_ballast_deseq_BWT)

#Filter based on adj p-value, removing non-significant values (p>0.05). This also removes p-values that are NA; NA p-values exist because when p-values do not pass the filter threshold (p-values less than 0.1) they are set to NA. More information in vignette("DESeq2") under "Note on p-values set to NA"
df_ballast_deseq_BWT_filt <- filter(df_ballast_deseq_BWT, padj <= 0.05)

#Order by adjusted p-value
df_ballast_deseq_BWT_filt = df_ballast_deseq_BWT_filt[order(df_ballast_deseq_BWT_filt$padj, na.last=NA), ]

#Add pathway descriptions to the DESeq results using the left_join function in dplyr; rownames in the DESeq output need to be a column labeled pathway in order to match that of the descriptions dataframe
df_ballast_deseq_BWT_filt <- rownames_to_column(df_ballast_deseq_BWT_filt, var = "pathway")

df_ballast_deseq_BWT_filt_descrip <- df_ballast_deseq_BWT_filt %>%
  left_join(ballast_metacyc_descriptions, by = "pathway")
```

### 7.5.1 Create dot plot of results - log2fold change

``` r
#Plot
log2foldplot <- ggplot(df_ballast_deseq_BWT_filt_descrip, aes(x=log2FoldChange, y=description)) +  geom_point(size=6) +
  geom_vline(xintercept = 0, linetype="dashed", color = "black", size=1) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) 

print(log2foldplot)
```

![](16S-sequence-analysis_files/figure-gfm/log2fold%20change-1.png)<!-- -->

# 8 Session Info

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.2.1 (2022-06-23 ucrt)
    ##  os       Windows 10 x64 (build 22000)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       America/New_York
    ##  date     2023-08-25
    ##  pandoc   2.18 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package              * version    date (UTC) lib source
    ##    ade4                   1.7-22     2023-02-06 [1] CRAN (R 4.2.3)
    ##    annotate               1.76.0     2022-11-03 [1] Bioconductor
    ##    AnnotationDbi          1.60.2     2023-03-10 [1] Bioconductor
    ##    ape                  * 5.7-1      2023-03-13 [1] CRAN (R 4.2.3)
    ##    apeglm                 1.20.0     2022-11-01 [1] Bioconductor
    ##    backports              1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
    ##    base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
    ##    bbmle                  1.0.25     2022-05-11 [1] CRAN (R 4.2.3)
    ##    bdsmatrix              1.3-6      2022-06-03 [1] CRAN (R 4.2.0)
    ##    Biobase              * 2.58.0     2022-11-01 [1] Bioconductor
    ##    BiocGenerics         * 0.44.0     2022-11-01 [1] Bioconductor
    ##    BiocParallel           1.32.6     2023-03-17 [1] Bioconductor
    ##    biomformat             1.26.0     2022-11-01 [1] Bioconductor
    ##    Biostrings             2.66.0     2022-11-01 [1] Bioconductor
    ##    bit                    4.0.5      2022-11-15 [1] CRAN (R 4.2.3)
    ##    bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.2.3)
    ##    bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
    ##    blob                   1.2.4      2023-03-17 [1] CRAN (R 4.2.3)
    ##    cachem                 1.0.8      2023-05-01 [1] CRAN (R 4.2.3)
    ##    callr                  3.7.3      2022-11-02 [1] CRAN (R 4.2.3)
    ##    checkmate              2.2.0      2023-04-27 [1] CRAN (R 4.2.3)
    ##    circlize               0.4.15     2022-05-10 [1] CRAN (R 4.2.3)
    ##    cli                    3.6.1      2023-03-23 [1] CRAN (R 4.2.3)
    ##    clue                   0.3-64     2023-01-31 [1] CRAN (R 4.2.3)
    ##    cluster                2.1.4      2022-08-22 [1] CRAN (R 4.2.3)
    ##    coda                   0.19-4     2020-09-30 [1] CRAN (R 4.2.3)
    ##    codetools              0.2-19     2023-02-01 [1] CRAN (R 4.2.2)
    ##    colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.2.3)
    ##    ComplexHeatmap       * 2.14.0     2022-11-01 [1] Bioconductor
    ##    crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.2.3)
    ##    crosstalk              1.2.0      2021-11-04 [1] CRAN (R 4.2.3)
    ##    data.table             1.14.8     2023-02-17 [1] CRAN (R 4.2.3)
    ##    DBI                    1.1.3      2022-06-18 [1] CRAN (R 4.2.3)
    ##    DelayedArray           0.24.0     2022-11-01 [1] Bioconductor
    ##    DESeq2               * 1.38.3     2023-01-19 [1] Bioconductor
    ##    devtools               2.4.5      2022-10-11 [1] CRAN (R 4.2.3)
    ##    digest                 0.6.31     2022-12-11 [1] CRAN (R 4.2.3)
    ##    doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.2.3)
    ##    dplyr                * 1.1.2      2023-04-20 [1] CRAN (R 4.2.3)
    ##    DT                     0.28       2023-05-18 [1] CRAN (R 4.2.3)
    ##    ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.2.3)
    ##    emdbook                1.3.13     2023-07-03 [1] CRAN (R 4.2.1)
    ##    evaluate               0.21       2023-05-05 [1] CRAN (R 4.2.3)
    ##    fansi                  1.0.4      2023-01-22 [1] CRAN (R 4.2.3)
    ##    farver                 2.1.1      2022-07-06 [1] CRAN (R 4.2.3)
    ##    fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.2.3)
    ##    forcats              * 1.0.0      2023-01-29 [1] CRAN (R 4.2.3)
    ##    foreach                1.5.2      2022-02-02 [1] CRAN (R 4.2.3)
    ##    foreign                0.8-84     2022-12-06 [1] CRAN (R 4.2.2)
    ##    Formula                1.2-5      2023-02-24 [1] CRAN (R 4.2.2)
    ##    fs                     1.6.2      2023-04-25 [1] CRAN (R 4.2.3)
    ##    geneplotter            1.76.0     2022-11-01 [1] Bioconductor
    ##    generics               0.1.3      2022-07-05 [1] CRAN (R 4.2.3)
    ##    GenomeInfoDb         * 1.34.9     2023-02-02 [1] Bioconductor
    ##    GenomeInfoDbData       1.2.9      2023-04-07 [1] Bioconductor
    ##    GenomicRanges        * 1.50.2     2022-12-27 [1] Bioconductor
    ##    GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.2.3)
    ##    ggplot2              * 3.4.2      2023-04-03 [1] CRAN (R 4.2.3)
    ##    GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.2.3)
    ##    glue                   1.6.2      2022-02-24 [1] CRAN (R 4.2.3)
    ##    gridExtra              2.3        2017-09-09 [1] CRAN (R 4.2.3)
    ##    gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.2.1)
    ##    highr                  0.10       2022-12-22 [1] CRAN (R 4.2.3)
    ##    Hmisc                  5.1-0      2023-05-08 [1] CRAN (R 4.2.3)
    ##    hms                    1.1.3      2023-03-21 [1] CRAN (R 4.2.3)
    ##    htmlTable              2.4.1      2022-07-07 [1] CRAN (R 4.2.3)
    ##    htmltools              0.5.5      2023-03-23 [1] CRAN (R 4.2.3)
    ##    htmlwidgets            1.6.2      2023-03-17 [1] CRAN (R 4.2.3)
    ##    httpuv                 1.6.11     2023-05-11 [1] CRAN (R 4.2.3)
    ##    httr                   1.4.7      2023-08-15 [1] CRAN (R 4.2.3)
    ##    igraph                 1.4.3      2023-05-22 [1] CRAN (R 4.2.3)
    ##    IRanges              * 2.32.0     2022-11-01 [1] Bioconductor
    ##    iterators              1.0.14     2022-02-05 [1] CRAN (R 4.2.3)
    ##    jsonlite               1.8.5      2023-06-05 [1] CRAN (R 4.2.3)
    ##    kableExtra             1.3.4      2021-02-20 [1] CRAN (R 4.2.3)
    ##    KEGGREST               1.38.0     2022-11-01 [1] Bioconductor
    ##    knitr                  1.43       2023-05-25 [1] CRAN (R 4.2.3)
    ##    labeling               0.4.2      2020-10-20 [1] CRAN (R 4.2.0)
    ##    later                  1.3.1      2023-05-02 [1] CRAN (R 4.2.3)
    ##    lattice              * 0.21-8     2023-04-05 [1] CRAN (R 4.2.3)
    ##    lazyeval               0.2.2      2019-03-15 [1] CRAN (R 4.2.3)
    ##    lifecycle              1.0.3      2022-10-07 [1] CRAN (R 4.2.3)
    ##    locfit                 1.5-9.8    2023-06-11 [1] CRAN (R 4.2.3)
    ##    lubridate            * 1.9.2      2023-02-10 [1] CRAN (R 4.2.3)
    ##    magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.2.3)
    ##    MASS                   7.3-60     2023-05-04 [1] CRAN (R 4.2.3)
    ##    Matrix                 1.5-4.1    2023-05-18 [1] CRAN (R 4.2.3)
    ##    MatrixGenerics       * 1.10.0     2022-11-01 [1] Bioconductor
    ##    matrixStats          * 1.0.0      2023-06-02 [1] CRAN (R 4.2.3)
    ##    memoise                2.0.1      2021-11-26 [1] CRAN (R 4.2.3)
    ##    mgcv                   1.8-42     2023-03-02 [1] CRAN (R 4.2.3)
    ##    microViz             * 0.10.8     2023-05-01 [1] Github (david-barnett/microViz@f37d835)
    ##    mime                   0.12       2021-09-28 [1] CRAN (R 4.2.0)
    ##    miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.2.3)
    ##    multtest               2.54.0     2022-11-01 [1] Bioconductor
    ##    munsell                0.5.0      2018-06-12 [1] CRAN (R 4.2.3)
    ##    mvtnorm                1.2-2      2023-06-08 [1] CRAN (R 4.2.3)
    ##    NADA                   1.6-1.1    2020-03-22 [1] CRAN (R 4.2.3)
    ##    nlme                 * 3.1-162    2023-01-31 [1] CRAN (R 4.2.3)
    ##    nnet                   7.3-19     2023-05-03 [1] CRAN (R 4.2.3)
    ##    numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.2.0)
    ##    patchwork            * 1.1.2.9000 2023-04-24 [1] Github (thomasp85/patchwork@c14c960)
    ##    permute              * 0.9-7      2022-01-27 [1] CRAN (R 4.2.3)
    ##    phyloseq             * 1.42.0     2022-11-01 [1] Bioconductor
    ##    picante              * 1.8.2      2020-06-10 [1] CRAN (R 4.2.3)
    ##    pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.2.3)
    ##    pkgbuild               1.4.2      2023-06-26 [1] CRAN (R 4.2.1)
    ##    pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.2.3)
    ##    pkgload                1.3.2.1    2023-07-08 [1] CRAN (R 4.2.3)
    ##    plotly               * 4.10.2     2023-06-03 [1] CRAN (R 4.2.3)
    ##    plyr                 * 1.8.8      2022-11-11 [1] CRAN (R 4.2.3)
    ##    png                    0.1-8      2022-11-29 [1] CRAN (R 4.2.2)
    ##    prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.2.3)
    ##    processx               3.8.1      2023-04-18 [1] CRAN (R 4.2.1)
    ##    profvis                0.3.8      2023-05-02 [1] CRAN (R 4.2.3)
    ##    promises               1.2.0.1    2021-02-11 [1] CRAN (R 4.2.3)
    ##    ps                     1.7.5      2023-04-18 [1] CRAN (R 4.2.3)
    ##    purrr                * 1.0.1      2023-01-10 [1] CRAN (R 4.2.3)
    ##    qiime2R              * 0.99.6     2023-04-07 [1] Github (jbisanz/qiime2R@2a3cee1)
    ##    R6                     2.5.1      2021-08-19 [1] CRAN (R 4.2.3)
    ##    RColorBrewer         * 1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
    ##    Rcpp                   1.0.10     2023-01-22 [1] CRAN (R 4.2.3)
    ##    RCurl                  1.98-1.12  2023-03-27 [1] CRAN (R 4.2.3)
    ##    readr                * 2.1.4      2023-02-10 [1] CRAN (R 4.2.3)
    ##    remotes                2.4.2.1    2023-07-18 [1] CRAN (R 4.2.1)
    ##    reshape2             * 1.4.4      2020-04-09 [1] CRAN (R 4.2.3)
    ##    rhdf5                  2.42.0     2022-11-01 [1] Bioconductor
    ##  D rhdf5filters           1.10.1     2023-03-24 [1] Bioconductor
    ##    Rhdf5lib               1.20.0     2022-11-01 [1] Bioconductor
    ##    rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.2.0)
    ##    rlang                  1.1.1      2023-04-28 [1] CRAN (R 4.2.3)
    ##    rmarkdown              2.24       2023-08-14 [1] CRAN (R 4.2.3)
    ##    rpart                  4.1.19     2022-10-21 [1] CRAN (R 4.2.3)
    ##    RSQLite                2.3.1      2023-04-03 [1] CRAN (R 4.2.3)
    ##    rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.2.3)
    ##    rvest                  1.0.3      2022-08-19 [1] CRAN (R 4.2.3)
    ##    S4Vectors            * 0.36.2     2023-02-26 [1] Bioconductor
    ##    scales               * 1.2.1      2022-08-20 [1] CRAN (R 4.2.3)
    ##    seecolor             * 0.2.0      2023-02-24 [1] CRAN (R 4.2.3)
    ##    sessioninfo            1.2.2      2021-12-06 [1] CRAN (R 4.2.3)
    ##    shape                  1.4.6      2021-05-19 [1] CRAN (R 4.2.0)
    ##    shiny                  1.7.4.1    2023-07-06 [1] CRAN (R 4.2.3)
    ##    speedyseq            * 0.5.3.9018 2023-05-01 [1] Github (mikemc/speedyseq@ceb941f)
    ##    stringi                1.7.12     2023-01-11 [1] CRAN (R 4.2.2)
    ##    stringr              * 1.5.0      2022-12-02 [1] CRAN (R 4.2.3)
    ##    SummarizedExperiment * 1.28.0     2022-11-01 [1] Bioconductor
    ##    survival               3.5-5      2023-03-12 [1] CRAN (R 4.2.3)
    ##    svglite                2.1.1      2023-01-10 [1] CRAN (R 4.2.3)
    ##    systemfonts            1.0.4      2022-02-11 [1] CRAN (R 4.2.3)
    ##    tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.2.3)
    ##    tidyr                * 1.3.0      2023-01-24 [1] CRAN (R 4.2.3)
    ##    tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.2.3)
    ##    tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.2.3)
    ##    timechange             0.2.0      2023-01-11 [1] CRAN (R 4.2.3)
    ##    truncnorm              1.0-9      2023-03-20 [1] CRAN (R 4.2.3)
    ##    tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.2.3)
    ##    urlchecker             1.0.1      2021-11-30 [1] CRAN (R 4.2.3)
    ##    usethis                2.2.2      2023-07-06 [1] CRAN (R 4.2.3)
    ##    utf8                   1.2.3      2023-01-31 [1] CRAN (R 4.2.3)
    ##    vctrs                  0.6.2      2023-04-19 [1] CRAN (R 4.2.3)
    ##    vegan                * 2.6-4      2022-10-11 [1] CRAN (R 4.2.3)
    ##    viridis              * 0.6.3      2023-05-03 [1] CRAN (R 4.2.3)
    ##    viridisLite          * 0.4.2      2023-05-02 [1] CRAN (R 4.2.3)
    ##    vroom                  1.6.3      2023-04-28 [1] CRAN (R 4.2.3)
    ##    webshot                0.5.5      2023-06-26 [1] CRAN (R 4.2.1)
    ##    withr                  2.5.0      2022-03-03 [1] CRAN (R 4.2.3)
    ##    xfun                   0.39       2023-04-20 [1] CRAN (R 4.2.3)
    ##    XML                    3.99-0.14  2023-03-19 [1] CRAN (R 4.2.3)
    ##    xml2                   1.3.4      2023-04-27 [1] CRAN (R 4.2.3)
    ##    xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.2.3)
    ##    XVector                0.38.0     2022-11-01 [1] Bioconductor
    ##    yaml                   2.3.7      2023-01-23 [1] CRAN (R 4.2.3)
    ##    zCompositions          1.4.0-1    2022-03-26 [1] CRAN (R 4.2.3)
    ##    zlibbioc               1.44.0     2022-11-01 [1] Bioconductor
    ## 
    ##  [1] C:/Program Files/R/R-4.2.1/library
    ## 
    ##  D ── DLL MD5 mismatch, broken installation.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
