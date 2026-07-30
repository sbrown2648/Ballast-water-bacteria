Ballast water bacteria 16S rRNA analysis
================
Sarah Brown
2026-07-30

- [1 Load Packages and Import Data](#1-load-packages-and-import-data)
- [2 Preparing the Data](#2-preparing-the-data)
  - [2.1 Examine the control samples](#21-examine-the-control-samples)
  - [2.2 Examine the total number of reads and
    ASV’s](#22-examine-the-total-number-of-reads-and-asvs)
  - [2.3 Create table showing the number of reads per
    sample](#23-create-table-showing-the-number-of-reads-per-sample)
  - [2.4 Create a rarefaction curve](#24-create-a-rarefaction-curve)
  - [2.5 Filter samples](#25-filter-samples)
  - [2.6 Re-check negative-control ASVs after
    filtering](#26-re-check-negative-control-asvs-after-filtering)
  - [2.7 Label and order variables](#27-label-and-order-variables)
  - [2.8 Group rarefaction curves by
    Sample_type](#28-group-rarefaction-curves-by-sample_type)
  - [2.9 Abundance transformation](#29-abundance-transformation)
- [3 Alpha Diversity](#3-alpha-diversity)
  - [3.1 Calculate Shannon Diversity per
    Sample](#31-calculate-shannon-diversity-per-sample)
  - [3.2 Create Alpha Diversity
    Boxplots](#32-create-alpha-diversity-boxplots)
    - [3.2.1 Shannon diversity: sample
      type](#321-shannon-diversity-sample-type)
    - [3.2.2 Shannon diversity: sample type and
      voyage](#322-shannon-diversity-sample-type-and-voyage)
    - [3.2.3 Shannon diversity: voyage and sample
      type](#323-shannon-diversity-voyage-and-sample-type)
    - [3.2.4 Richness: sample type](#324-richness-sample-type)
    - [3.2.5 Richness: voyage and sample
      type](#325-richness-voyage-and-sample-type)
    - [3.2.6 Reads: voyage and sample
      type](#326-reads-voyage-and-sample-type)
    - [3.2.7 Combine plots for Fig. 2](#327-combine-plots-for-fig-2)
  - [3.3 Alpha Diversity statistics](#33-alpha-diversity-statistics)
    - [3.3.1 Subset alpha diversity
      dataframe](#331-subset-alpha-diversity-dataframe)
    - [3.3.2 Shannon diversity: sample
      type](#332-shannon-diversity-sample-type)
    - [3.3.3 Shannon diversity: voyage and sample
      type](#333-shannon-diversity-voyage-and-sample-type)
    - [3.3.4 Richness: sample type](#334-richness-sample-type)
    - [3.3.5 Richness: voyage and sample
      type](#335-richness-voyage-and-sample-type)
- [4 Ordination plots](#4-ordination-plots)
  - [4.1 Find plot color hex codes](#41-find-plot-color-hex-codes)
  - [4.2 Create NMDS with the Bray-Curtis dissimilarity
    matrix](#42-create-nmds-with-the-bray-curtis-dissimilarity-matrix)
  - [4.3 Create PCoA with the Bray-Curtis dissimilarity
    matrix](#43-create-pcoa-with-the-bray-curtis-dissimilarity-matrix)
    - [4.3.1 PCoA with labels and
      arrows](#431-pcoa-with-labels-and-arrows)
    - [4.3.2 PCoA without labels and
      arrows](#432-pcoa-without-labels-and-arrows)
- [5 PERMANOVA Analysis](#5-permanova-analysis)
  - [5.1 PERMANOVA analysis](#51-permanova-analysis)
  - [5.2 Pairwise PERMANOVA](#52-pairwise-permanova)
- [6 Preliminary taxonomy fix](#6-preliminary-taxonomy-fix)
- [7 Bacterial community composition bar
  charts](#7-bacterial-community-composition-bar-charts)
  - [7.1 Find color hex codes](#71-find-color-hex-codes)
  - [7.2 Family level, individual
    samples](#72-family-level-individual-samples)
  - [7.3 Family level, grouped by sample
    type](#73-family-level-grouped-by-sample-type)
  - [7.4 Family level, grouped by sample type and
    voyage](#74-family-level-grouped-by-sample-type-and-voyage)
  - [7.5 Genus level, grouped by voyage and sample
    type](#75-genus-level-grouped-by-voyage-and-sample-type)
  - [7.6 Species level, grouped by voyage and sample
    type](#76-species-level-grouped-by-voyage-and-sample-type)
    - [7.6.1 Create Figure 3: Faceted grouped bar
      plot](#761-create-figure-3-faceted-grouped-bar-plot)
  - [7.7 Species level, grouped by voyage, sample type, and
    tank](#77-species-level-grouped-by-voyage-sample-type-and-tank)
    - [7.7.1 Create plot with
      fantaxtic](#771-create-plot-with-fantaxtic)
- [8 Bacterial community composition dot
  plots](#8-bacterial-community-composition-dot-plots)
  - [8.1 Genus level, grouped by sample type and
    voyage](#81-genus-level-grouped-by-sample-type-and-voyage)
  - [8.2 Genus level, grouped by voyage and sample
    type](#82-genus-level-grouped-by-voyage-and-sample-type)
  - [8.3 Species level, grouped by sample type and
    voyage](#83-species-level-grouped-by-sample-type-and-voyage)
  - [8.4 Species level, grouped by voyage and sample
    type](#84-species-level-grouped-by-voyage-and-sample-type)
- [9 Filter ASVs to determine how many were given species-level
  assignments](#9-filter-asvs-to-determine-how-many-were-given-species-level-assignments)
- [10 Session Info](#10-session-info)

# 1 Load Packages and Import Data

Here, we load the required packages and then import the .qza files
exported from QIIME2 using the qiime2R package. The final object is a
phyloseq object called Ballast_physeq.

``` r
library(ggplot2)
library(vegan)
```

    ## Loading required package: permute

``` r
library(plyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, mutate, rename, summarise, summarize

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

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     explain

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

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

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

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

    ## Loading required package: Seqinfo

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
library(microViz)
```

    ## microViz version 0.13.1 - Copyright (C) 2021-2026 David Barnett
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
    ## ComplexHeatmap version 2.26.1
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
    ## ✔ forcats   1.0.1     ✔ readr     2.2.0
    ## ✔ lubridate 1.9.5     ✔ stringr   1.6.0
    ## ✔ purrr     1.2.2

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
    ## ✖ plotly::filter()         masks dplyr::filter(), stats::filter()
    ## ✖ S4Vectors::first()       masks dplyr::first()
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
library(ggh4x)
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

``` r
library(colorspace)

set.seed(13745)

# Reusable color palette for sample types
sample_type_colors <- c(
  "Port uptake" = "#E7298A",
  "Ocean uptake" = "#7570B3",
  "BWT" = "#1B9E77",
  "BWT+BWE" = "#D95F02"
)

# Reusable theme for clean plots
theme_ballast <- function() {
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.ticks.x = element_line(colour = "#000000", linewidth = 0.1),
    axis.ticks.y = element_line(colour = "#000000", linewidth = 0.1),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )
}

# Voyage shape mapping
voyage_shapes <- c("Voyage 1" = 16, "Voyage 2" = 17)

# Samples to exclude (low read counts)
samples_to_exclude <- c("104", "105", "106", "107", "108", "109")

# Factor levels and labels for sample types
sample_type_levels <- c("port_uptake", "BWT", "ocean_uptake", "BWT_BWE")
sample_type_labels <- c("Port uptake", "BWT", "Ocean uptake", "BWT+BWE")

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
    ## [5] "Sample_site"   "Sample_date"   "Type_voyage"   "Sample_group" 
    ## [9] "notes"

``` r
#Remove chloroplast sequences and any contaminant sequences
Ballast_physeq <- subset_taxa(Ballast_physeq, Kingdom != "d__Archaea")
Ballast_physeq <- subset_taxa(Ballast_physeq, Order != "Chloroplast")
Ballast_physeq <- subset_taxa(Ballast_physeq, Family != "Mitochondria")

#Check that contaminant sequences are removed (easiest to save as data frame and search)
taxtabl <- as.data.frame(tax_table(Ballast_physeq))
```

## 2.1 Examine the control samples

Before removing the blanks and controls, we examine how many reads they
contain, the identity (taxonomy) of the ASVs detected in the negative
controls (`Sample_type == "negative"`), and whether any of those ASVs
also occur in the true samples. ASVs that appear in both the negative
controls and the samples are potential contaminants.

``` r
#Examine the number of reads in the control samples before removing them
control_reads <- sort(sample_sums(subset_samples(Ballast_physeq, Sample_type %in% c("control", "negative"))))
knitr::kable(control_reads, col.names = "Number of reads") %>%
    kableExtra::kable_styling("striped", latex_options = "scale_down") %>%
    kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Number of reads
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Control-NPCR1-Negative-PCR
</td>

<td style="text-align:right;">

0
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR6-Negative-PCR-EMP
</td>

<td style="text-align:right;">

0
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR7-Negative-PCR
</td>

<td style="text-align:right;">

0
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR8-Negative-PCR
</td>

<td style="text-align:right;">

0
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR9-Negative-PCR-EMP
</td>

<td style="text-align:right;">

0
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR5-Negative-PCR-EMP
</td>

<td style="text-align:right;">

2
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR4-Negative-PCR-EMP
</td>

<td style="text-align:right;">

4
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR2-Negative-PCR
</td>

<td style="text-align:right;">

5
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-NPCR3-Negative-PCR
</td>

<td style="text-align:right;">

5
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM5-Mock-Microbiome-EMP
</td>

<td style="text-align:right;">

5344
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM4-Mock-Microbiome-EMP
</td>

<td style="text-align:right;">

6713
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM6-Mock-Microbiome-EMP
</td>

<td style="text-align:right;">

8250
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM2-Mock-Microbiome
</td>

<td style="text-align:right;">

8616
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM9-Mock-Microbiome-EMP
</td>

<td style="text-align:right;">

8809
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM1-Mock-Microbiome
</td>

<td style="text-align:right;">

9547
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM3-Mock-Microbiome
</td>

<td style="text-align:right;">

11487
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM7-Mock-Microbiome
</td>

<td style="text-align:right;">

12098
</td>

</tr>

<tr>

<td style="text-align:left;">

Control-MM8-Mock-Microbiome
</td>

<td style="text-align:right;">

12799
</td>

</tr>

</tbody>

</table>

</div>

``` r
#Identify the ASVs found in the negative controls
#Subset to the negative controls and keep only ASVs that actually occur in them
neg_controls <- subset_samples(Ballast_physeq, Sample_type == "negative")
neg_controls <- prune_taxa(taxa_sums(neg_controls) > 0, neg_controls)

#Reads for each ASV summed across the true (non-control) samples
true_samples <- subset_samples(Ballast_physeq, !(Sample_type %in% c("control", "negative")))
sample_taxa_sums <- taxa_sums(true_samples)

#Build a table of the negative-control ASV identities, their read counts in the
#negative controls, and whether/how much they also occur in the true samples
neg_control_taxa <- as.data.frame(tax_table(neg_controls))
neg_control_taxa$Reads_in_negatives <- as.integer(taxa_sums(neg_controls))
neg_control_taxa$Reads_in_samples <- as.integer(sample_taxa_sums[rownames(neg_control_taxa)])
neg_control_taxa$Present_in_samples <- neg_control_taxa$Reads_in_samples > 0

#Order by abundance in the negative controls
neg_control_taxa <- neg_control_taxa[order(-neg_control_taxa$Reads_in_negatives), ]

knitr::kable(neg_control_taxa, row.names = FALSE) %>%
    kableExtra::kable_styling("striped", latex_options = "scale_down") %>%
    kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

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

<th style="text-align:right;">

Reads_in_negatives
</th>

<th style="text-align:right;">

Reads_in_samples
</th>

<th style="text-align:left;">

Present_in_samples
</th>

</tr>

</thead>

<tbody>

<tr>

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

Enterobacterales
</td>

<td style="text-align:left;">

Enterobacteriaceae
</td>

<td style="text-align:left;">

Escherichia-Shigella
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:right;">

12
</td>

<td style="text-align:right;">

10
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:left;">

d\_\_Bacteria
</td>

<td style="text-align:left;">

Bacteroidota
</td>

<td style="text-align:left;">

Bacteroidia
</td>

<td style="text-align:left;">

Bacteroidales
</td>

<td style="text-align:left;">

Prevotellaceae
</td>

<td style="text-align:left;">

Prevotella
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

<tr>

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

Pseudomonadales
</td>

<td style="text-align:left;">

Pseudomonadaceae
</td>

<td style="text-align:left;">

Pseudomonas
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

</tbody>

</table>

</div>

``` r
#How many of the negative-control ASVs were also detected in the samples?
cat(sprintf("%d of %d ASVs detected in the negative controls were also present in the samples.\n",
            sum(neg_control_taxa$Present_in_samples), nrow(neg_control_taxa)))
```

    ## 2 of 3 ASVs detected in the negative controls were also present in the samples.

``` r
#At this point, blanks and positive controls should also be removed
Ballast_physeq = subset_samples(Ballast_physeq, Sample_type != "control")
Ballast_physeq = subset_samples(Ballast_physeq, Sample_type != "negative")
```

## 2.2 Examine the total number of reads and ASV’s

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
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank())

print(ASV_occurance)
```

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](16S-sequence-analysis_files/figure-gfm/ASV%20abundance-1.png)<!-- -->

## 2.3 Create table showing the number of reads per sample

``` r
#Sum reads and sort by samples with least to greatest number of reads
samples_reads <- sort(sample_sums(Ballast_physeq))

#Create table 
knitr::kable(samples_reads) %>% 
    kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
    kableExtra::scroll_box(width = "100%") 
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

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

## 2.4 Create a rarefaction curve

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
    ##    0.24    0.00    0.24

## 2.5 Filter samples

Based on the table above, there are some samples with very few to no
reads, so we’ll filter these out here. In this case, we keep all samples
with more than 1,000 reads.

Removing these samples with low sequencing depths is important because
samples with low reads typically are lower quality and have a greater
probability of containing contaminant sequences ([Weiss et al.,
2017](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y)).

``` r
#Remove samples with few reads
Ballast_physeq <- subset_samples(Ballast_physeq,
                                 !(Sample_number %in% samples_to_exclude))

#Filter out low abundance ASVs
Ballast_physeq = prune_taxa(taxa_sums(Ballast_physeq) > 2, Ballast_physeq)
```

## 2.6 Re-check negative-control ASVs after filtering

Now that the low-read samples and low-abundance ASVs have been removed,
we check whether the ASVs that were detected in the negative controls
are still present in the retained samples. ASVs that dropped out here
were only supported by the samples/reads that were filtered out.

``` r
#taxa_sums over the filtered object; ASVs that were pruned out are treated as 0 reads
filtered_taxa_sums <- taxa_sums(Ballast_physeq)

#Add post-filtering read counts to the negative-control ASV table built earlier
neg_control_after_filter <- neg_control_taxa
neg_control_after_filter$Reads_after_filtering <-
  as.integer(filtered_taxa_sums[rownames(neg_control_after_filter)])
neg_control_after_filter$Reads_after_filtering[is.na(neg_control_after_filter$Reads_after_filtering)] <- 0L
neg_control_after_filter$Present_after_filtering <- neg_control_after_filter$Reads_after_filtering > 0

knitr::kable(neg_control_after_filter, row.names = FALSE) %>%
    kableExtra::kable_styling("striped", latex_options = "scale_down") %>%
    kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

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

<th style="text-align:right;">

Reads_in_negatives
</th>

<th style="text-align:right;">

Reads_in_samples
</th>

<th style="text-align:left;">

Present_in_samples
</th>

<th style="text-align:right;">

Reads_after_filtering
</th>

<th style="text-align:left;">

Present_after_filtering
</th>

</tr>

</thead>

<tbody>

<tr>

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

Enterobacterales
</td>

<td style="text-align:left;">

Enterobacteriaceae
</td>

<td style="text-align:left;">

Escherichia-Shigella
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:right;">

12
</td>

<td style="text-align:right;">

10
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:right;">

10
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:left;">

d\_\_Bacteria
</td>

<td style="text-align:left;">

Bacteroidota
</td>

<td style="text-align:left;">

Bacteroidia
</td>

<td style="text-align:left;">

Bacteroidales
</td>

<td style="text-align:left;">

Prevotellaceae
</td>

<td style="text-align:left;">

Prevotella
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

<tr>

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

Pseudomonadales
</td>

<td style="text-align:left;">

Pseudomonadaceae
</td>

<td style="text-align:left;">

Pseudomonas
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

</tbody>

</table>

</div>

``` r
#Of the negative-control ASVs that overlapped with the samples, how many remain?
n_overlap <- sum(neg_control_after_filter$Present_in_samples)
n_overlap_retained <- sum(neg_control_after_filter$Present_in_samples &
                          neg_control_after_filter$Present_after_filtering)
cat(sprintf("Of the %d negative-control ASVs that were present in the samples, %d were still present after low-read samples and low-abundance ASVs were removed.\n",
            n_overlap, n_overlap_retained))
```

    ## Of the 2 negative-control ASVs that were present in the samples, 1 were still present after low-read samples and low-abundance ASVs were removed.

## 2.7 Label and order variables

Now we can relabel the variables; if you type in Ballast_physeq you will
see that sample_data is the matrix that holds the information that we
want to change, so we need to include sample_data in our code here.

``` r
#Order factors
sample_data(Ballast_physeq)$Sample_type <- factor(
  sample_data(Ballast_physeq)$Sample_type,
  levels = sample_type_levels,
  labels = sample_type_labels)

sample_numbers <- c("80", "81", "82", "83", "84", "85", "86", "87", "88",
                    "89", "90", "91", "92", "93", "94", "95", "96", "97",
                    "98", "100", "101", "102", "103", "110", "111", "112",
                    "113", "114", "115")
sample_data(Ballast_physeq)$Sample_number <- factor(
  sample_data(Ballast_physeq)$Sample_number,
  levels = sample_numbers)

sample_data(Ballast_physeq)$Voyage <- factor(
  sample_data(Ballast_physeq)$Voyage,
  levels = c("1", "2"),
  labels = c("Voyage 1", "Voyage 2"))
```

## 2.8 Group rarefaction curves by Sample_type

Comparing the sequencing depth and species richness between different
sample types. Using the ggrare function in the phyloseq.extended package
which wraps vegan’s rarecurve function so that it can be easily used
with phyloseq objects.

Note: I load the phyloseq-extended package here instead of above in the
‘Load packages’ section because it seems to create conflicts with the
phyloseq package and results in several warning messages.

``` r
library(phyloseq.extended)

p <- phyloseq.extended:::ggrare(Ballast_physeq, step = 1000, color = "Sample_type", label = "Sample_number", se = FALSE)
```

    ## rarefying sample Ship-Ballast-100-02-01-17-SA-2-tank-6-uptake-3

    ## rarefying sample Ship-Ballast-101-02-01-17-SA-2-tank-2-uptake-1

    ## rarefying sample Ship-Ballast-102-02-01-17-SA-2-tank-2-uptake-2

    ## rarefying sample Ship-Ballast-103-02-01-17-SA-2-tank-2-uptake-3

    ## rarefying sample Ship-Ballast-110-02-05-17-SA-2-Dis6-1

    ## rarefying sample Ship-Ballast-111-02-05-17-SA-2-Dis6-2

    ## rarefying sample Ship-Ballast-112-02-05-17-SA-2-Dis6-3

    ## rarefying sample Ship-Ballast-113-02-05-17-SA-2-Dis2B-1

    ## rarefying sample Ship-Ballast-114-02-05-17-SA-2-Dis2B-2

    ## rarefying sample Ship-Ballast-115-02-05-17-SA-2-Dis2B-3

    ## rarefying sample Ship-Ballast-80-10-18-16-SA-1-UP-5-1

    ## rarefying sample Ship-Ballast-81-10-18-16-SA-1-UP-5-2

    ## rarefying sample Ship-Ballast-82-10-18-16-SA-1-UP-5-3

    ## rarefying sample Ship-Ballast-83-10-17-16-SA-1-UP-6-1

    ## rarefying sample Ship-Ballast-84-10-17-16-SA-1-UP-6-2

    ## rarefying sample Ship-Ballast-85-10-17-16-SA-1-UP-6-3

    ## rarefying sample Ship-Ballast-86-10-20-16-EX-Dis-5-1

    ## rarefying sample Ship-Ballast-87-10-20-16-EX-Dis-5-2

    ## rarefying sample Ship-Ballast-88-10-20-16-EX-Dis-5-3

    ## rarefying sample Ship-Ballast-89-10-20-16-EX-UP-5-1

    ## rarefying sample Ship-Ballast-90-10-20-16-EX-UP-5-2

    ## rarefying sample Ship-Ballast-91-10-20-16-EX-UP-5-3

    ## rarefying sample Ship-Ballast-92-10-23-16-SA-1-Dis-6-1

    ## rarefying sample Ship-Ballast-93-10-23-16-SA-1-Dis-6-2

    ## rarefying sample Ship-Ballast-94-10-23-16-SA-1-Dis-6-3

    ## rarefying sample Ship-Ballast-95-10-23-16-SA-1-Dis-5-1

    ## rarefying sample Ship-Ballast-96-10-23-16-SA-1-Dis-5-2

    ## rarefying sample Ship-Ballast-97-10-23-16-SA-1-Dis-5-3

    ## rarefying sample Ship-Ballast-98-02-01-17-SA-2-tank-6-uptake-1

    ## rarefying sample Ship-Ballast-99-02-01-17-SA-2-tank-6-uptake-2

![](16S-sequence-analysis_files/figure-gfm/rarefaction%202-1.png)<!-- -->

``` r
p <- p + facet_wrap(~Sample_type, scale="free") +
  scale_colour_manual(values = sample_type_colors, "Sample Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = "Number of Reads")


print(p)
```

![](16S-sequence-analysis_files/figure-gfm/rarefaction%202-2.png)<!-- -->

``` r
#Save as svg
ggsave(filename="R output/rarecurve_sampletype.svg", plot=p, width=8, height=6, device=svg)
```

## 2.9 Abundance transformation

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
#Change the - in the sample names to _ to prevent the estimate_richness function from changing the - to a . during the calculation of Observed richness
sample_names(Ballast_physeq) <- paste(gsub("-", "_", sample_names(Ballast_physeq)), sep = "_")

#Calculating richness - Shannon diversity and number of observations (Observed) in a new dataframe
richness <- data.frame(estimate_richness(Ballast_physeq, measures = c("Shannon", "Observed")))

richness <- setNames(cbind(rownames(richness), richness, row.names = NULL), 
                     c("sample-id", "Observed", "Shannon"))

#Add the sample metadata to the dataframe
s <- data.frame(sample_data(Ballast_physeq))
s <- setNames(cbind(rownames(s), s, row.names = NULL), 
              c("sample-id", "Sample_type", "Voyage", "Tank", 
                "Sample_number", "Sample_site", "Sample_date"))

alphadiv <- merge(s, richness, by = "sample-id")

#Order factors
alphadiv$Sample_type <- factor(alphadiv$Sample_type, 
                      levels = c("Port uptake", "BWT", "Ocean uptake", "BWT+BWE"),
                      labels = c("Port uptake", "BWT", "Ocean uptake", "BWT+BWE")) 

#Shows the calculated indices
knitr::kable(head(alphadiv)) %>% 
  kableExtra::kable_styling("striped", latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

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

<th style="text-align:left;">

NA
</th>

<th style="text-align:left;">

NA
</th>

<th style="text-align:right;">

Observed
</th>

<th style="text-align:right;">

Shannon
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3
</td>

<td style="text-align:left;">

Port uptake
</td>

<td style="text-align:left;">

Voyage 2
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

port_uptake_2
</td>

<td style="text-align:left;">

V2_UP_T6
</td>

<td style="text-align:left;">

</td>

<td style="text-align:right;">

93
</td>

<td style="text-align:right;">

4.053663
</td>

</tr>

<tr>

<td style="text-align:left;">

Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1
</td>

<td style="text-align:left;">

Port uptake
</td>

<td style="text-align:left;">

Voyage 2
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

port_uptake_2
</td>

<td style="text-align:left;">

V2_UP_T2
</td>

<td style="text-align:left;">

</td>

<td style="text-align:right;">

76
</td>

<td style="text-align:right;">

3.984224
</td>

</tr>

<tr>

<td style="text-align:left;">

Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2
</td>

<td style="text-align:left;">

Port uptake
</td>

<td style="text-align:left;">

Voyage 2
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

port_uptake_2
</td>

<td style="text-align:left;">

V2_UP_T2
</td>

<td style="text-align:left;">

</td>

<td style="text-align:right;">

78
</td>

<td style="text-align:right;">

3.999169
</td>

</tr>

<tr>

<td style="text-align:left;">

Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3
</td>

<td style="text-align:left;">

Port uptake
</td>

<td style="text-align:left;">

Voyage 2
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

port_uptake_2
</td>

<td style="text-align:left;">

V2_UP_T2
</td>

<td style="text-align:left;">

</td>

<td style="text-align:right;">

59
</td>

<td style="text-align:right;">

3.755243
</td>

</tr>

<tr>

<td style="text-align:left;">

Ship_Ballast_110_02_05_17_SA_2_Dis6_1
</td>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:left;">

Voyage 2
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

BWT_2
</td>

<td style="text-align:left;">

V2_BWT_T6
</td>

<td style="text-align:left;">

treatment only
</td>

<td style="text-align:right;">

27
</td>

<td style="text-align:right;">

1.770741
</td>

</tr>

<tr>

<td style="text-align:left;">

Ship_Ballast_111_02_05_17_SA_2_Dis6_2
</td>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:left;">

Voyage 2
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

BWT_2
</td>

<td style="text-align:left;">

V2_BWT_T6
</td>

<td style="text-align:left;">

treatment only
</td>

<td style="text-align:right;">

24
</td>

<td style="text-align:right;">

1.681808
</td>

</tr>

</tbody>

</table>

</div>

## 3.2 Create Alpha Diversity Boxplots

### 3.2.1 Shannon diversity: sample type

``` r
#Plot
ggplot(data=alphadiv, aes(x=Sample_type, y=Shannon), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Sample_type)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(), 
        axis.title.x = element_blank())
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

![](16S-sequence-analysis_files/figure-gfm/boxplot-1.png)<!-- -->

### 3.2.2 Shannon diversity: sample type and voyage

``` r
#Plot
ggplot(data=alphadiv, aes(x=Voyage, y=Shannon), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Sample_type)) +
  facet_wrap(~Sample_type, scale="free") +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(),
        axis.title.x = element_blank())
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

![](16S-sequence-analysis_files/figure-gfm/boxplot2-1.png)<!-- -->

### 3.2.3 Shannon diversity: voyage and sample type

``` r
#Plot
ggplot(data=alphadiv, aes(x=Sample_type, y=Shannon), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Voyage)) +
  facet_wrap(~Voyage, scale="free") +
  scale_y_continuous(limits = c(0, 4.5)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

![](16S-sequence-analysis_files/figure-gfm/boxplot3-1.png)<!-- -->

### 3.2.4 Richness: sample type

``` r
#Plot
ggplot(data=alphadiv, aes(x=Sample_type, y=Observed), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Sample_type)) +
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(), 
        axis.title.x = element_blank())
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

![](16S-sequence-analysis_files/figure-gfm/boxplot%20richness%201-1.png)<!-- -->

### 3.2.5 Richness: voyage and sample type

``` r
#Plot
ggplot(data=alphadiv, aes(x=Sample_type, y=Observed), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Voyage)) +
  facet_wrap(~Voyage, scale="free") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

![](16S-sequence-analysis_files/figure-gfm/boxplot%20richness%202-1.png)<!-- -->

### 3.2.6 Reads: voyage and sample type

``` r
#Sum reads and sort by samples with least to greatest number of reads
samples_reads <- sort(sample_sums(Ballast_physeq))

samples_reads_df <- data.frame(sample_sums(Ballast_physeq))

#Make column out of the row names
samples_reads_df <- tibble::rownames_to_column(samples_reads_df, "sample-id")

#Add the sample metadata to the dataframe
s <- data.frame(sample_data(Ballast_physeq))
s <- setNames(cbind(rownames(s), s, row.names = NULL), 
              c("sample-id", "Sample_type", "Voyage", "Tank", 
                "Sample_number", "Sample_site", "Sample_date"))

reads_df <- merge(s, samples_reads_df, by = "sample-id")
```

    ## Warning in merge.data.frame(s, samples_reads_df, by = "sample-id"): column
    ## names 'NA', 'NA' are duplicated in the result

``` r
#Order factors
reads_df$Sample_type <- factor(reads_df$Sample_type, 
                      levels = c("Port uptake", "BWT", "Ocean uptake", "BWT+BWE"),
                      labels = c("Port uptake", "BWT", "Ocean uptake", "BWT+BWE"))


#Plot
reads_box <- ggplot(data=reads_df, aes(x=Sample_type, y=sample_sums.Ballast_physeq.), 
       alpha=0.1) + 
  labs(y = "Reads") +
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Voyage)) +
  facet_wrap(~Voyage, scale="free") +
  #scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

### 3.2.7 Combine plots for Fig. 2

``` r
#Observed without x-axis labels
richness_box <- ggplot(data=alphadiv, aes(x=Sample_type, y=Observed), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Voyage)) +
  facet_wrap(~Voyage, scale="free") +
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

``` r
#Shannon without x-axis labels
shannon_box <- ggplot(data=alphadiv, aes(x=Sample_type, y=Shannon), alpha=0.1) + 
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  geom_point(position=position_dodge(width=0.75),aes(group=Voyage)) +
  facet_wrap(~Voyage, scale="free") +
  scale_y_continuous(limits = c(0, 4.5)) +
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • alpha = 0.1
    ## ℹ Did you misspell an argument name?

``` r
#Combine with the sample coverage plot using patchwork
boxplot_combo <- richness_box / shannon_box / reads_box

#Save as an png file
ggsave(filename="R output/Fig2_combined_boxplot.png", 
       plot=boxplot_combo, width=5, height=8, device=png)
```

## 3.3 Alpha Diversity statistics

Perform Wilcoxon test.

Use section 7.2
<https://microbiome.github.io/course_2021_radboud/alpha-diversity.html>

Some information from the tutorial linked above: “To further investigate
if patient status could explain the variation of Shannon index, let’s do
a Wilcoxon test. This is a non-parametric test that doesn’t make
specific assumptions about the distribution, unlike popular parametric
tests, such as the t test, which assumes normally distributed
observations. Wilcoxon test can be used to estimate whether the
differences between two groups is statistically significant.”

### 3.3.1 Subset alpha diversity dataframe

``` r
#Subset the alphadiversity dataframe by voyage
v1.alpha <- subset(alphadiv, Voyage == "Voyage 1")
v2.alpha <- subset(alphadiv, Voyage == "Voyage 2")
```

### 3.3.2 Shannon diversity: sample type

Voyages are grouped together and differences in sample type are tested
using a pairwise Wilcoxon test.

``` r
sample.shan.wil <- pairwise.wilcox.test(alphadiv$Shannon, alphadiv$Sample_type, p.adjust.method="fdr", paired=FALSE) 
sample.shan.wil
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum exact test 
    ## 
    ## data:  alphadiv$Shannon and alphadiv$Sample_type 
    ## 
    ##              Port uptake BWT   Ocean uptake
    ## BWT          0.116       -     -           
    ## Ocean uptake 0.280       0.053 -           
    ## BWT+BWE      0.280       0.840 0.200       
    ## 
    ## P value adjustment method: fdr

``` r
#Convert output to dataframe
sample.shan.wil <- as.data.frame(sample.shan.wil$p.value)

#Create table of p-values
knitr::kable(sample.shan.wil) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Port uptake
</th>

<th style="text-align:right;">

BWT
</th>

<th style="text-align:right;">

Ocean uptake
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:right;">

0.1161634
</td>

<td style="text-align:right;">

NA
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

Ocean uptake
</td>

<td style="text-align:right;">

0.2795604
</td>

<td style="text-align:right;">

0.0527473
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

BWT+BWE
</td>

<td style="text-align:right;">

0.2795604
</td>

<td style="text-align:right;">

0.8395604
</td>

<td style="text-align:right;">

0.2
</td>

</tr>

</tbody>

</table>

</div>

### 3.3.3 Shannon diversity: voyage and sample type

Testing for differences in sample types in the two voyages separately
using a pairwise Wilcoxon test.

``` r
#Voyage 1
v1.sample.shan.wil <- pairwise.wilcox.test(v1.alpha$Shannon, v1.alpha$Sample_type, p.adjust.method="fdr", paired=FALSE) 
v1.sample.shan.wil
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum exact test 
    ## 
    ## data:  v1.alpha$Shannon and v1.alpha$Sample_type 
    ## 
    ##              Port uptake BWT  Ocean uptake
    ## BWT          1.00        -    -           
    ## Ocean uptake 1.00        0.29 -           
    ## BWT+BWE      1.00        1.00 0.30        
    ## 
    ## P value adjustment method: fdr

``` r
#Voyage 2
v2.sample.shan.wil <- pairwise.wilcox.test(v2.alpha$Shannon, v2.alpha$Sample_type, p.adjust.method="fdr", paired=FALSE) 
v2.sample.shan.wil
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum exact test 
    ## 
    ## data:  v2.alpha$Shannon and v2.alpha$Sample_type 
    ## 
    ##     Port uptake
    ## BWT 0.0022     
    ## 
    ## P value adjustment method: fdr

``` r
#Convert output to dataframe
v1.sample.shan.wil.df <- as.data.frame(v1.sample.shan.wil$p.value)

#Create table of p-values
knitr::kable(v1.sample.shan.wil.df) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Port uptake
</th>

<th style="text-align:right;">

BWT
</th>

<th style="text-align:right;">

Ocean uptake
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

NA
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

Ocean uptake
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

0.2857143
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

BWT+BWE
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1.0000000
</td>

<td style="text-align:right;">

0.3
</td>

</tr>

</tbody>

</table>

</div>

``` r
#Convert output to dataframe
v2.sample.shan.wil.df <- as.data.frame(v2.sample.shan.wil$p.value)

#Create table of p-values
knitr::kable(v2.sample.shan.wil.df) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Port uptake
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:right;">

0.0021645
</td>

</tr>

</tbody>

</table>

</div>

### 3.3.4 Richness: sample type

Voyages are grouped together and differences in sample type are tested
using a pairwise Wilcoxon test.

``` r
sample.richness.wil <- pairwise.wilcox.test(alphadiv$Observed, alphadiv$Sample_type, p.adjust.method="fdr", paired=FALSE) 
```

    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties

``` r
sample.richness.wil
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  alphadiv$Observed and alphadiv$Sample_type 
    ## 
    ##              Port uptake BWT   Ocean uptake
    ## BWT          0.025       -     -           
    ## Ocean uptake 0.828       0.025 -           
    ## BWT+BWE      0.025       0.025 0.092       
    ## 
    ## P value adjustment method: fdr

``` r
#Convert output to dataframe
sample.richness.wil <- as.data.frame(sample.richness.wil$p.value)

#Create table of p-values
knitr::kable(sample.richness.wil) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Port uptake
</th>

<th style="text-align:right;">

BWT
</th>

<th style="text-align:right;">

Ocean uptake
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:right;">

0.0251165
</td>

<td style="text-align:right;">

NA
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

Ocean uptake
</td>

<td style="text-align:right;">

0.8279872
</td>

<td style="text-align:right;">

0.0251165
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

BWT+BWE
</td>

<td style="text-align:right;">

0.0251165
</td>

<td style="text-align:right;">

0.0251165
</td>

<td style="text-align:right;">

0.091827
</td>

</tr>

</tbody>

</table>

</div>

### 3.3.5 Richness: voyage and sample type

Testing for differences in sample types in the two voyages separately
using a pairwise Wilcoxon test.

``` r
#Voyage 1
v1.sample.richness.wil <- pairwise.wilcox.test(v1.alpha$Observed, v1.alpha$Sample_type, p.adjust.method="fdr", paired=FALSE) 
```

    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties
    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot compute
    ## exact p-value with ties

``` r
v1.sample.richness.wil
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  v1.alpha$Observed and v1.alpha$Sample_type 
    ## 
    ##              Port uptake BWT   Ocean uptake
    ## BWT          0.747       -     -           
    ## Ocean uptake 0.185       0.083 -           
    ## BWT+BWE      0.102       0.083 0.115       
    ## 
    ## P value adjustment method: fdr

``` r
#Voyage 2
v2.sample.richness.wil <- pairwise.wilcox.test(v2.alpha$Observed, v2.alpha$Sample_type, p.adjust.method="fdr", paired=FALSE) 
v2.sample.richness.wil
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum exact test 
    ## 
    ## data:  v2.alpha$Observed and v2.alpha$Sample_type 
    ## 
    ##     Port uptake
    ## BWT 0.0022     
    ## 
    ## P value adjustment method: fdr

``` r
#Convert output to dataframe
v1.sample.richness.wil.df <- as.data.frame(v1.sample.richness.wil$p.value)

#Create table of p-values
knitr::kable(v1.sample.richness.wil.df) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Port uptake
</th>

<th style="text-align:right;">

BWT
</th>

<th style="text-align:right;">

Ocean uptake
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:right;">

0.7474911
</td>

<td style="text-align:right;">

NA
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

Ocean uptake
</td>

<td style="text-align:right;">

0.1846253
</td>

<td style="text-align:right;">

0.0825957
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

BWT+BWE
</td>

<td style="text-align:right;">

0.1016786
</td>

<td style="text-align:right;">

0.0825957
</td>

<td style="text-align:right;">

0.1147838
</td>

</tr>

</tbody>

</table>

</div>

``` r
#Convert output to dataframe
v2.sample.richness.wil.df <- as.data.frame(v2.sample.richness.wil$p.value)

#Create table of p-values
knitr::kable(v2.sample.richness.wil.df) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Port uptake
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:right;">

0.0021645
</td>

</tr>

</tbody>

</table>

</div>

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
  scale_colour_manual(values = sample_type_colors, "Sample Type") +
  scale_shape_manual(values = voyage_shapes, name = "Voyage") +
  geom_point(mapping = aes(colour = factor(Sample_type), shape = factor(Voyage), size = 5)) +
  guides(size = "none") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.justification = c("right", "top")) 

print(all.nmds.sample.type)
```

![](16S-sequence-analysis_files/figure-gfm/NMDS%20Bray-1.png)<!-- -->

## 4.3 Create PCoA with the Bray-Curtis dissimilarity matrix

### 4.3.1 PCoA with labels and arrows

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
  scale_colour_manual(values = sample_type_colors, "Sample Type") +
  scale_shape_manual(values = voyage_shapes, name = "Voyage") +
  geom_point(mapping = aes(colour = factor(Sample_type), shape = factor(Voyage), size = 5)) +
  geom_text(aes(label = Sample_number, colour = factor(Sample_type)), nudge_x=0.05, nudge_y=0.05, check_overlap = TRUE) +
  guides(size = "none") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
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

    ## Warning in geom_segment(aes(x = -0.45, y = -0.3, xend = 0.3, yend = -0.3), : All aesthetics have length 1, but the data has 30 rows.
    ## ℹ Please consider using `annotate()` or provide this layer with data containing
    ##   a single row.

    ## Warning in geom_segment(aes(x = -0.08, y = 0.03, xend = 0.29, yend = -0.08), : All aesthetics have length 1, but the data has 30 rows.
    ## ℹ Please consider using `annotate()` or provide this layer with data containing
    ##   a single row.

    ## Warning in geom_segment(aes(x = -0.14, y = 0.35, xend = 0.04, yend = 0.34), : All aesthetics have length 1, but the data has 30 rows.
    ## ℹ Please consider using `annotate()` or provide this layer with data containing
    ##   a single row.

    ## Warning in geom_curve(aes(x = -0.14, y = 0.43, xend = 0.25, yend = 0.21), : All aesthetics have length 1, but the data has 30 rows.
    ## ℹ Please consider using `annotate()` or provide this layer with data containing
    ##   a single row.

    ## Warning in geom_segment(aes(x = 0.06, y = 0.34, xend = 0.29, yend = -0.08), : All aesthetics have length 1, but the data has 30 rows.
    ## ℹ Please consider using `annotate()` or provide this layer with data containing
    ##   a single row.

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text()`).

![](16S-sequence-analysis_files/figure-gfm/PCoA%20Bray%20labels-1.png)<!-- -->

### 4.3.2 PCoA without labels and arrows

``` r
all.pcoa.source.ord <- ordinate(
  physeq = Ballast_physeqRA, 
  method = "PCoA", 
  distance = "bray"
)

#Plot, color coding by Sample_type
all.pcoa.sample.type.no.label <- plot_ordination(
  physeq = Ballast_physeqRA,
  ordination = all.pcoa.source.ord) + 
  scale_colour_manual(values = sample_type_colors, "Sample Type") +
  scale_shape_manual(values = voyage_shapes, name = "Voyage") +
  geom_point(mapping = aes(colour = factor(Sample_type), shape = factor(Voyage), size = 5)) +
  guides(size = "none") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        panel.background = element_blank(), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.justification = c("right", "top"))  
  

print(all.pcoa.sample.type.no.label)
```

![](16S-sequence-analysis_files/figure-gfm/PCoA%20Bray%20no%20labels-1.png)<!-- -->

# 5 PERMANOVA Analysis

Using the relative abundance-transformed phyloseq object here
(Ballast_physeqRA). The adonis PERMANOVA function requires the input to
be an abundance table; the output of this function will tell you if
environmental factors are influencing differences in bacterial community
composition. Useful sources with information on how to interpret
PERMANOVA results:

<https://forum.qiime2.org/t/adonis-vs-anosim-vs-permanova/9744>

<https://forum.qiime2.org/t/adonis-correct-interpretation-of-interaction-among-factors/10327>

<https://sites.google.com/site/mb3gustame/hypothesis-tests/manova/npmanova>

## 5.1 PERMANOVA analysis

``` r
#Run PERMANOVA with Bray-Curtis dissimilarity
permanova_bray <- phyloseq::distance(Ballast_physeqRA, method = "bray")
permanova_bray_df <- data.frame(sample_data(Ballast_physeqRA))
all <- adonis2(permanova_bray ~ Voyage*Sample_type, data = permanova_bray_df)
all
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = permanova_bray ~ Voyage * Sample_type, data = permanova_bray_df)
    ##          Df SumOfSqs      R2     F Pr(>F)    
    ## Model     5   9.5201 0.84915 27.02  0.001 ***
    ## Residual 24   1.6912 0.15085                 
    ## Total    29  11.2113 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Convert output to dataframe
all.pvalue <- as.data.frame(all)

#Create table of p-values
knitr::kable(all.pvalue) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Df
</th>

<th style="text-align:right;">

SumOfSqs
</th>

<th style="text-align:right;">

R2
</th>

<th style="text-align:right;">

F
</th>

<th style="text-align:right;">

Pr(\>F)
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Model
</td>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

9.520118
</td>

<td style="text-align:right;">

0.8491519
</td>

<td style="text-align:right;">

27.0201
</td>

<td style="text-align:right;">

0.001
</td>

</tr>

<tr>

<td style="text-align:left;">

Residual
</td>

<td style="text-align:right;">

24
</td>

<td style="text-align:right;">

1.691207
</td>

<td style="text-align:right;">

0.1508481
</td>

<td style="text-align:right;">

NA
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

Total
</td>

<td style="text-align:right;">

29
</td>

<td style="text-align:right;">

11.211324
</td>

<td style="text-align:right;">

1.0000000
</td>

<td style="text-align:right;">

NA
</td>

<td style="text-align:right;">

NA
</td>

</tr>

</tbody>

</table>

</div>

## 5.2 Pairwise PERMANOVA

Examine whether BWT from different tank was statistically the same

``` r
library(pairwiseAdonis)
```

    ## Loading required package: cluster

``` r
#Run pairwise PERMANOVA with bray-curtis dissimilarity
Ballast_physeq_df <- data.frame(sample_data(Ballast_physeq))

Ballast_physeq_pairwise <- phyloseq::distance(Ballast_physeq, method = "bray")

Ballast.physeq.pairwise <- pairwise.adonis(Ballast_physeq_pairwise, Ballast_physeq_df$Sample_type, p.adjust.m='BH')
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

``` r
Ballast.physeq.pairwise
```

    ##                         pairs Df SumsOfSqs    F.Model        R2 p.value
    ## 1          Port uptake vs BWT  1 2.4414627   8.875528 0.2874616   0.001
    ## 2 Port uptake vs Ocean uptake  1 1.6976855   6.522544 0.3341032   0.004
    ## 3      Port uptake vs BWT+BWE  1 1.6290724   6.216411 0.3234949   0.003
    ## 4         BWT vs Ocean uptake  1 1.8221514   8.799639 0.4036598   0.003
    ## 5              BWT vs BWT+BWE  1 0.9703614   4.646179 0.2632966   0.004
    ## 6     Ocean uptake vs BWT+BWE  1 1.2606365 107.315202 0.9640660   0.100
    ##   p.adjusted sig
    ## 1     0.0048   *
    ## 2     0.0048   *
    ## 3     0.0048   *
    ## 4     0.0048   *
    ## 5     0.0048   *
    ## 6     0.1000

# 6 Preliminary taxonomy fix

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

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

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

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

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

# 7 Bacterial community composition bar charts

## 7.1 Find color hex codes

I’ve included here a list of hex codes for colors that can be used to
identify taxa in plots of bacterial community composition; these hex
codes are visualized with the seecolor package.

``` r
bcc_hex <- c("black", "#440154FF", "#450659FF","#460B5EFF","#472D7AFF","#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF", "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF","#2E6E8EFF", "#2D718EFF", "#2B748EFF","#2A778EFF", "#297B8EFF","#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF","#1FA287FF", "#21A585FF", "#23A983FF","#25AC82FF", "#29AF7FFF", "#2DB27DFF","#32B67AFF", "#37B878FF", "#3CBC74FF","#57C766FF", "#5EC962FF","#A2DA37FF","#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF","darkorchid1", "darkorchid2", "darkorchid3","#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF","#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#FF6699","#F68D45FF", "#FCA537FF","#F6E726FF", "#F4ED27FF", "#00489C", "#CCCCCC", "#999999", "#A1C299","#300018")

print_color(bcc_hex, type = "r")
```

    ## 
    ##  ------ bcc_hex ------
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

## 7.2 Family level, individual samples

``` r
BallastSeqR_family <- Ballast_physeq %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  group_by(Sample, Kingdom, Phylum, Class, Order, Family) %>%
  filter(Abundance > 0.01) %>%
  arrange(Class)

BallastSeqR_family$Class_family <- paste(BallastSeqR_family$Class,
                                         BallastSeqR_family$Family, sep = "_")

n_family_taxa <- length(unique(BallastSeqR_family$Class_family))
family_colors <- colorRampPalette(c(bcc_hex))(n_family_taxa)

BallastSeqR_family_bar <- ggplot(BallastSeqR_family,
                                 aes(x = Sample, y = Abundance, fill = Class_family)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, position = "fill") +
  scale_fill_manual(values = family_colors) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 1%) \n") + xlab("Sample ID") +
  theme_minimal() + theme_ballast()

print(BallastSeqR_family_bar)
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20family%20sample%20id-1.png)<!-- -->

``` r
BallastSeqR_family_bar_interactive <- ggplotly(BallastSeqR_family_bar)
BallastSeqR_family_bar_interactive
```

<div id="htmlwidget-36780052264c44c40182"
class="plotly html-widget html-fill-item"
style="width:672px;height:480px;">

</div>

<script type="application/json" data-for="htmlwidget-36780052264c44c40182">{"x":{"data":[{"orientation":"v","width":0.89999999999999858,"base":0.98325820991629109,"x":[30],"y":[0.016741790083708907],"text":"Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.01674179<br />Class_family: Actinobacteria_Microbacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,0,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Microbacteriaceae","legendgroup":"Actinobacteria_Microbacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.89999999999999991,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.86477987421383651,0.9039481437831467,0.92507374631268446,0.92330040674026725,0.9315843621399178,0.92337411461687058,0.98740399385560684],"x":[2,4,1,3,29,30,13],"y":[0.13522012578616349,0.096051856216853304,0.074926253687315536,0.076699593259732746,0.068415637860082201,0.059884095299420514,0.012596006144393157],"text":["Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.13522013<br />Class_family: Actinobacteria_Sporichthyaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.09605186<br />Class_family: Actinobacteria_Sporichthyaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.07492625<br />Class_family: Actinobacteria_Sporichthyaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.07669959<br />Class_family: Actinobacteria_Sporichthyaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.06841564<br />Class_family: Actinobacteria_Sporichthyaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.05988410<br />Class_family: Actinobacteria_Sporichthyaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.01259601<br />Class_family: Actinobacteria_Sporichthyaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(68,4,87,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Sporichthyaceae","legendgroup":"Actinobacteria_Sporichthyaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.90000000000000036,0.90000000000000013],"base":[0.9053444945267225,0.90703079604880876,0.84853249475890979],"x":[30,3,2],"y":[0.01802962009014808,0.016269610691458491,0.016247379454926714],"text":["Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.01802962<br />Class_family: Alphaproteobacteria_Beijerinckiaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.01626961<br />Class_family: Alphaproteobacteria_Beijerinckiaceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.01624738<br />Class_family: Alphaproteobacteria_Beijerinckiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(70,26,106,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Beijerinckiaceae","legendgroup":"Alphaproteobacteria_Beijerinckiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.89999999999999991,0.89999999999999858],"base":[0.83176100628930816,0.88862698880377144,0.91091445427728623,0.89311010946555058],"x":[2,4,1,30],"y":[0.016771488469601636,0.015321154979375251,0.01415929203539823,0.012234385061171915],"text":["Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.01677149<br />Class_family: Alphaproteobacteria_Caulobacteraceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.01532115<br />Class_family: Alphaproteobacteria_Caulobacteraceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.01415929<br />Class_family: Alphaproteobacteria_Caulobacteraceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.01223439<br />Class_family: Alphaproteobacteria_Caulobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(58,81,139,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Caulobacteraceae","legendgroup":"Alphaproteobacteria_Caulobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858],"base":[0.97266944941363542,0.98938899253731349],"x":[20,22],"y":[0.027330550586364577,0.010611007462686506],"text":["Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.02733055<br />Class_family: Alphaproteobacteria_Clade_I","Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.01061101<br />Class_family: Alphaproteobacteria_Clade_I"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(56,88,140,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Clade_I","legendgroup":"Alphaproteobacteria_Clade_I","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999991,0.90000000000000013],"base":[0.86735350933676758,0.90637860082304522,0.89262536873156351,0.81603773584905659],"x":[30,29,1,2],"y":[0.025756600128783003,0.025205761316872577,0.018289085545722727,0.015723270440251569],"text":["Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.02575660<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.02520576<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.01828909<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.01572327<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(52,94,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","legendgroup":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999858,0.89999999999999991,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.90000000000000036],"base":[0.84796700058927521,0.83322601416613007,0.868731563421829,0.88271604938271608,0.9789892106757524,0.9824645688205621,0.98279087167976065,0.89192330040674017],"x":[4,30,1,29,26,27,28,3],"y":[0.040659988214496234,0.034127495170637512,0.023893805309734506,0.023662551440329138,0.021010789324247603,0.017535431179437899,0.017209128320239353,0.015107495642068591],"text":["Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.04065999<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.03412750<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.02389381<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.02366255<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample: Ship_Ballast_95_10_23_16_SA_1_Dis_5_1<br />Abundance: 0.02101079<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample: Ship_Ballast_96_10_23_16_SA_1_Dis_5_2<br />Abundance: 0.01753543<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample: Ship_Ballast_97_10_23_16_SA_1_Dis_5_3<br />Abundance: 0.01720913<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.01510750<br />Class_family: Alphaproteobacteria_Rhodobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(50,100,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhodobacteraceae","legendgroup":"Alphaproteobacteria_Rhodobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.90000000000000013,0.89999999999999858,0.89999999999999858,0.89999999999999991,0.90000000000000036,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.79008370895041857,0.96006144393241166,0.78878406708595383,0.97447226313205693,0.97498511018463374,0.84719764011799414,0.8710052295177223,0.82852091926929883,0.98290180752320466,0.98349120433017589,0.96894874672652442,0.96763202725724018],"x":[30,13,2,11,12,1,3,4,18,17,28,26],"y":[0.043142305215711496,0.027342549923195181,0.027253668763102756,0.025527736867943074,0.025014889815366259,0.021533923303834857,0.02091807088901787,0.019446081319976383,0.017098192476795337,0.016508795669824106,0.013842124953236223,0.011357183418512218],"text":["Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.04314231<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.02734255<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.02725367<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.02552774<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.02501489<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.02153392<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.02091807<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.01944608<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_87_10_20_16_EX_Dis_5_2<br />Abundance: 0.01709819<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_86_10_20_16_EX_Dis_5_1<br />Abundance: 0.01650880<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_97_10_23_16_SA_1_Dis_5_3<br />Abundance: 0.01384212<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample: Ship_Ballast_95_10_23_16_SA_1_Dis_5_1<br />Abundance: 0.01135718<br />Class_family: Alphaproteobacteria_Sphingomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(47,106,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Sphingomonadaceae","legendgroup":"Alphaproteobacteria_Sphingomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999991],"base":[0.84156378600823045,0.7630392788151964,0.82094395280235999],"x":[29,30,1],"y":[0.041152263374485631,0.027044430135222175,0.026253687315634155],"text":["Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.04115226<br />Class_family: Alphaproteobacteria_Tistrellaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.02704443<br />Class_family: Alphaproteobacteria_Tistrellaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.02625369<br />Class_family: Alphaproteobacteria_Tistrellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(45,112,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Tistrellaceae","legendgroup":"Alphaproteobacteria_Tistrellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.90000000000000036,0.89999999999999991,0.90000000000000036,0.90000000000000013,0.89999999999999858],"base":[0.76286008230452673,0.68898905344494532,0.80185938407902391,0.76548672566371689,0.78668238067177365,0.74580712788259951,0.95652173913043481],"x":[29,30,3,1,4,2,12],"y":[0.07870370370370372,0.074050225370251077,0.069145845438698395,0.055457227138643095,0.041838538597525177,0.042976939203354325,0.018463371054198929],"text":["Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.07870370<br />Class_family: Bacteroidia_Chitinophagaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.07405023<br />Class_family: Bacteroidia_Chitinophagaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.06914585<br />Class_family: Bacteroidia_Chitinophagaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.05545723<br />Class_family: Bacteroidia_Chitinophagaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.04183854<br />Class_family: Bacteroidia_Chitinophagaceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.04297694<br />Class_family: Bacteroidia_Chitinophagaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.01846337<br />Class_family: Bacteroidia_Chitinophagaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(42,117,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Chitinophagaceae","legendgroup":"Bacteroidia_Chitinophagaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999858,0.89999999999999991,0.90000000000000036,0.89999999999999858,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.68867924528301883,0.7119341563786008,0.71474926253687321,0.75653689715281813,0.66065679330328397,0.76075427224513847,0.93609831029185875,0.93508040500297807,0.95360824742268036],"x":[2,29,1,3,30,4,13,12,11],"y":[0.057127882599580682,0.05092592592592593,0.050737463126843685,0.045322486926205774,0.028332260141661347,0.025928108426635177,0.023963133640552914,0.021441334127456746,0.020864015709376571],"text":["Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.05712788<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.05092593<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.05073746<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.04532249<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.02833226<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.02592811<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.02396313<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.02144133<br />Class_family: Bacteroidia_Crocinitomicaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.02086402<br />Class_family: Bacteroidia_Crocinitomicaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(38,129,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Crocinitomicaceae","legendgroup":"Bacteroidia_Crocinitomicaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.89247311827956988,0.91237113402061853,0.89338892197736752],"x":[13,11,12],"y":[0.043625192012288871,0.04123711340206182,0.041691483025610543],"text":["Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.04362519<br />Class_family: Bacteroidia_Cryomorphaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.04123711<br />Class_family: Bacteroidia_Cryomorphaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.04169148<br />Class_family: Bacteroidia_Cryomorphaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,151,139,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cryomorphaceae","legendgroup":"Bacteroidia_Cryomorphaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.86881720430107534,0.87224538415723651,0.89175257731958768],"x":[13,12,11],"y":[0.023655913978494536,0.021143537820131009,0.020618556701030855],"text":["Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.02365591<br />Class_family: Bacteroidia_Cyclobacteriaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.02114354<br />Class_family: Bacteroidia_Cyclobacteriaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.02061856<br />Class_family: Bacteroidia_Cyclobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(30,158,136,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cyclobacteriaceae","legendgroup":"Bacteroidia_Cyclobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991,0.89999999999999858,0.89999999999999858,0.90000000000000036,0.90000000000000036],"base":[0.65618448637316573,0.6846607669616519,0.63168061815840315,0.68827160493827166,0.73678094131318994,0.74307601649970534],"x":[2,1,30,29,3,4],"y":[0.032494758909853094,0.030088495575221308,0.028976175144880822,0.023662551440329138,0.019755955839628192,0.017678255745433136],"text":["Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.03249476<br />Class_family: Bacteroidia_env.OPS_17","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.03008850<br />Class_family: Bacteroidia_env.OPS_17","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.02897618<br />Class_family: Bacteroidia_env.OPS_17","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.02366255<br />Class_family: Bacteroidia_env.OPS_17","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.01975596<br />Class_family: Bacteroidia_env.OPS_17","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.01767826<br />Class_family: Bacteroidia_env.OPS_17"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,163,134,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_env.OPS_17","legendgroup":"Bacteroidia_env.OPS_17","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999947,0.89999999999999858,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.90000000000000036,0.90000000000000013,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999991,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.76446700507614218,0.64212076583210598,0.76796863656163783,0.62626563430613469,0.77299548767788961,0.63563748079877125,0.85503107676801615,0.88499844382197312,0.90601813685078325,0.92116820377689945,0.92886482449589247,0.66414875072632196,0.58438155136268333,0.61985596707818924,0.56922086284610429,0.94099733536353247,0.63274336283185839,0.69416617560400706,0.95458159411647781,0.98366784730421097,0.97493003731343286,0.98840000000000006,0.97239863214460187],"x":[10,11,9,12,8,13,7,5,6,14,15,3,2,29,30,16,1,4,20,21,22,19,18],"y":[0.23553299492385782,0.2496318114874817,0.23203136343836217,0.24597974985110183,0.22700451232211039,0.23317972350230409,0.14496892323198385,0.11500155617802688,0.093981863149216749,0.078831796223100548,0.071135175504107528,0.072632190586867984,0.071802935010482405,0.068415637860082423,0.062459755312298859,0.059002664636467528,0.05191740412979351,0.048909840895698276,0.018087855297157618,0.016332152695789026,0.014458955223880632,0.011599999999999944,0.010503175378602791],"text":["Sample: Ship_Ballast_115_02_05_17_SA_2_Dis2B_3<br />Abundance: 0.23553299<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.24963181<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_114_02_05_17_SA_2_Dis2B_2<br />Abundance: 0.23203136<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.24597975<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_113_02_05_17_SA_2_Dis2B_1<br />Abundance: 0.22700451<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.23317972<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_112_02_05_17_SA_2_Dis6_3<br />Abundance: 0.14496892<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_110_02_05_17_SA_2_Dis6_1<br />Abundance: 0.11500156<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_111_02_05_17_SA_2_Dis6_2<br />Abundance: 0.09398186<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_83_10_17_16_SA_1_UP_6_1<br />Abundance: 0.07883180<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_84_10_17_16_SA_1_UP_6_2<br />Abundance: 0.07113518<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.07263219<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.07180294<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.06841564<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.06245976<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_85_10_17_16_SA_1_UP_6_3<br />Abundance: 0.05900266<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.05191740<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.04890984<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.01808786<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_90_10_20_16_EX_UP_5_2<br />Abundance: 0.01633215<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.01445896<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_88_10_20_16_EX_Dis_5_3<br />Abundance: 0.01160000<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample: Ship_Ballast_87_10_20_16_EX_Dis_5_2<br />Abundance: 0.01050318<br />Class_family: Bacteroidia_Flavobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(35,169,130,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Flavobacteriaceae","legendgroup":"Bacteroidia_Flavobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036,0.90000000000000013,0.89999999999999858,0.89999999999999991,0.89999999999999858],"base":[0.59209761766414881,0.63052445492044784,0.52672955974842772,0.57458847736625518,0.59380530973451329,0.53187379265936896],"x":[3,4,2,29,1,30],"y":[0.072051133062173145,0.063641720683559222,0.057651991614255604,0.045267489711934061,0.038938053097345104,0.037347070186735332],"text":["Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.07205113<br />Class_family: Bacteroidia_NS11-12_marine_group","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.06364172<br />Class_family: Bacteroidia_NS11-12_marine_group","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.05765199<br />Class_family: Bacteroidia_NS11-12_marine_group","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.04526749<br />Class_family: Bacteroidia_NS11-12_marine_group","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.03893805<br />Class_family: Bacteroidia_NS11-12_marine_group","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.03734707<br />Class_family: Bacteroidia_NS11-12_marine_group"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(40,174,127,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_NS11-12_marine_group","legendgroup":"Bacteroidia_NS11-12_marine_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036,0.89999999999999991],"base":[0.58809664113140836,0.55142359093550264,0.56283185840707972],"x":[4,3,1],"y":[0.042427813789039481,0.040674026728646173,0.030973451327433565],"text":["Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.04242781<br />Class_family: Bacteroidia_Sphingobacteriaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.04067403<br />Class_family: Bacteroidia_Sphingobacteriaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.03097345<br />Class_family: Bacteroidia_Sphingobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,180,123,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Sphingobacteriaceae","legendgroup":"Bacteroidia_Sphingobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000013,0.90000000000000036,0.89999999999999991,0.89999999999999858,0.89999999999999858,0.89999999999999947,0.89999999999999858,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.89999999999999858],"base":[0.43114468332364903,0.40985324947589102,0.49440188568061288,0.46961651917404129,0.48611111111111116,0.4481648422408242,0.74036792780284622,0.73805720923479012,0.83117755753401645,0.88252267106347893,0.8630563336445688,0.74289340101522849],"x":[3,2,4,1,29,30,8,9,7,6,5,10],"y":[0.12027890761185361,0.1168763102725367,0.093694755450795475,0.09321533923303843,0.088477366255144019,0.083708950418544759,0.032627559875043399,0.029911427326847706,0.023853519233999698,0.023495465787304326,0.021942110177404328,0.021573604060913687],"text":["Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.12027891<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.11687631<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.09369476<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.09321534<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.08847737<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.08370895<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_113_02_05_17_SA_2_Dis2B_1<br />Abundance: 0.03262756<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_114_02_05_17_SA_2_Dis2B_2<br />Abundance: 0.02991143<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_112_02_05_17_SA_2_Dis6_3<br />Abundance: 0.02385352<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_111_02_05_17_SA_2_Dis6_2<br />Abundance: 0.02349547<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_110_02_05_17_SA_2_Dis6_1<br />Abundance: 0.02194211<br />Class_family: Bacteroidia_Spirosomaceae","Sample: Ship_Ballast_115_02_05_17_SA_2_Dis2B_3<br />Abundance: 0.02157360<br />Class_family: Bacteroidia_Spirosomaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(56,185,118,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Spirosomaceae","legendgroup":"Bacteroidia_Spirosomaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.95452425373134331,0.93450606241303924,0.96497441951987406,0.62346588119783997],"x":[22,20,21,11],"y":[0.020405783582089554,0.020075531703438565,0.018693427784336913,0.018654884634266011],"text":["Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.02040578<br />Class_family: Cyanobacteriia_Cyanobiaceae","Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.02007553<br />Class_family: Cyanobacteriia_Cyanobiaceae","Sample: Ship_Ballast_90_10_20_16_EX_UP_5_2<br />Abundance: 0.01869343<br />Class_family: Cyanobacteriia_Cyanobiaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.01865488<br />Class_family: Cyanobacteriia_Cyanobiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(87,199,101,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Cyanobiaceae","legendgroup":"Cyanobacteriia_Cyanobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000036,"base":0.48025928108426635,"x":[4],"y":[0.014142604596346531],"text":"Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.01414260<br />Class_family: Fibrobacteria_Fibrobacteraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(146,214,64,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Fibrobacteria_Fibrobacteraceae","legendgroup":"Fibrobacteria_Fibrobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.081652074609821088,0.091299477221807326,0.11242863416776459,0.3450076804915515,0.36688505062537224,0.3703976435935199],"x":[16,15,14,13,12,11],"y":[0.85934526075371143,0.83756534727408516,0.80873956960913485,0.29062980030721974,0.25938058368076244,0.25306823760432007],"text":["Sample: Ship_Ballast_85_10_17_16_SA_1_UP_6_3<br />Abundance: 0.85934526<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Sample: Ship_Ballast_84_10_17_16_SA_1_UP_6_2<br />Abundance: 0.83756535<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Sample: Ship_Ballast_83_10_17_16_SA_1_UP_6_1<br />Abundance: 0.80873957<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.29062980<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.25938058<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.25306824<br />Class_family: Gammaproteobacteria_Aeromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(222,227,25,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Aeromonadaceae","legendgroup":"Gammaproteobacteria_Aeromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999991],"base":[0.38626834381551367,0.40848343986054614,0.43077913715389571,0.46862139917695478,0.45575221238938052],"x":[2,3,30,29,1],"y":[0.023584905660377353,0.022661243463102887,0.017385705086928493,0.017489711934156382,0.013864306784660774],"text":["Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.02358491<br />Class_family: Gammaproteobacteria_Alcaligenaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.02266124<br />Class_family: Gammaproteobacteria_Alcaligenaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.01738571<br />Class_family: Gammaproteobacteria_Alcaligenaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.01748971<br />Class_family: Gammaproteobacteria_Alcaligenaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.01386431<br />Class_family: Gammaproteobacteria_Alcaligenaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(238,229,27,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alcaligenaceae","legendgroup":"Gammaproteobacteria_Alcaligenaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.65962506994963621,0.70009768804949524,0.71738578680203047,0.87932992617830774,0.91256305548883021,0.91058735503179944,0.94316644113667125,0.95300000000000007,0.94137762579384465,0.33621203097081598,0.3546882670594011],"x":[25,23,24,26,27,28,17,19,18,12,11],"y":[0.34037493005036379,0.29990231195050476,0.28261421319796953,0.088302101078932438,0.06990151333173189,0.058361391694724984,0.040324763193504642,0.035399999999999987,0.031021006350757219,0.030673019654556266,0.015709376534118802],"text":["Sample: Ship_Ballast_94_10_23_16_SA_1_Dis_6_3<br />Abundance: 0.34037493<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_92_10_23_16_SA_1_Dis_6_1<br />Abundance: 0.29990231<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_93_10_23_16_SA_1_Dis_6_2<br />Abundance: 0.28261421<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_95_10_23_16_SA_1_Dis_5_1<br />Abundance: 0.08830210<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_96_10_23_16_SA_1_Dis_5_2<br />Abundance: 0.06990151<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_97_10_23_16_SA_1_Dis_5_3<br />Abundance: 0.05836139<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_86_10_20_16_EX_Dis_5_1<br />Abundance: 0.04032476<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_88_10_20_16_EX_Dis_5_3<br />Abundance: 0.03540000<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_87_10_20_16_EX_Dis_5_2<br />Abundance: 0.03102101<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.03067302<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.01570938<br />Class_family: Gammaproteobacteria_Alteromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(193,70,243,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alteromonadaceae","legendgroup":"Gammaproteobacteria_Alteromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999858,0.89999999999999991,0.90000000000000036,0.90000000000000036,0.89999999999999858],"base":[0.32442348008385746,0.41923868312757201,0.40766961651917405,0.36083672283556067,0.43606364172068357,0.39343206696716038],"x":[2,29,1,3,4,30],"y":[0.061844863731656208,0.049382716049382769,0.048082595870206468,0.047646717024985463,0.044195639363582784,0.037347070186735332],"text":["Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.06184486<br />Class_family: Gammaproteobacteria_Burkholderiaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.04938272<br />Class_family: Gammaproteobacteria_Burkholderiaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.04808260<br />Class_family: Gammaproteobacteria_Burkholderiaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.04764672<br />Class_family: Gammaproteobacteria_Burkholderiaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.04419564<br />Class_family: Gammaproteobacteria_Burkholderiaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.03734707<br />Class_family: Gammaproteobacteria_Burkholderiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(161,52,215,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Burkholderiaceae","legendgroup":"Gammaproteobacteria_Burkholderiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.31797235023041476,0.3117927337701013,0.33136966126656847],"x":[13,12,11],"y":[0.027035330261136747,0.024419297200714674,0.023318605792832625],"text":["Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.02703533<br />Class_family: Gammaproteobacteria_Chromobacteriaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.02441930<br />Class_family: Gammaproteobacteria_Chromobacteriaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.02331861<br />Class_family: Gammaproteobacteria_Chromobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(163,29,153,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Chromobacteriaceae","legendgroup":"Gammaproteobacteria_Chromobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.8519283746556473,0.84853078358208966,0.8363148479427549],"x":[21,22,20],"y":[0.11304604486422676,0.10599347014925364,0.098191214470284338],"text":["Sample: Ship_Ballast_90_10_20_16_EX_UP_5_2<br />Abundance: 0.11304604<br />Class_family: Gammaproteobacteria_Colwelliaceae","Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.10599347<br />Class_family: Gammaproteobacteria_Colwelliaceae","Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.09819121<br />Class_family: Gammaproteobacteria_Colwelliaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(171,35,147,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Colwelliaceae","legendgroup":"Gammaproteobacteria_Colwelliaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.89999999999999858,0.89999999999999858,0.90000000000000036,0.90000000000000036,0.90000000000000013,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.89999999999999858,0.89999999999999858],"base":[0.11150442477876109,0.15020576131687244,0.13650998068254991,0.20329994107248087,0.14294015107495642,0.11949685534591195,0.85142083897158316,0.87560000000000004,0.25527736867943052,0.87933561309233021,0.24627754615842765,0.25591397849462366,0.83286647992530338,0.85346248969497107,0.80446833529312956,0.7238578680203045,0.72658632205604756],"x":[1,29,30,4,3,2,17,19,11,18,12,13,5,6,7,10,9],"y":[0.29616519174041295,0.2690329218106996,0.25692208628461044,0.2327637006482027,0.21789657176060426,0.20492662473794551,0.091745602165088092,0.077400000000000024,0.076092292587137955,0.062042012701514437,0.065515187611673659,0.062058371735791096,0.030189853719265414,0.029060181368507854,0.026709222240886898,0.019035532994923998,0.011470887178742561],"text":["Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.29616519<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.26903292<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.25692209<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.23276370<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.21789657<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.20492662<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_86_10_20_16_EX_Dis_5_1<br />Abundance: 0.09174560<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_88_10_20_16_EX_Dis_5_3<br />Abundance: 0.07740000<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.07609229<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_87_10_20_16_EX_Dis_5_2<br />Abundance: 0.06204201<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.06551519<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.06205837<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_110_02_05_17_SA_2_Dis6_1<br />Abundance: 0.03018985<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_111_02_05_17_SA_2_Dis6_2<br />Abundance: 0.02906018<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_112_02_05_17_SA_2_Dis6_3<br />Abundance: 0.02670922<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_115_02_05_17_SA_2_Dis2B_3<br />Abundance: 0.01903553<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample: Ship_Ballast_114_02_05_17_SA_2_Dis2B_2<br />Abundance: 0.01147089<br />Class_family: Gammaproteobacteria_Comamonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(177,42,143,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Comamonadaceae","legendgroup":"Gammaproteobacteria_Comamonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036,0.89999999999999858,0.90000000000000036,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999947],"base":[0.14908662345315263,0.092969203951191168,0.075811209439528029,0.084905660377358486,0.77641525281370738,0.12037037037037038,0.82811211871393242,0.81061313414254588,0.11461687057308435,0.70761421319796958,0.71061420066792502,0.72544255466851781],"x":[4,3,1,2,7,29,6,5,30,10,9,8],"y":[0.054213317619328238,0.04997094712376525,0.035693215339233059,0.034591194968553465,0.028053082479422176,0.029835390946502061,0.02535037098103865,0.022253345782757505,0.021893110109465555,0.016243654822334919,0.01597212138812254,0.014925373134328401],"text":["Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.05421332<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.04997095<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.03569322<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.03459119<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_112_02_05_17_SA_2_Dis6_3<br />Abundance: 0.02805308<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.02983539<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_111_02_05_17_SA_2_Dis6_2<br />Abundance: 0.02535037<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_110_02_05_17_SA_2_Dis6_1<br />Abundance: 0.02225335<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.02189311<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_115_02_05_17_SA_2_Dis2B_3<br />Abundance: 0.01624365<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_114_02_05_17_SA_2_Dis2B_2<br />Abundance: 0.01597212<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample: Ship_Ballast_113_02_05_17_SA_2_Dis2B_1<br />Abundance: 0.01492537<br />Class_family: Gammaproteobacteria_Methylophilaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(184,49,137,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylophilaceae","legendgroup":"Gammaproteobacteria_Methylophilaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.7696629213483146,0.74451962110960757,0.77100000000000002,0.20004909180166913,0.19029184038117927,0.20614439324116746],"x":[18,17,19,11,12,13],"y":[0.10967269174401562,0.10690121786197559,0.10460000000000003,0.055228276877761384,0.055985705777248373,0.049769585253456206],"text":["Sample: Ship_Ballast_87_10_20_16_EX_Dis_5_2<br />Abundance: 0.10967269<br />Class_family: Gammaproteobacteria_Moraxellaceae","Sample: Ship_Ballast_86_10_20_16_EX_Dis_5_1<br />Abundance: 0.10690122<br />Class_family: Gammaproteobacteria_Moraxellaceae","Sample: Ship_Ballast_88_10_20_16_EX_Dis_5_3<br />Abundance: 0.10460000<br />Class_family: Gammaproteobacteria_Moraxellaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.05522828<br />Class_family: Gammaproteobacteria_Moraxellaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.05598571<br />Class_family: Gammaproteobacteria_Moraxellaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.04976959<br />Class_family: Gammaproteobacteria_Moraxellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(190,56,132,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Moraxellaceae","legendgroup":"Gammaproteobacteria_Moraxellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999947,0.89999999999999858,0.90000000000000036,0.90000000000000036,0.90000000000000036,0.89999999999999991,0.90000000000000013,0.89999999999999858,0.90000000000000036,0.90000000000000036],"base":[0.49368375199651515,0.50798333911836169,0.50228426395939085,0.64639677473542745,0.70115416323165702,0.68689698101462804,0.03716814159292036,0.051362683438155136,0.093621399176954737,0.12669416617560403,0.080185938407902377],"x":[9,8,10,7,6,5,1,2,29,4,3],"y":[0.21693044867140987,0.21745921555015613,0.20532994923857872,0.13001847807827993,0.1269579554822754,0.12371615312791784,0.038643067846607669,0.033542976939203349,0.026748971193415641,0.0223924572775486,0.012783265543288791],"text":["Sample: Ship_Ballast_114_02_05_17_SA_2_Dis2B_2<br />Abundance: 0.21693045<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_113_02_05_17_SA_2_Dis2B_1<br />Abundance: 0.21745922<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_115_02_05_17_SA_2_Dis2B_3<br />Abundance: 0.20532995<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_112_02_05_17_SA_2_Dis6_3<br />Abundance: 0.13001848<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_111_02_05_17_SA_2_Dis6_2<br />Abundance: 0.12695796<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_110_02_05_17_SA_2_Dis6_1<br />Abundance: 0.12371615<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.03864307<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.03354298<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.02674897<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.02239246<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.01278327<br />Class_family: Gammaproteobacteria_Oxalobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(197,63,125,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Oxalobacteraceae","legendgroup":"Gammaproteobacteria_Oxalobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.47585875570502045,0.44946053378762063,0.48821548821548821,0.70534048507462688,0.7157622739018088,0.73278236914600547],"x":[27,26,28,22,20,21],"y":[0.43670429978380976,0.42986939239068711,0.42237186681631123,0.14319029850746279,0.1205525740409461,0.11914600550964183],"text":["Sample: Ship_Ballast_96_10_23_16_SA_1_Dis_5_2<br />Abundance: 0.43670430<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","Sample: Ship_Ballast_95_10_23_16_SA_1_Dis_5_1<br />Abundance: 0.42986939<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","Sample: Ship_Ballast_97_10_23_16_SA_1_Dis_5_3<br />Abundance: 0.42237187<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.14319030<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.12055257<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","Sample: Ship_Ballast_90_10_20_16_EX_UP_5_2<br />Abundance: 0.11914601<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(203,69,120,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudoalteromonadaceae","legendgroup":"Gammaproteobacteria_Pseudoalteromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.90000000000000036,0.90000000000000036,0.89999999999999858,0.90000000000000036,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999947,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0,0,0.064086294416243653,0,0.065125366330185605,0.063933967543368775,0,0,0,0.096894874672652453,0.11602209944751381,0.10902896081771721,0.025282277859597448,0.034408602150537641,0.028290649195949973,0.015810276679841896,0.014563106796116504,0.017320137038446898],"x":[18,19,17,6,5,24,7,23,25,10,8,9,28,27,26,11,13,12,14,15,16],"y":[0.7696629213483146,0.77100000000000002,0.74451962110960757,0.70115416323165702,0.68689698101462804,0.65329949238578677,0.64639677473542745,0.63497232171930962,0.59569110240626744,0.50228426395939085,0.50798333911836169,0.49368375199651515,0.39132061354283576,0.35983665625750666,0.34043157296990345,0.17476681394207169,0.17173579109062981,0.16200119118522929,0.09661835748792269,0.076736370425690828,0.06433193757137419],"text":["Sample: Ship_Ballast_87_10_20_16_EX_Dis_5_2<br />Abundance: 0.76966292<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_88_10_20_16_EX_Dis_5_3<br />Abundance: 0.77100000<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_86_10_20_16_EX_Dis_5_1<br />Abundance: 0.74451962<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_111_02_05_17_SA_2_Dis6_2<br />Abundance: 0.70115416<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_110_02_05_17_SA_2_Dis6_1<br />Abundance: 0.68689698<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_93_10_23_16_SA_1_Dis_6_2<br />Abundance: 0.65329949<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_112_02_05_17_SA_2_Dis6_3<br />Abundance: 0.64639677<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_92_10_23_16_SA_1_Dis_6_1<br />Abundance: 0.63497232<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_94_10_23_16_SA_1_Dis_6_3<br />Abundance: 0.59569110<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_115_02_05_17_SA_2_Dis2B_3<br />Abundance: 0.50228426<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_113_02_05_17_SA_2_Dis2B_1<br />Abundance: 0.50798334<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_114_02_05_17_SA_2_Dis2B_2<br />Abundance: 0.49368375<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_97_10_23_16_SA_1_Dis_5_3<br />Abundance: 0.39132061<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_96_10_23_16_SA_1_Dis_5_2<br />Abundance: 0.35983666<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_95_10_23_16_SA_1_Dis_5_1<br />Abundance: 0.34043157<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.17476681<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.17173579<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.16200119<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_83_10_17_16_SA_1_UP_6_1<br />Abundance: 0.09661836<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_84_10_17_16_SA_1_UP_6_2<br />Abundance: 0.07673637<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample: Ship_Ballast_85_10_17_16_SA_1_UP_6_3<br />Abundance: 0.06433194<br />Class_family: Gammaproteobacteria_Pseudomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(207,76,115,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudomonadaceae","legendgroup":"Gammaproteobacteria_Pseudomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0.49294374875770225,0.48787313432835822,0.51554506099960651,0.036512130674993995,0.036342986939239069,0.034044145155256271],"x":[20,22,21,27,26,28],"y":[0.22281852514410655,0.21746735074626866,0.21723730814639897,0.079509968772519812,0.072685973878478138,0.062850729517396176],"text":["Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.22281853<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.21746735<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","Sample: Ship_Ballast_90_10_20_16_EX_UP_5_2<br />Abundance: 0.21723731<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","Sample: Ship_Ballast_96_10_23_16_SA_1_Dis_5_2<br />Abundance: 0.07950997<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","Sample: Ship_Ballast_95_10_23_16_SA_1_Dis_5_1<br />Abundance: 0.07268597<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","Sample: Ship_Ballast_97_10_23_16_SA_1_Dis_5_3<br />Abundance: 0.06285073<br />Class_family: Gammaproteobacteria_Saccharospirillaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(213,83,110,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Saccharospirillaceae","legendgroup":"Gammaproteobacteria_Saccharospirillaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0,0,0,0,0.45918843283582089,0.48740653286107827,0.4694891671635858,0,0,0,0.01935483870967742,0.011291114383897889],"x":[23,25,24,26,27,28,22,21,20,16,14,15,13,11],"y":[0.065125366330185605,0.063933967543368775,0.064086294416243653,0.036342986939239069,0.036512130674993995,0.034044145155256271,0.028684701492537323,0.02813852813852824,0.023454581594116453,0.017320137038446898,0.015810276679841896,0.014563106796116504,0.015053763440860221,0.013991163475699559],"text":["Sample: Ship_Ballast_92_10_23_16_SA_1_Dis_6_1<br />Abundance: 0.06512537<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_94_10_23_16_SA_1_Dis_6_3<br />Abundance: 0.06393397<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_93_10_23_16_SA_1_Dis_6_2<br />Abundance: 0.06408629<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_95_10_23_16_SA_1_Dis_5_1<br />Abundance: 0.03634299<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_96_10_23_16_SA_1_Dis_5_2<br />Abundance: 0.03651213<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_97_10_23_16_SA_1_Dis_5_3<br />Abundance: 0.03404415<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.02868470<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_90_10_20_16_EX_UP_5_2<br />Abundance: 0.02813853<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.02345458<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_85_10_17_16_SA_1_UP_6_3<br />Abundance: 0.01732014<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_83_10_17_16_SA_1_UP_6_1<br />Abundance: 0.01581028<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_84_10_17_16_SA_1_UP_6_2<br />Abundance: 0.01456311<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.01505376<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.01399116<br />Class_family: Gammaproteobacteria_Shewanellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(217,90,105,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Shewanellaceae","legendgroup":"Gammaproteobacteria_Shewanellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0,0],"x":[21,20,22],"y":[0.48740653286107827,0.4694891671635858,0.45918843283582089],"text":["Sample: Ship_Ballast_90_10_20_16_EX_UP_5_2<br />Abundance: 0.48740653<br />Class_family: Gammaproteobacteria_Vibrionaceae","Sample: Ship_Ballast_89_10_20_16_EX_UP_5_1<br />Abundance: 0.46948917<br />Class_family: Gammaproteobacteria_Vibrionaceae","Sample: Ship_Ballast_91_10_20_16_EX_UP_5_3<br />Abundance: 0.45918843<br />Class_family: Gammaproteobacteria_Vibrionaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(222,97,100,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Vibrionaceae","legendgroup":"Gammaproteobacteria_Vibrionaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.90000000000000036,0.90000000000000036,0.89999999999999858,0.90000000000000013,0.89999999999999991],"base":[0.081777205408886028,0.097230406599882149,0.054619407321324816,0.072530864197530867,0.030922431865828093,0.01858407079646018],"x":[30,4,3,29,2,1],"y":[0.032839665164198326,0.029463759575721879,0.025566531086577561,0.021090534979423869,0.020440251572327043,0.01858407079646018],"text":["Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.03283967<br />Class_family: Planctomycetes_Isosphaeraceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.02946376<br />Class_family: Planctomycetes_Isosphaeraceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.02556653<br />Class_family: Planctomycetes_Isosphaeraceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.02109053<br />Class_family: Planctomycetes_Isosphaeraceae","Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.02044025<br />Class_family: Planctomycetes_Isosphaeraceae","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.01858407<br />Class_family: Planctomycetes_Isosphaeraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(251,116,121,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Planctomycetes_Isosphaeraceae","legendgroup":"Planctomycetes_Isosphaeraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.89999999999999858,0.89999999999999858],"base":[0,0.011911852293031567,0],"x":[13,12,11],"y":[0.01935483870967742,0.016378796902918408,0.011291114383897889],"text":["Sample: Ship_Ballast_82_10_18_16_SA_1_UP_5_3<br />Abundance: 0.01935484<br />Class_family: Planctomycetes_Rubinisphaeraceae","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.01637880<br />Class_family: Planctomycetes_Rubinisphaeraceae","Sample: Ship_Ballast_80_10_18_16_SA_1_UP_5_1<br />Abundance: 0.01129111<br />Class_family: Planctomycetes_Rubinisphaeraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(251,171,53,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Planctomycetes_Rubinisphaeraceae","legendgroup":"Planctomycetes_Rubinisphaeraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999858,0.90000000000000036,0.89999999999999858],"base":[0.054526748971193424,0.081909251620506787,0.06889890534449454],"x":[29,4,30],"y":[0.018004115226337443,0.015321154979375362,0.012878300064391487],"text":["Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.01800412<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.01532115<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.01287830<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(244,235,38,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Chthoniobacteraceae","legendgroup":"Verrucomicrobiae_Chthoniobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000036,"base":0.068945197407189157,"x":[4],"y":[0.01296405421331763],"text":"Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.01296405<br />Class_family: Verrucomicrobiae_Rubritaleaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(112,144,182,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Rubritaleaceae","legendgroup":"Verrucomicrobiae_Rubritaleaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999858,0.90000000000000036,0.89999999999999858,0.89999999999999991,0.89999999999999858],"base":[0.018856806128461991,0.028332260141661302,0.019174898314933179,0.024176954732510289,0,0],"x":[4,30,3,29,1,12],"y":[0.050088391278727162,0.040566645202833235,0.035444509006391636,0.030349794238683135,0.01858407079646018,0.011911852293031567],"text":["Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.05008839<br />Class_family: Verrucomicrobiae_uncultured","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.04056665<br />Class_family: Verrucomicrobiae_uncultured","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.03544451<br />Class_family: Verrucomicrobiae_uncultured","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.03034979<br />Class_family: Verrucomicrobiae_uncultured","Sample: Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3<br />Abundance: 0.01858407<br />Class_family: Verrucomicrobiae_uncultured","Sample: Ship_Ballast_81_10_18_16_SA_1_UP_5_2<br />Abundance: 0.01191185<br />Class_family: Verrucomicrobiae_uncultured"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(155,164,153,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_uncultured","legendgroup":"Verrucomicrobiae_uncultured","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999858,0.89999999999999858,0.90000000000000036,0.90000000000000036],"base":[0,0,0,0,0],"x":[2,30,29,4,3],"y":[0.030922431865828093,0.028332260141661302,0.024176954732510289,0.018856806128461991,0.019174898314933179],"text":["Sample: Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1<br />Abundance: 0.03092243<br />Class_family: Verrucomicrobiae_Verrucomicrobiaceae","Sample: Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2<br />Abundance: 0.02833226<br />Class_family: Verrucomicrobiae_Verrucomicrobiaceae","Sample: Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1<br />Abundance: 0.02417695<br />Class_family: Verrucomicrobiae_Verrucomicrobiaceae","Sample: Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3<br />Abundance: 0.01885681<br />Class_family: Verrucomicrobiae_Verrucomicrobiaceae","Sample: Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2<br />Abundance: 0.01917490<br />Class_family: Verrucomicrobiae_Verrucomicrobiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Verrucomicrobiaceae","legendgroup":"Verrucomicrobiae_Verrucomicrobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":187.54763495747218,"l":53.466168534661698},"paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,30.600000000000001],"tickmode":"array","ticktext":["Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3","Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1","Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2","Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3","Ship_Ballast_110_02_05_17_SA_2_Dis6_1","Ship_Ballast_111_02_05_17_SA_2_Dis6_2","Ship_Ballast_112_02_05_17_SA_2_Dis6_3","Ship_Ballast_113_02_05_17_SA_2_Dis2B_1","Ship_Ballast_114_02_05_17_SA_2_Dis2B_2","Ship_Ballast_115_02_05_17_SA_2_Dis2B_3","Ship_Ballast_80_10_18_16_SA_1_UP_5_1","Ship_Ballast_81_10_18_16_SA_1_UP_5_2","Ship_Ballast_82_10_18_16_SA_1_UP_5_3","Ship_Ballast_83_10_17_16_SA_1_UP_6_1","Ship_Ballast_84_10_17_16_SA_1_UP_6_2","Ship_Ballast_85_10_17_16_SA_1_UP_6_3","Ship_Ballast_86_10_20_16_EX_Dis_5_1","Ship_Ballast_87_10_20_16_EX_Dis_5_2","Ship_Ballast_88_10_20_16_EX_Dis_5_3","Ship_Ballast_89_10_20_16_EX_UP_5_1","Ship_Ballast_90_10_20_16_EX_UP_5_2","Ship_Ballast_91_10_20_16_EX_UP_5_3","Ship_Ballast_92_10_23_16_SA_1_Dis_6_1","Ship_Ballast_93_10_23_16_SA_1_Dis_6_2","Ship_Ballast_94_10_23_16_SA_1_Dis_6_3","Ship_Ballast_95_10_23_16_SA_1_Dis_5_1","Ship_Ballast_96_10_23_16_SA_1_Dis_5_2","Ship_Ballast_97_10_23_16_SA_1_Dis_5_3","Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1","Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2"],"tickvals":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30],"categoryorder":"array","categoryarray":["Ship_Ballast_100_02_01_17_SA_2_tank_6_uptake_3","Ship_Ballast_101_02_01_17_SA_2_tank_2_uptake_1","Ship_Ballast_102_02_01_17_SA_2_tank_2_uptake_2","Ship_Ballast_103_02_01_17_SA_2_tank_2_uptake_3","Ship_Ballast_110_02_05_17_SA_2_Dis6_1","Ship_Ballast_111_02_05_17_SA_2_Dis6_2","Ship_Ballast_112_02_05_17_SA_2_Dis6_3","Ship_Ballast_113_02_05_17_SA_2_Dis2B_1","Ship_Ballast_114_02_05_17_SA_2_Dis2B_2","Ship_Ballast_115_02_05_17_SA_2_Dis2B_3","Ship_Ballast_80_10_18_16_SA_1_UP_5_1","Ship_Ballast_81_10_18_16_SA_1_UP_5_2","Ship_Ballast_82_10_18_16_SA_1_UP_5_3","Ship_Ballast_83_10_17_16_SA_1_UP_6_1","Ship_Ballast_84_10_17_16_SA_1_UP_6_2","Ship_Ballast_85_10_17_16_SA_1_UP_6_3","Ship_Ballast_86_10_20_16_EX_Dis_5_1","Ship_Ballast_87_10_20_16_EX_Dis_5_2","Ship_Ballast_88_10_20_16_EX_Dis_5_3","Ship_Ballast_89_10_20_16_EX_UP_5_1","Ship_Ballast_90_10_20_16_EX_UP_5_2","Ship_Ballast_91_10_20_16_EX_UP_5_3","Ship_Ballast_92_10_23_16_SA_1_Dis_6_1","Ship_Ballast_93_10_23_16_SA_1_Dis_6_2","Ship_Ballast_94_10_23_16_SA_1_Dis_6_3","Ship_Ballast_95_10_23_16_SA_1_Dis_5_1","Ship_Ballast_96_10_23_16_SA_1_Dis_5_2","Ship_Ballast_97_10_23_16_SA_1_Dis_5_3","Ship_Ballast_98_02_01_17_SA_2_tank_6_uptake_1","Ship_Ballast_99_02_01_17_SA_2_tank_6_uptake_2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"Sample ID","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.13283520132835203,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835205},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"Relative Abundance (Family > 1%) <br />","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"title":{"text":"Class_family","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"1d386e043767":{"x":{},"y":{},"fill":{},"type":"bar"}},"cur_data":"1d386e043767","visdat":{"1d386e043767":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

## 7.3 Family level, grouped by sample type

``` r
BallastSeqR_family_sample <- Ballast_physeq %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  group_by(Sample_type, Kingdom, Phylum, Class, Order, Family) %>%
  dplyr::summarize(Mean = mean(Abundance, na.rm = TRUE), .groups = "keep") %>%
  filter(Mean > 0.01) %>%
  arrange(Class)

BallastSeqR_family_sample$Class_family <- paste(BallastSeqR_family_sample$Class,
                                                BallastSeqR_family_sample$Family, sep = "_")

n_family_sample_taxa <- length(unique(BallastSeqR_family_sample$Class_family))
family_sample_colors <- colorRampPalette(c(bcc_hex))(n_family_sample_taxa)

BallastSeqR_sample_family_bar <- ggplot(BallastSeqR_family_sample,
                                        aes(x = Sample_type, y = Mean, fill = Class_family)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, position = "fill") +
  scale_fill_manual(values = family_sample_colors) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 1%) \n") + xlab("Sample Type") +
  theme_minimal() + theme_ballast()

print(BallastSeqR_sample_family_bar)
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20family%20sample%20type-1.png)<!-- -->

``` r
BallastSeqR_sample_family_bar_interactive <- ggplotly(BallastSeqR_sample_family_bar)
BallastSeqR_sample_family_bar_interactive
```

<div id="htmlwidget-ea957a0f5dc9f2692740"
class="plotly html-widget html-fill-item"
style="width:672px;height:480px;">

</div>

<script type="application/json" data-for="htmlwidget-ea957a0f5dc9f2692740">{"x":{"data":[{"orientation":"v","width":0.89999999999999991,"base":0.95191228951387408,"x":[1],"y":[0.048087710486125923],"text":"Sample_type: Port uptake<br />Mean: 0.04808771<br />Class_family: Actinobacteria_Sporichthyaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,0,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Sporichthyaceae","legendgroup":"Actinobacteria_Sporichthyaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000036,"base":0.98393148707175926,"x":[3],"y":[0.016068512928240741],"text":"Sample_type: Ocean uptake<br />Mean: 0.01606851<br />Class_family: Alphaproteobacteria_Clade_I","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(69,9,92,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Clade_I","legendgroup":"Alphaproteobacteria_Clade_I","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000036],"base":[0.93659072346700623,0.98146259887823217],"x":[1,4],"y":[0.015321566046867852,0.018537401121767827],"text":["Sample_type: Port uptake<br />Mean: 0.01532157<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Sample_type: BWT+BWE<br />Mean: 0.01853740<br />Class_family: Alphaproteobacteria_Rhodobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(58,82,139,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhodobacteraceae","legendgroup":"Alphaproteobacteria_Rhodobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000036],"base":[0.91551278621031129,0.97035754985952827],"x":[1,4],"y":[0.021077937256694934,0.011105049018703905],"text":["Sample_type: Port uptake<br />Mean: 0.02107794<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Sample_type: BWT+BWE<br />Mean: 0.01110505<br />Class_family: Alphaproteobacteria_Sphingomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(53,93,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Sphingomonadaceae","legendgroup":"Alphaproteobacteria_Sphingomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.87954507066954601,"x":[1],"y":[0.035967715540765277],"text":"Sample_type: Port uptake<br />Mean: 0.03596772<br />Class_family: Bacteroidia_Chitinophagaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,103,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Chitinophagaceae","legendgroup":"Bacteroidia_Chitinophagaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.8506869414683168,"x":[1],"y":[0.02885812920122921],"text":"Sample_type: Port uptake<br />Mean: 0.02885813<br />Class_family: Bacteroidia_Crocinitomicaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(45,112,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Crocinitomicaceae","legendgroup":"Bacteroidia_Crocinitomicaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.834942682315305,"x":[1],"y":[0.015744259153011808],"text":"Sample_type: Port uptake<br />Mean: 0.01574426<br />Class_family: Bacteroidia_env.OPS_17","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(41,121,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_env.OPS_17","legendgroup":"Bacteroidia_env.OPS_17","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[0.71725101949574266,0.90993475730050999,0.96769512599305452],"x":[1,2,3],"y":[0.11769166281956234,0.090065242699490011,0.016236361078704742],"text":["Sample_type: Port uptake<br />Mean: 0.11769166<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample_type: BWT<br />Mean: 0.09006524<br />Class_family: Bacteroidia_Flavobacteriaceae","Sample_type: Ocean uptake<br />Mean: 0.01623636<br />Class_family: Bacteroidia_Flavobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(30,153,138,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Flavobacteriaceae","legendgroup":"Bacteroidia_Flavobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.68863102682742838,"x":[1],"y":[0.028619992668314276],"text":"Sample_type: Port uptake<br />Mean: 0.02861999<br />Class_family: Bacteroidia_NS11-12_marine_group","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,162,134,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_NS11-12_marine_group","legendgroup":"Bacteroidia_NS11-12_marine_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.63322834703267106,0.89670463141553536],"x":[1,2],"y":[0.055402679794757326,0.013230125884974631],"text":["Sample_type: Port uptake<br />Mean: 0.05540268<br />Class_family: Bacteroidia_Spirosomaceae","Sample_type: BWT<br />Mean: 0.01323013<br />Class_family: Bacteroidia_Spirosomaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(36,171,130,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Spirosomaceae","legendgroup":"Bacteroidia_Spirosomaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000036,"base":0.94803014017525467,"x":[3],"y":[0.019664985817799852],"text":"Sample_type: Ocean uptake<br />Mean: 0.01966499<br />Class_family: Cyanobacteriia_Cyanobiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,180,123,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Cyanobiaceae","legendgroup":"Cyanobacteriia_Cyanobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.32642955285480008,"x":[1],"y":[0.30679879417787098],"text":"Sample_type: Port uptake<br />Mean: 0.30679879<br />Class_family: Gammaproteobacteria_Aeromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(69,191,110,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Aeromonadaceae","legendgroup":"Gammaproteobacteria_Aeromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036],"base":[0.80884102954926973,0.89836879895678778],"x":[2,4],"y":[0.087863601866265628,0.071988750902740484],"text":["Sample_type: BWT<br />Mean: 0.08786360<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Sample_type: BWT+BWE<br />Mean: 0.07198875<br />Class_family: Gammaproteobacteria_Alteromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(168,219,51,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alteromonadaceae","legendgroup":"Gammaproteobacteria_Alteromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.30016709802667979,"x":[1],"y":[0.026262454828120285],"text":"Sample_type: Port uptake<br />Mean: 0.02626245<br />Class_family: Gammaproteobacteria_Burkholderiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(235,228,26,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Burkholderiaceae","legendgroup":"Gammaproteobacteria_Burkholderiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000036,"base":0.84269002092840517,"x":[3],"y":[0.1053401192468495],"text":"Sample_type: Ocean uptake<br />Mean: 0.10534012<br />Class_family: Gammaproteobacteria_Colwelliaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(182,59,244,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Colwelliaceae","legendgroup":"Gammaproteobacteria_Colwelliaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.1495943632433101,0.77912403976480993],"x":[1,2],"y":[0.15057273478336969,0.029716989784459802],"text":["Sample_type: Port uptake<br />Mean: 0.15057273<br />Class_family: Gammaproteobacteria_Comamonadaceae","Sample_type: BWT<br />Mean: 0.02971699<br />Class_family: Gammaproteobacteria_Comamonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(163,29,153,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Comamonadaceae","legendgroup":"Gammaproteobacteria_Comamonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.12951330118815496,0.76887714110232086],"x":[1,2],"y":[0.020081062055155141,0.010246898662489068],"text":["Sample_type: Port uptake<br />Mean: 0.02008106<br />Class_family: Gammaproteobacteria_Methylophilaceae","Sample_type: BWT<br />Mean: 0.01024690<br />Class_family: Gammaproteobacteria_Methylophilaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(174,40,145,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylophilaceae","legendgroup":"Gammaproteobacteria_Methylophilaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.11283520717950167,0.74112283835449111],"x":[1,2],"y":[0.016678094008653291,0.027754302747829751],"text":["Sample_type: Port uptake<br />Mean: 0.01667809<br />Class_family: Gammaproteobacteria_Moraxellaceae","Sample_type: BWT<br />Mean: 0.02775430<br />Class_family: Gammaproteobacteria_Moraxellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(185,50,136,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Moraxellaceae","legendgroup":"Gammaproteobacteria_Moraxellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.099050169911882688,0.65610956744426974],"x":[1,2],"y":[0.013785037267618983,0.085013270910221372],"text":["Sample_type: Port uptake<br />Mean: 0.01378504<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Sample_type: BWT<br />Mean: 0.08501327<br />Class_family: Gammaproteobacteria_Oxalobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(195,61,127,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Oxalobacteraceae","legendgroup":"Gammaproteobacteria_Oxalobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036],"base":[0.71542956841954009,0.46990897041686136],"x":[3,4],"y":[0.12726045250886509,0.42845982853992642],"text":["Sample_type: Ocean uptake<br />Mean: 0.12726045<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","Sample_type: BWT+BWE<br />Mean: 0.42845983<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,72,118,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudoalteromonadaceae","legendgroup":"Gammaproteobacteria_Pseudoalteromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[0.03156133722687253,0.016106571071065758,0.10698613425597461],"x":[1,2,4],"y":[0.067488832685010158,0.64000299637320401,0.36292283616088672],"text":["Sample_type: Port uptake<br />Mean: 0.06748883<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample_type: BWT<br />Mean: 0.64000300<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Sample_type: BWT+BWE<br />Mean: 0.36292284<br />Class_family: Gammaproteobacteria_Pseudomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(213,83,110,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudomonadaceae","legendgroup":"Gammaproteobacteria_Pseudomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036],"base":[0.49698650050018384,0.035532549168628705],"x":[3,4],"y":[0.21844306791935625,0.071453585087345908],"text":["Sample_type: Ocean uptake<br />Mean: 0.21844307<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","Sample_type: BWT+BWE<br />Mean: 0.07145359<br />Class_family: Gammaproteobacteria_Saccharospirillaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(220,93,102,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Saccharospirillaceae","legendgroup":"Gammaproteobacteria_Saccharospirillaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[0,0.47032518614457058,0],"x":[2,3,4],"y":[0.016106571071065758,0.026661314355613253,0.035532549168628705],"text":["Sample_type: BWT<br />Mean: 0.01610657<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample_type: Ocean uptake<br />Mean: 0.02666131<br />Class_family: Gammaproteobacteria_Shewanellaceae","Sample_type: BWT+BWE<br />Mean: 0.03553255<br />Class_family: Gammaproteobacteria_Shewanellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(248,130,92,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Shewanellaceae","legendgroup":"Gammaproteobacteria_Shewanellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000036,"base":0,"x":[3],"y":[0.47032518614457058],"text":"Sample_type: Ocean uptake<br />Mean: 0.47032519<br />Class_family: Gammaproteobacteria_Vibrionaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(245,233,38,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Vibrionaceae","legendgroup":"Gammaproteobacteria_Vibrionaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.017980068663630357,"x":[1],"y":[0.013581268563242172],"text":"Sample_type: Port uptake<br />Mean: 0.01358127<br />Class_family: Planctomycetes_Isosphaeraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(191,191,191,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Planctomycetes_Isosphaeraceae","legendgroup":"Planctomycetes_Isosphaeraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0,"x":[1],"y":[0.017980068663630357],"text":"Sample_type: Port uptake<br />Mean: 0.01798007<br />Class_family: Verrucomicrobiae_uncultured","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_uncultured","legendgroup":"Verrucomicrobiae_uncultured","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":69.345406856584162,"l":53.466168534661698},"paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,4.5999999999999996],"tickmode":"array","ticktext":["Port uptake","BWT","Ocean uptake","BWT+BWE"],"tickvals":[1,2,3,4],"categoryorder":"array","categoryarray":["Port uptake","BWT","Ocean uptake","BWT+BWE"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"Sample Type","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.13283520132835203,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835205},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"Relative Abundance (Family > 1%) <br />","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"title":{"text":"Class_family","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"1d381316114a":{"x":{},"y":{},"fill":{},"type":"bar"}},"cur_data":"1d381316114a","visdat":{"1d381316114a":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

## 7.4 Family level, grouped by sample type and voyage

Creating the same bar plot as above, but with subsets within sample type
showing the different voyages, since similar sample types that were
collected during different voyages cluster separately on the NMDS and
PCoA plots.

``` r
BallastSeqR_family_voyage <- Ballast_physeq %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  group_by(Sample_type, Voyage, Kingdom, Phylum, Class, Order, Family) %>%
  dplyr::summarize(Mean = mean(Abundance, na.rm = TRUE), .groups = "keep") %>%
  filter(Mean > 0.01) %>%
  arrange(Class)

BallastSeqR_family_voyage$Class_family <- paste(BallastSeqR_family_voyage$Class,
                                                BallastSeqR_family_voyage$Family, sep = "_")

n_family_voyage_taxa <- length(unique(BallastSeqR_family_voyage$Class_family))
family_voyage_colors <- colorRampPalette(c(bcc_hex))(n_family_voyage_taxa)

BallastSeqR_voyage_family_bar <- ggplot(BallastSeqR_family_voyage,
                                        aes(x = Voyage, y = Mean, fill = Class_family)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, position = "fill") +
  scale_fill_manual(values = family_voyage_colors) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 1%) \n") + xlab("Voyage") +
  theme_minimal() + theme_ballast() +
  facet_grid(. ~ Sample_type, margins = TRUE, scale = "free")

print(BallastSeqR_voyage_family_bar)
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20family%20sample%20type%20and%20voyage-1.png)<!-- -->

``` r
BallastSeqR_voyage_family_bar_interactive <- ggplotly(BallastSeqR_voyage_family_bar)
BallastSeqR_voyage_family_bar_interactive
```

<div id="htmlwidget-7d85057e286831d57b36"
class="plotly html-widget html-fill-item"
style="width:672px;height:480px;">

</div>

<script type="application/json" data-for="htmlwidget-7d85057e286831d57b36">{"x":{"data":[{"orientation":"v","width":0.90000000000000013,"base":0.91291991277266515,"x":[2],"y":[0.087080087227334846],"text":"Voyage: Voyage 2<br />Mean: 0.087080087<br />Class_family: Actinobacteria_Sporichthyaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,0,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Sporichthyaceae","legendgroup":"Actinobacteria_Sporichthyaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.95852328094837203,"x":[2],"y":[0.041476719051627975],"text":"Voyage: Voyage 2<br />Mean: 0.041476719<br />Class_family: Actinobacteria_Sporichthyaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,0,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Sporichthyaceae","legendgroup":"Actinobacteria_Sporichthyaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.98393148707175926,"x":[1],"y":[0.016068512928240741],"text":"Voyage: Voyage 1<br />Mean: 0.016068513<br />Class_family: Alphaproteobacteria_Clade_I","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(69,6,89,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Clade_I","legendgroup":"Alphaproteobacteria_Clade_I","showlegend":true,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.99598889259153922,"x":[1],"y":[0.004011107408460779],"text":"Voyage: Voyage 1<br />Mean: 0.004011107<br />Class_family: Alphaproteobacteria_Clade_I","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(69,6,89,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Clade_I","legendgroup":"Alphaproteobacteria_Clade_I","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.89828689209290558,"x":[2],"y":[0.014633020679759579],"text":"Voyage: Voyage 2<br />Mean: 0.014633021<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(68,51,125,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","legendgroup":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.95155349363033037,"x":[2],"y":[0.0069697873180416536],"text":"Voyage: Voyage 2<br />Mean: 0.006969787<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(68,51,125,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","legendgroup":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.87295354890163501,"x":[2],"y":[0.025333343191270563],"text":"Voyage: Voyage 2<br />Mean: 0.025333343<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(57,85,140,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhodobacteraceae","legendgroup":"Alphaproteobacteria_Rhodobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.98146259887823217,"x":[1],"y":[0.018537401121767827],"text":"Voyage: Voyage 1<br />Mean: 0.018537401<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(57,85,140,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhodobacteraceae","legendgroup":"Alphaproteobacteria_Rhodobacteraceae","showlegend":false,"xaxis":"x4","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[0.93948708451615837,0.99122835327255354],"x":[2,1],"y":[0.012066409114172005,0.0047605393189856793],"text":["Voyage: Voyage 2<br />Mean: 0.012066409<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Voyage: Voyage 1<br />Mean: 0.004760539<br />Class_family: Alphaproteobacteria_Rhodobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(57,85,140,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhodobacteraceae","legendgroup":"Alphaproteobacteria_Rhodobacteraceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.98337705819456411,0.84971699914334664],"x":[1,2],"y":[0.016622941805435887,0.023236549758288372],"text":["Voyage: Voyage 1<br />Mean: 0.016622942<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Voyage: Voyage 2<br />Mean: 0.023236550<br />Class_family: Alphaproteobacteria_Sphingomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(53,93,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Sphingomonadaceae","legendgroup":"Alphaproteobacteria_Sphingomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.97035754985952827,"x":[1],"y":[0.011105049018703905],"text":"Voyage: Voyage 1<br />Mean: 0.011105049<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(53,93,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Sphingomonadaceae","legendgroup":"Alphaproteobacteria_Sphingomonadaceae","showlegend":false,"xaxis":"x4","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.89999999999999991],"base":[0.98435480462370195,0.92841938951560521,0.98837649598438571],"x":[1,2,1],"y":[0.0040216913606837545,0.011067695000553157,0.0028518572881678361],"text":["Voyage: Voyage 1<br />Mean: 0.004021691<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Voyage: Voyage 2<br />Mean: 0.011067695<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Voyage: Voyage 1<br />Mean: 0.002851857<br />Class_family: Alphaproteobacteria_Sphingomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(53,93,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Sphingomonadaceae","legendgroup":"Alphaproteobacteria_Sphingomonadaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.83332800128177154,"x":[2],"y":[0.016388997861575105],"text":"Voyage: Voyage 2<br />Mean: 0.016388998<br />Class_family: Alphaproteobacteria_Tistrellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(50,100,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Tistrellaceae","legendgroup":"Alphaproteobacteria_Tistrellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.92061322071759211,"x":[2],"y":[0.0078061687980131023],"text":"Voyage: Voyage 2<br />Mean: 0.007806169<br />Class_family: Alphaproteobacteria_Tistrellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(50,100,141,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Tistrellaceae","legendgroup":"Alphaproteobacteria_Tistrellaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.7712731857956906,"x":[2],"y":[0.062054815486080939],"text":"Voyage: Voyage 2<br />Mean: 0.062054815<br />Class_family: Bacteroidia_Chitinophagaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(46,108,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Chitinophagaceae","legendgroup":"Bacteroidia_Chitinophagaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.89105617467332898,"x":[2],"y":[0.029557046044263124],"text":"Voyage: Voyage 2<br />Mean: 0.029557046<br />Class_family: Bacteroidia_Chitinophagaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(46,108,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Chitinophagaceae","legendgroup":"Bacteroidia_Chitinophagaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.97228404521621026,0.72730823351906415],"x":[1,2],"y":[0.011093012978353856,0.043964952276626446],"text":["Voyage: Voyage 1<br />Mean: 0.011093013<br />Class_family: Bacteroidia_Crocinitomicaceae","Voyage: Voyage 2<br />Mean: 0.043964952<br />Class_family: Bacteroidia_Crocinitomicaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(43,114,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Crocinitomicaceae","legendgroup":"Bacteroidia_Crocinitomicaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.98167100322553991,0.87011542884087134],"x":[1,2],"y":[0.0026838013981620445,0.020940745832457641],"text":["Voyage: Voyage 1<br />Mean: 0.002683801<br />Class_family: Bacteroidia_Crocinitomicaceae","Voyage: Voyage 2<br />Mean: 0.020940746<br />Class_family: Bacteroidia_Crocinitomicaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(43,114,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Crocinitomicaceae","legendgroup":"Bacteroidia_Crocinitomicaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.95153003172503059,"x":[1],"y":[0.020754013491179668],"text":"Voyage: Voyage 1<br />Mean: 0.020754013<br />Class_family: Bacteroidia_Cryomorphaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(41,121,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cryomorphaceae","legendgroup":"Bacteroidia_Cryomorphaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.97664985608239741,"x":[1],"y":[0.0050211471431425014],"text":"Voyage: Voyage 1<br />Mean: 0.005021147<br />Class_family: Bacteroidia_Cryomorphaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(41,121,142,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cryomorphaceae","legendgroup":"Bacteroidia_Cryomorphaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.94053404691487763,"x":[1],"y":[0.010995984810152959],"text":"Voyage: Voyage 1<br />Mean: 0.010995985<br />Class_family: Bacteroidia_Cyclobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,151,139,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cyclobacteriaceae","legendgroup":"Bacteroidia_Cyclobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.97398952931046956,"x":[1],"y":[0.0026603267719278456],"text":"Voyage: Voyage 1<br />Mean: 0.002660327<br />Class_family: Bacteroidia_Cyclobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,151,139,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cyclobacteriaceae","legendgroup":"Bacteroidia_Cyclobacteriaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.70120645305388218,"x":[2],"y":[0.026101780465181967],"text":"Voyage: Voyage 2<br />Mean: 0.026101780<br />Class_family: Bacteroidia_env.OPS_17","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(30,158,136,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_env.OPS_17","legendgroup":"Bacteroidia_env.OPS_17","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.85768300887469351,"x":[2],"y":[0.012432419966177832],"text":"Voyage: Voyage 2<br />Mean: 0.012432420<br />Class_family: Bacteroidia_env.OPS_17","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(30,158,136,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_env.OPS_17","legendgroup":"Bacteroidia_env.OPS_17","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.78477861630860202,0.63702129719421319],"x":[1,2],"y":[0.15575543060627561,0.064185155859668996],"text":["Voyage: Voyage 1<br />Mean: 0.155755431<br />Class_family: Bacteroidia_Flavobacteriaceae","Voyage: Voyage 2<br />Mean: 0.064185156<br />Class_family: Bacteroidia_Flavobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(33,165,133,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Flavobacteriaceae","legendgroup":"Bacteroidia_Flavobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.82564185546932345,"x":[2],"y":[0.17435814453067655],"text":"Voyage: Voyage 2<br />Mean: 0.174358145<br />Class_family: Bacteroidia_Flavobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(33,165,133,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Flavobacteriaceae","legendgroup":"Bacteroidia_Flavobacteriaceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.96769512599305452,"x":[1],"y":[0.016236361078704742],"text":"Voyage: Voyage 1<br />Mean: 0.016236361<br />Class_family: Bacteroidia_Flavobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(33,165,133,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Flavobacteriaceae","legendgroup":"Bacteroidia_Flavobacteriaceae","showlegend":false,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000013,0.89999999999999991],"base":[0.9322536458638877,0.73580085792284677,0.7663725966709537,0.96993652275659115],"x":[1,2,2,1],"y":[0.037682876892703443,0.030571738748106925,0.09131041220373981,0.0040530065538784132],"text":["Voyage: Voyage 1<br />Mean: 0.037682877<br />Class_family: Bacteroidia_Flavobacteriaceae","Voyage: Voyage 2<br />Mean: 0.030571739<br />Class_family: Bacteroidia_Flavobacteriaceae","Voyage: Voyage 2<br />Mean: 0.091310412<br />Class_family: Bacteroidia_Flavobacteriaceae","Voyage: Voyage 1<br />Mean: 0.004053007<br />Class_family: Bacteroidia_Flavobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(33,165,133,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Flavobacteriaceae","legendgroup":"Bacteroidia_Flavobacteriaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.58333650670381543,"x":[2],"y":[0.05368479049039776],"text":"Voyage: Voyage 2<br />Mean: 0.053684790<br />Class_family: Bacteroidia_NS11-12_marine_group","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(37,172,129,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_NS11-12_marine_group","legendgroup":"Bacteroidia_NS11-12_marine_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.71023050029673762,"x":[2],"y":[0.025570357626109153],"text":"Voyage: Voyage 2<br />Mean: 0.025570358<br />Class_family: Bacteroidia_NS11-12_marine_group","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(37,172,129,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_NS11-12_marine_group","legendgroup":"Bacteroidia_NS11-12_marine_group","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.56094524210123753,"x":[2],"y":[0.022391264602577898],"text":"Voyage: Voyage 2<br />Mean: 0.022391265<br />Class_family: Bacteroidia_Sphingobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(45,178,124,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Sphingobacteriaceae","legendgroup":"Bacteroidia_Sphingobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.69956541921215809,"x":[2],"y":[0.010665081084579531],"text":"Voyage: Voyage 2<br />Mean: 0.010665081<br />Class_family: Bacteroidia_Sphingobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(45,178,124,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Sphingobacteriaceae","legendgroup":"Bacteroidia_Sphingobacteriaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.45927155236432982,"x":[2],"y":[0.10167368973690771],"text":"Voyage: Voyage 2<br />Mean: 0.101673690<br />Class_family: Bacteroidia_Spirosomaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(56,185,118,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Spirosomaceae","legendgroup":"Bacteroidia_Spirosomaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.80013479127981502,"x":[2],"y":[0.025507064189508433],"text":"Voyage: Voyage 2<br />Mean: 0.025507064<br />Class_family: Bacteroidia_Spirosomaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(56,185,118,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Spirosomaceae","legendgroup":"Bacteroidia_Spirosomaceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000013],"base":[0.63777977711203448,0.68620750839290479],"x":[2,2],"y":[0.048427731280870301,0.013357910819253305],"text":["Voyage: Voyage 2<br />Mean: 0.048427731<br />Class_family: Bacteroidia_Spirosomaceae","Voyage: Voyage 2<br />Mean: 0.013357911<br />Class_family: Bacteroidia_Spirosomaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(56,185,118,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Spirosomaceae","legendgroup":"Bacteroidia_Spirosomaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.94803014017525467,"x":[1],"y":[0.019664985817799852],"text":"Voyage: Voyage 1<br />Mean: 0.019664986<br />Class_family: Cyanobacteriia_Cyanobiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(89,199,100,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Cyanobiaceae","legendgroup":"Cyanobacteriia_Cyanobiaceae","showlegend":true,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.92734476782064179,"x":[1],"y":[0.0049088780432459123],"text":"Voyage: Voyage 1<br />Mean: 0.004908878<br />Class_family: Cyanobacteriia_Cyanobiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(89,199,100,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Cyanobiaceae","legendgroup":"Cyanobacteriia_Cyanobiaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.21812525834657009,"x":[1],"y":[0.56665335796203187],"text":"Voyage: Voyage 1<br />Mean: 0.566653358<br />Class_family: Gammaproteobacteria_Aeromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(187,222,41,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Aeromonadaceae","legendgroup":"Gammaproteobacteria_Aeromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.79025080810569925,"x":[1],"y":[0.13709395971494254],"text":"Voyage: Voyage 1<br />Mean: 0.137093960<br />Class_family: Gammaproteobacteria_Aeromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(187,222,41,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Aeromonadaceae","legendgroup":"Gammaproteobacteria_Aeromonadaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.44316874887244789,"x":[2],"y":[0.016102803491881923],"text":"Voyage: Voyage 2<br />Mean: 0.016102803<br />Class_family: Gammaproteobacteria_Alcaligenaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(232,228,26,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alcaligenaceae","legendgroup":"Gammaproteobacteria_Alcaligenaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.6301099242485122,"x":[2],"y":[0.0076698528635222862],"text":"Voyage: Voyage 2<br />Mean: 0.007669853<br />Class_family: Gammaproteobacteria_Alcaligenaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(232,228,26,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alcaligenaceae","legendgroup":"Gammaproteobacteria_Alcaligenaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.82631436930453606,"x":[1],"y":[0.17368563069546394],"text":"Voyage: Voyage 1<br />Mean: 0.173685631<br />Class_family: Gammaproteobacteria_Alteromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(210,123,173,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alteromonadaceae","legendgroup":"Gammaproteobacteria_Alteromonadaceae","showlegend":true,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.89836879895678778,"x":[1],"y":[0.071988750902740484],"text":"Voyage: Voyage 1<br />Mean: 0.071988751<br />Class_family: Gammaproteobacteria_Alteromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(210,123,173,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alteromonadaceae","legendgroup":"Gammaproteobacteria_Alteromonadaceae","showlegend":false,"xaxis":"x4","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.89999999999999991],"base":[0.72805884471313309,0.77176357373309146],"x":[1,1],"y":[0.043704729019958366,0.018487234372607797],"text":["Voyage: Voyage 1<br />Mean: 0.043704729<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Voyage: Voyage 1<br />Mean: 0.018487234<br />Class_family: Gammaproteobacteria_Alteromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(210,123,173,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alteromonadaceae","legendgroup":"Gammaproteobacteria_Alteromonadaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.39396198452817205,"x":[2],"y":[0.049206764344275844],"text":"Voyage: Voyage 2<br />Mean: 0.049206764<br />Class_family: Gammaproteobacteria_Burkholderiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(160,52,213,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Burkholderiaceae","legendgroup":"Gammaproteobacteria_Burkholderiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.60667247482544207,"x":[2],"y":[0.023437449423070134],"text":"Voyage: Voyage 2<br />Mean: 0.023437449<br />Class_family: Gammaproteobacteria_Burkholderiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(160,52,213,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Burkholderiaceae","legendgroup":"Gammaproteobacteria_Burkholderiaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.20499435725628257,"x":[1],"y":[0.013130901090287528],"text":"Voyage: Voyage 1<br />Mean: 0.013130901<br />Class_family: Gammaproteobacteria_Chromobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(165,31,152,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Chromobacteriaceae","legendgroup":"Gammaproteobacteria_Chromobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.7248820044106824,"x":[1],"y":[0.0031768403024506897],"text":"Voyage: Voyage 1<br />Mean: 0.003176840<br />Class_family: Gammaproteobacteria_Chromobacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(165,31,152,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Chromobacteriaceae","legendgroup":"Gammaproteobacteria_Chromobacteriaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.84269002092840517,"x":[1],"y":[0.1053401192468495],"text":"Voyage: Voyage 1<br />Mean: 0.105340119<br />Class_family: Gammaproteobacteria_Colwelliaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(173,39,146,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Colwelliaceae","legendgroup":"Gammaproteobacteria_Colwelliaceae","showlegend":true,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.69858644522570168,"x":[1],"y":[0.026295559184980721],"text":"Voyage: Voyage 1<br />Mean: 0.026295559<br />Class_family: Gammaproteobacteria_Colwelliaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(173,39,146,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Colwelliaceae","legendgroup":"Gammaproteobacteria_Colwelliaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.17017189044074787,0.1405490961407827],"x":[1,2],"y":[0.034822466815534692,0.25341288838738935],"text":["Voyage: Voyage 1<br />Mean: 0.034822467<br />Class_family: Gammaproteobacteria_Comamonadaceae","Voyage: Voyage 2<br />Mean: 0.253412888<br />Class_family: Gammaproteobacteria_Comamonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(182,48,138,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Comamonadaceae","legendgroup":"Gammaproteobacteria_Comamonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.78748733738039556,0.77908162090921684],"x":[1,2],"y":[0.038827031924140498,0.021053170370598173],"text":["Voyage: Voyage 1<br />Mean: 0.038827032<br />Class_family: Gammaproteobacteria_Comamonadaceae","Voyage: Voyage 2<br />Mean: 0.021053170<br />Class_family: Gammaproteobacteria_Comamonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(182,48,138,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Comamonadaceae","legendgroup":"Gammaproteobacteria_Comamonadaceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.89999999999999991,0.90000000000000013],"base":[0.68039153612556413,0.47494510728806988,0.68881635139526154,0.59564704398073232],"x":[1,2,1,2],"y":[0.0084248152696974055,0.12070193669266244,0.0097700938304401408,0.011025430844709749],"text":["Voyage: Voyage 1<br />Mean: 0.008424815<br />Class_family: Gammaproteobacteria_Comamonadaceae","Voyage: Voyage 2<br />Mean: 0.120701937<br />Class_family: Gammaproteobacteria_Comamonadaceae","Voyage: Voyage 1<br />Mean: 0.009770094<br />Class_family: Gammaproteobacteria_Comamonadaceae","Voyage: Voyage 2<br />Mean: 0.011025431<br />Class_family: Gammaproteobacteria_Comamonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(182,48,138,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Comamonadaceae","legendgroup":"Gammaproteobacteria_Comamonadaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.10191535976421803,"x":[2],"y":[0.03863373637656467],"text":"Voyage: Voyage 2<br />Mean: 0.038633736<br />Class_family: Gammaproteobacteria_Methylophilaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(190,56,132,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylophilaceae","legendgroup":"Gammaproteobacteria_Methylophilaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.75863019118817243,"x":[2],"y":[0.020451429721044412],"text":"Voyage: Voyage 2<br />Mean: 0.020451430<br />Class_family: Gammaproteobacteria_Methylophilaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(190,56,132,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylophilaceae","legendgroup":"Gammaproteobacteria_Methylophilaceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000013],"base":[0.44583334603918989,0.46423480474951634],"x":[2,2],"y":[0.018401458710326446,0.01071030253855354],"text":["Voyage: Voyage 2<br />Mean: 0.018401459<br />Class_family: Gammaproteobacteria_Methylophilaceae","Voyage: Voyage 2<br />Mean: 0.010710303<br />Class_family: Gammaproteobacteria_Methylophilaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(190,56,132,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylophilaceae","legendgroup":"Gammaproteobacteria_Methylophilaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.13936766906258086,"x":[1],"y":[0.030804221378167018],"text":"Voyage: Voyage 1<br />Mean: 0.030804221<br />Class_family: Gammaproteobacteria_Moraxellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(197,63,125,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Moraxellaceae","legendgroup":"Gammaproteobacteria_Moraxellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.73126949471499314,"x":[1],"y":[0.056217842665402418],"text":"Voyage: Voyage 1<br />Mean: 0.056217843<br />Class_family: Gammaproteobacteria_Moraxellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(197,63,125,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Moraxellaceae","legendgroup":"Gammaproteobacteria_Moraxellaceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.89999999999999991],"base":[0.65879271494316693,0.66624537118411264],"x":[1,1],"y":[0.0074526562409457098,0.014146164941451489],"text":["Voyage: Voyage 1<br />Mean: 0.007452656<br />Class_family: Gammaproteobacteria_Moraxellaceae","Voyage: Voyage 1<br />Mean: 0.014146165<br />Class_family: Gammaproteobacteria_Moraxellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(197,63,125,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Moraxellaceae","legendgroup":"Gammaproteobacteria_Moraxellaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.077381111256526219,"x":[2],"y":[0.024534248507691814],"text":"Voyage: Voyage 2<br />Mean: 0.024534249<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,71,119,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Oxalobacteraceae","legendgroup":"Gammaproteobacteria_Oxalobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.58895515168171342,"x":[2],"y":[0.16967503950645901],"text":"Voyage: Voyage 2<br />Mean: 0.169675040<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,71,119,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Oxalobacteraceae","legendgroup":"Gammaproteobacteria_Oxalobacteraceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000013],"base":[0.34528965450061871,0.35697545047185802],"x":[2,2],"y":[0.01168579597123931,0.088857895567331868],"text":["Voyage: Voyage 2<br />Mean: 0.011685796<br />Class_family: Gammaproteobacteria_Oxalobacteraceae","Voyage: Voyage 2<br />Mean: 0.088857896<br />Class_family: Gammaproteobacteria_Oxalobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,71,119,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Oxalobacteraceae","legendgroup":"Gammaproteobacteria_Oxalobacteraceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.71542956841954009,"x":[1],"y":[0.12726045250886509],"text":"Voyage: Voyage 1<br />Mean: 0.127260453<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(210,79,112,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudoalteromonadaceae","legendgroup":"Gammaproteobacteria_Pseudoalteromonadaceae","showlegend":true,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.46990897041686136,"x":[1],"y":[0.42845982853992642],"text":"Voyage: Voyage 1<br />Mean: 0.428459829<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(210,79,112,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudoalteromonadaceae","legendgroup":"Gammaproteobacteria_Pseudoalteromonadaceae","showlegend":false,"xaxis":"x4","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.89999999999999991],"base":[0.51699368834792359,0.54876111736528133],"x":[1,1],"y":[0.03176742901735774,0.11003159757788561],"text":["Voyage: Voyage 1<br />Mean: 0.031767429<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","Voyage: Voyage 1<br />Mean: 0.110031598<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(210,79,112,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudoalteromonadaceae","legendgroup":"Gammaproteobacteria_Pseudoalteromonadaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.014716678277633564,"x":[1],"y":[0.12465099078494729],"text":"Voyage: Voyage 1<br />Mean: 0.124650991<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(216,88,106,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudomonadaceae","legendgroup":"Gammaproteobacteria_Pseudomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[0.032624731616545513,0],"x":[1,2],"y":[0.69864476309844759,0.58895515168171342],"text":["Voyage: Voyage 1<br />Mean: 0.698644763<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Voyage: Voyage 2<br />Mean: 0.588955152<br />Class_family: Gammaproteobacteria_Pseudomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(216,88,106,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudomonadaceae","legendgroup":"Gammaproteobacteria_Pseudomonadaceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.10698613425597461,"x":[1],"y":[0.36292283616088672],"text":"Voyage: Voyage 1<br />Mean: 0.362922836<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(216,88,106,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudomonadaceae","legendgroup":"Gammaproteobacteria_Pseudomonadaceae","showlegend":false,"xaxis":"x4","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.89999999999999991,0.90000000000000013,0.89999999999999991],"base":[0.21783403811280602,0.24799162505548197,0.036857044057739749,0.42379246666896203],"x":[1,1,2,1],"y":[0.030157586942675951,0.17580084161348006,0.30843261044287895,0.093201221678961554],"text":["Voyage: Voyage 1<br />Mean: 0.030157587<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Voyage: Voyage 1<br />Mean: 0.175800842<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Voyage: Voyage 2<br />Mean: 0.308432610<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Voyage: Voyage 1<br />Mean: 0.093201222<br />Class_family: Gammaproteobacteria_Pseudomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(216,88,106,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudomonadaceae","legendgroup":"Gammaproteobacteria_Pseudomonadaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.49698650050018384,"x":[1],"y":[0.21844306791935625],"text":"Voyage: Voyage 1<br />Mean: 0.218443068<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(222,96,100,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Saccharospirillaceae","legendgroup":"Gammaproteobacteria_Saccharospirillaceae","showlegend":true,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.035532549168628705,"x":[1],"y":[0.071453585087345908],"text":"Voyage: Voyage 1<br />Mean: 0.071453585<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(222,96,100,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Saccharospirillaceae","legendgroup":"Gammaproteobacteria_Saccharospirillaceae","showlegend":false,"xaxis":"x4","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.89999999999999991],"base":[0.14495532119315716,0.1994842382106306],"x":[1,1],"y":[0.054528917017473449,0.018349799902175418],"text":["Voyage: Voyage 1<br />Mean: 0.054528917<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","Voyage: Voyage 1<br />Mean: 0.018349800<br />Class_family: Gammaproteobacteria_Saccharospirillaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(222,96,100,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Saccharospirillaceae","legendgroup":"Gammaproteobacteria_Saccharospirillaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0,"x":[1],"y":[0.014716678277633564],"text":"Voyage: Voyage 1<br />Mean: 0.014716678<br />Class_family: Gammaproteobacteria_Shewanellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(249,126,99,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Shewanellaceae","legendgroup":"Gammaproteobacteria_Shewanellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0,"x":[1],"y":[0.032624731616545513],"text":"Voyage: Voyage 1<br />Mean: 0.032624732<br />Class_family: Gammaproteobacteria_Shewanellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(249,126,99,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Shewanellaceae","legendgroup":"Gammaproteobacteria_Shewanellaceae","showlegend":false,"xaxis":"x2","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0.47032518614457058,"x":[1],"y":[0.026661314355613253],"text":"Voyage: Voyage 1<br />Mean: 0.026661314<br />Class_family: Gammaproteobacteria_Shewanellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(249,126,99,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Shewanellaceae","legendgroup":"Gammaproteobacteria_Shewanellaceae","showlegend":false,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0,"x":[1],"y":[0.035532549168628705],"text":"Voyage: Voyage 1<br />Mean: 0.035532549<br />Class_family: Gammaproteobacteria_Shewanellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(249,126,99,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Shewanellaceae","legendgroup":"Gammaproteobacteria_Shewanellaceae","showlegend":false,"xaxis":"x4","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.89999999999999991,0.89999999999999991,0.89999999999999991],"base":[0.1174050670995567,0.12096556430756381,0.12917496567796299,0.13583030435234345],"x":[1,1,1,1],"y":[0.0035604972080071073,0.0082094013703991886,0.006655338674380451,0.00912501684081371],"text":["Voyage: Voyage 1<br />Mean: 0.003560497<br />Class_family: Gammaproteobacteria_Shewanellaceae","Voyage: Voyage 1<br />Mean: 0.008209401<br />Class_family: Gammaproteobacteria_Shewanellaceae","Voyage: Voyage 1<br />Mean: 0.006655339<br />Class_family: Gammaproteobacteria_Shewanellaceae","Voyage: Voyage 1<br />Mean: 0.009125017<br />Class_family: Gammaproteobacteria_Shewanellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(249,126,99,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Shewanellaceae","legendgroup":"Gammaproteobacteria_Shewanellaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0,"x":[1],"y":[0.47032518614457058],"text":"Voyage: Voyage 1<br />Mean: 0.470325186<br />Class_family: Gammaproteobacteria_Vibrionaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(247,213,42,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Vibrionaceae","legendgroup":"Gammaproteobacteria_Vibrionaceae","showlegend":true,"xaxis":"x3","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.89999999999999991,"base":0,"x":[1],"y":[0.1174050670995567],"text":"Voyage: Voyage 1<br />Mean: 0.117405067<br />Class_family: Gammaproteobacteria_Vibrionaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(247,213,42,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Vibrionaceae","legendgroup":"Gammaproteobacteria_Vibrionaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.051972418448809801,"x":[2],"y":[0.025408692807716418],"text":"Voyage: Voyage 2<br />Mean: 0.025408693<br />Class_family: Planctomycetes_Isosphaeraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(44,101,134,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Planctomycetes_Isosphaeraceae","legendgroup":"Planctomycetes_Isosphaeraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.02475474551153482,"x":[2],"y":[0.012102298546204929],"text":"Voyage: Voyage 2<br />Mean: 0.012102299<br />Class_family: Planctomycetes_Isosphaeraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(44,101,134,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Planctomycetes_Isosphaeraceae","legendgroup":"Planctomycetes_Isosphaeraceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.021710572514301505,"x":[2],"y":[0.030261845934508296],"text":"Voyage: Voyage 2<br />Mean: 0.030261846<br />Class_family: Verrucomicrobiae_uncultured","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(157,157,157,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_uncultured","legendgroup":"Verrucomicrobiae_uncultured","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0.010340863741613396,"x":[2],"y":[0.014413881769921423],"text":"Voyage: Voyage 2<br />Mean: 0.014413882<br />Class_family: Verrucomicrobiae_uncultured","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(157,157,157,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_uncultured","legendgroup":"Verrucomicrobiae_uncultured","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0,"x":[2],"y":[0.021710572514301505],"text":"Voyage: Voyage 2<br />Mean: 0.021710573<br />Class_family: Verrucomicrobiae_Verrucomicrobiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Verrucomicrobiaceae","legendgroup":"Verrucomicrobiae_Verrucomicrobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.90000000000000013,"base":0,"x":[2],"y":[0.010340863741613396],"text":"Voyage: Voyage 2<br />Mean: 0.010340864<br />Class_family: Verrucomicrobiae_Verrucomicrobiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.1338582677165354,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Verrucomicrobiaceae","legendgroup":"Verrucomicrobiae_Verrucomicrobiaceae","showlegend":false,"xaxis":"x5","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":39.246160232461605,"r":23.246160232461605,"b":55.670675597908364,"l":53.466168534661698},"paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,0.19184605348988912],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,2.6000000000000001],"tickmode":"array","ticktext":["Voyage 1","Voyage 2"],"tickvals":[1,2],"categoryorder":"array","categoryarray":["Voyage 1","Voyage 2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"annotations":[{"text":"Voyage","x":0.5,"y":0,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"top","annotationType":"axis","yshift":-38.999857831200181},{"text":"Relative Abundance (Family > 1%) <br />","x":0,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis","xshift":-36.795350767953508},{"text":"Port uptake","x":0.095923026744944559,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"BWT","x":0.29999999999999999,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"Ocean uptake","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"BWT+BWE","x":0.70000000000000007,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"(all)","x":0.90407697325505554,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.050000000000000003,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.13283520132835203,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835205},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.19184605348988912,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.19184605348988912,"y0":0,"y1":27.629721876297225,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0.2081539465101109,"x1":0.3918460534898891,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.2081539465101109,"x1":0.3918460534898891,"y0":0,"y1":27.629721876297225,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0.40815394651011094,"x1":0.59184605348988917,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.40815394651011094,"x1":0.59184605348988917,"y0":0,"y1":27.629721876297225,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0.60815394651011101,"x1":0.79184605348988912,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.60815394651011101,"x1":0.79184605348988912,"y0":0,"y1":27.629721876297225,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0.80815394651011097,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.80815394651011097,"x1":1,"y0":0,"y1":27.629721876297225,"yanchor":1,"ysizemode":"pixel"}],"xaxis2":{"type":"linear","autorange":false,"range":[0.40000000000000002,2.6000000000000001],"tickmode":"array","ticktext":["Voyage 1","Voyage 2"],"tickvals":[1,2],"categoryorder":"array","categoryarray":["Voyage 1","Voyage 2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.2081539465101109,0.3918460534898891],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"xaxis3":{"type":"linear","autorange":false,"range":[0.40000000000000002,1.6000000000000001],"tickmode":"array","ticktext":["Voyage 1"],"tickvals":[1],"categoryorder":"array","categoryarray":["Voyage 1"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.40815394651011094,0.59184605348988917],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"xaxis4":{"type":"linear","autorange":false,"range":[0.40000000000000002,1.6000000000000001],"tickmode":"array","ticktext":["Voyage 1"],"tickvals":[1],"categoryorder":"array","categoryarray":["Voyage 1"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.60815394651011101,0.79184605348988912],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"xaxis5":{"type":"linear","autorange":false,"range":[0.40000000000000002,2.6000000000000001],"tickmode":"array","ticktext":["Voyage 1","Voyage 2"],"tickvals":[1,2],"categoryorder":"array","categoryarray":["Voyage 1","Voyage 2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.80815394651011097,1],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"title":{"text":"Class_family","font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"1d387964da8":{"x":{},"y":{},"fill":{},"type":"bar"}},"cur_data":"1d387964da8","visdat":{"1d387964da8":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

## 7.5 Genus level, grouped by voyage and sample type

Using `ggnested` for class-grouped legend colors.

``` r
library(ggnested)

BallastSeqR_genus_voyage <- Ballast_physeq %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  group_by(Sample_type, Voyage, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarize(Mean = mean(Abundance, na.rm = TRUE), .groups = "keep") %>%
  filter(Mean > 0.01) %>%
  arrange(Class)

BallastSeqR_genus_voyage_filt <- data.frame(BallastSeqR_genus_voyage)

pal_genus <- c("#440154FF", "blue", "green", "yellow", "purple", "orange", "#FF6699")
names(pal_genus) <- unique(BallastSeqR_genus_voyage_filt$Class)

BallastSeqR_voyage_genus_bar_ggnested <- ggnested::ggnested(
  BallastSeqR_genus_voyage_filt,
  aes(x = Sample_type, y = Mean, main_group = Class, sub_group = Genus),
  main_palette = pal_genus) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, position = "fill") +
  facet_wrap(~ Voyage, scale = "free") +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 1%) \n") + xlab("Sample Type") +
  theme_minimal() + theme_ballast()

print(BallastSeqR_voyage_genus_bar_ggnested)
```

![](16S-sequence-analysis_files/figure-gfm/bcc%20genus%20voyage%20and%20sample%20type-1.png)<!-- -->

``` r
ggsave(filename = "R output/BallastSeqR_voyage_sample_genus_bar.svg",
       plot = BallastSeqR_voyage_genus_bar_ggnested, width = 10, height = 6, device = svg)
```

## 7.6 Species level, grouped by voyage and sample type

``` r
#Examining BCC by sample type & voyage
BallastSeqR_species_voyage <- Ballast_physeq %>%
  tax_glom(taxrank = "Species") %>%                   #Agglomerate at species level
  transform_sample_counts(function(x) {x/sum(x)}) %>% #Transform to rel. abundance
  psmelt() %>%                                        #Melt to long format
  group_by(Sample_type, Voyage, 
           Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::summarize(Mean = 
                     mean(Abundance, na.rm=TRUE)) %>% #Calculate average
  arrange(Class)
```

    ## `summarise()` has regrouped the output.
    ## ℹ Summaries were computed grouped by Sample_type, Voyage, Kingdom, Phylum,
    ##   Class, Order, Family, Genus, and Species.
    ## ℹ Output is grouped by Sample_type, Voyage, Kingdom, Phylum, Class, Order,
    ##   Family, and Genus.
    ## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
    ## ℹ Use `summarise(.by = c(Sample_type, Voyage, Kingdom, Phylum, Class, Order,
    ##   Family, Genus, Species))` for per-operation grouping (`?dplyr::dplyr_by`)
    ##   instead.

``` r
#Convert to a dataframe
BallastSeqR_species_voyage_filt <- data.frame(BallastSeqR_species_voyage) 

#Custom color palette
pal <- c("#440154FF", "blue", "green", "yellow", "purple", "orange", "#FF6699",
         "#440154FF", "blue", "green", "yellow", "purple", "orange", "#FF6699",
         "#440154FF", "blue", "green", "yellow", "purple", "orange", "#FF6699",
         "#440154FF", "blue", "green", "yellow", "purple", "orange", "#FF6699",
          "#440154FF", "blue", "green", "yellow")
names(pal) <- unique(BallastSeqR_species_voyage_filt$Class)

#Use ggnested to create the bar plot
BallastSeqR_voyage_species_bar_ggnested <- ggnested::ggnested(BallastSeqR_species_voyage_filt, 
         aes(x = Sample_type, 
             y = Mean, 
             main_group = Class, 
             sub_group = Species),
          main_palette = pal) + 
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, position = "fill") +
  facet_wrap(~ Voyage, scale="free") +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Species > 1%) \n") + xlab("Sample Type") +
  theme_minimal()+theme(
    strip.text = element_text(size=12),
    panel.border = element_rect(fill=NA, colour = "black"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), 
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, 
                              colour="black"), 
    axis.text.y = element_text(size=10, colour="black"), 
    plot.title = element_text(hjust = 0.5), 
    axis.ticks.x = element_line(colour="#000000", linewidth=0.1), 
    axis.ticks.y = element_line(colour="#000000", linewidth=0.1),
    legend.text = element_text(size = 10)) +
  facet_wrap(~ Voyage, scale="free")

print(BallastSeqR_voyage_species_bar_ggnested)
```

![](16S-sequence-analysis_files/figure-gfm/ggnested%20species-1.png)<!-- -->

### 7.6.1 Create Figure 3: Faceted grouped bar plot

Top 10 species from the top 3 classes, shown as a faceted grouped bar
plot with a common baseline. Each panel shows one species, with sample
types on the x-axis and mean relative abundance (+/- SE) on the y-axis.
This layout allows direct comparison of each taxon’s abundance across
treatments.

``` r
library("fantaxtic")

# Merge samples for identifying top taxa (used only for ranking)
Ballast_physeq_merged <- speedyseq::merge_samples2(Ballast_physeq, "Type_voyage",
                                                     fun_otu = mean)

top_nested <- nested_top_taxa(Ballast_physeq_merged,
                              top_tax_level = "Class",
                              nested_tax_level = "Species",
                              n_top_taxa = 3,
                              n_nested_taxa = 10)

# Extract class-species pairs from fantaxtic (preserves which species belongs to which class)
top_class_species <- tax_table(top_nested$ps_obj) %>%
  as.data.frame() %>%
  filter(Species != "Other", Class != "Other") %>%
  distinct(Class, Species)

top_class_names <- unique(top_class_species$Class)
top_species_names <- unique(top_class_species$Species)

# Clean species names (used for matching before gsub is applied to fig3_data)
top_species_names_clean <- gsub("_", " ", top_species_names)

# Now compute per-sample relative abundances at the species level (not merged)
fig3_data <- Ballast_physeq %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x / sum(x)}) %>%
  psmelt() %>%
  # Join on class-species pairs to avoid cross-class contamination
  inner_join(top_class_species, by = c("Class", "Species")) %>%
  group_by(Sample_type, Voyage, Class, Species) %>%
  dplyr::summarize(
    Mean = mean(Abundance, na.rm = TRUE),
    SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Order Class as a factor so panels group by class
fig3_data$Class <- factor(fig3_data$Class, levels = top_class_names)

# Replace underscores with spaces in species names
fig3_data$Species <- gsub("_", " ", fig3_data$Species)

# Remove sample types that don't exist in Voyage 2
fig3_data <- fig3_data %>%
  filter(!(Voyage == "Voyage 2" & Sample_type %in% c("Ocean uptake", "BWT+BWE"))) %>%
  droplevels()
```

#### 7.6.1.1 Calculate proportion of community represented by top taxa

``` r
# Per-sample relative abundances at species level
fig3_all_species <- Ballast_physeq %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x / sum(x)}) %>%
  psmelt()

# Clean species names to match
fig3_all_species$Species <- gsub("_", " ", fig3_all_species$Species)

# Total relative abundance accounted for by the top species per sample
fig3_top_proportion <- fig3_all_species %>%
  mutate(is_top = Species %in% gsub("_", " ", top_species_names) &
                  Class %in% top_class_names) %>%
  group_by(Sample, Sample_type, Voyage) %>%
  dplyr::summarize(
    top_taxa_abundance = sum(Abundance[is_top]),
    .groups = "drop"
  )

# Summary across all samples
cat("Proportion of community composition represented by top 10 species in top 3 classes:\n")
```

    ## Proportion of community composition represented by top 10 species in top 3 classes:

``` r
cat(sprintf("  Mean:   %.1f%%\n", mean(fig3_top_proportion$top_taxa_abundance) * 100))
```

    ##   Mean:   76.9%

``` r
cat(sprintf("  Median: %.1f%%\n", median(fig3_top_proportion$top_taxa_abundance) * 100))
```

    ##   Median: 80.0%

``` r
cat(sprintf("  Range:  %.1f%% - %.1f%%\n",
            min(fig3_top_proportion$top_taxa_abundance) * 100,
            max(fig3_top_proportion$top_taxa_abundance) * 100))
```

    ##   Range:  36.3% - 98.1%

``` r
# Breakdown by sample type
fig3_proportion_by_type <- fig3_top_proportion %>%
  group_by(Sample_type) %>%
  dplyr::summarize(
    Mean_pct = mean(top_taxa_abundance) * 100,
    SD_pct = sd(top_taxa_abundance) * 100,
    .groups = "drop"
  )

knitr::kable(fig3_proportion_by_type,
             col.names = c("Sample type", "Mean %", "SD %"),
             digits = 1) %>%
  kableExtra::kable_styling("striped", latex_options = "scale_down")
```

<table class="table table-striped" style="color: black; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Sample type
</th>

<th style="text-align:right;">

Mean %
</th>

<th style="text-align:right;">

SD %
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Port uptake
</td>

<td style="text-align:right;">

63.6
</td>

<td style="text-align:right;">

24.9
</td>

</tr>

<tr>

<td style="text-align:left;">

BWT
</td>

<td style="text-align:right;">

84.4
</td>

<td style="text-align:right;">

8.5
</td>

</tr>

<tr>

<td style="text-align:left;">

Ocean uptake
</td>

<td style="text-align:right;">

79.8
</td>

<td style="text-align:right;">

0.8
</td>

</tr>

<tr>

<td style="text-align:left;">

BWT+BWE
</td>

<td style="text-align:right;">

96.9
</td>

<td style="text-align:right;">

0.9
</td>

</tr>

</tbody>

</table>

``` r
# Create class-based color palette
n_classes <- length(top_class_names)
class_base_colors <- RColorBrewer::brewer.pal(max(3, n_classes), "Dark2")[1:n_classes]
names(class_base_colors) <- top_class_names

# Map each species to a shade of its parent class color
species_colors <- fig3_data %>%
  distinct(Class, Species) %>%
  arrange(Class, Species) %>%
  group_by(Class) %>%
  mutate(
    n = n(),
    idx = row_number(),
    color = colorRampPalette(
      c(class_base_colors[as.character(Class[1])],
        colorspace::lighten(class_base_colors[as.character(Class[1])], 0.6))
    )(n[1])[idx]
  ) %>%
  ungroup()

species_pal <- setNames(species_colors$color, species_colors$Species)

# Create a combined Class_Species factor for faceting, ordered by Class then Species
fig3_data <- fig3_data %>%
  arrange(Class, Species) %>%
  mutate(Class_Species = factor(Species, levels = unique(Species)))
```

``` r
fig3 <- ggplot(fig3_data, aes(x = Sample_type, y = Mean, fill = Sample_type)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(Mean - SE, 0), ymax = Mean + SE),
                width = 0.25, linewidth = 0.3) +
  facet_nested(Class + Class_Species ~ Voyage, scales = "free_y",
               nest_line = element_line(colour = "black")) +
  scale_fill_manual(values = sample_type_colors, "Sample Type") +
  labs(x = "Sample type",
       y = "Mean relative abundance",
       caption = "Relative abundance calculated as proportion of total reads per sample.\nBars show mean across biological replicates; error bars show ± 1 SE.") +
  theme_minimal() +
  theme_ballast() +
  theme(
    strip.text.y = element_text(size = 9, angle = 0),
    strip.text.x = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.caption = element_text(size = 8, hjust = 0, colour = "grey40"),
    ggh4x.facet.nestline = element_line(linetype = "solid")
  )

print(fig3)
```

![](16S-sequence-analysis_files/figure-gfm/fig3%20grouped%20bar%20plot-1.png)<!-- -->

``` r
ggsave(filename = "R output/Fig3_grouped_bar_species.svg",
       plot = fig3, width = 12, height = 22, device = svg)

ggsave(filename = "R output/Fig3_grouped_bar_species.png",
       plot = fig3, width = 12, height = 22, dpi = 300)
```

##### 7.6.1.1.1 Create this plot with patterns overlayed

``` r
#Load ggpattern and magick packages so that patterns can be added over colors
library(ggpattern)
library(magick)
```

    ## Warning: package 'magick' was built under R version 4.5.3

    ## Linking to ImageMagick 6.9.13.29
    ## Enabled features: cairo, freetype, fftw, ghostscript, heic, lcms, pango, raw, rsvg, webp
    ## Disabled features: fontconfig, x11

``` r
#This code prints a figure displaying possible magick patterns and their names that
#can be used in the following manual pattern scale function.
#From https://coolbutuseless.github.io/package/ggpattern/articles/pattern-magick.html

df1 <- data.frame(
  x    = rep(1:6, 9),
  y    = rep(1:9, each=6),
  name = gridpattern::names_magick,
  stringsAsFactors = FALSE
)


ggplot(df1) + 
  geom_tile_pattern(
    aes(x, y, pattern_type = I(name)),
    pattern       = 'magick',
    pattern_scale = 1.5,
    pattern_fill  = 'black', 
    width         = 0.9, 
    height        = 0.9
  ) + 
  geom_label(aes(x+0.4, y+0.4, label = name), hjust = 1, vjust = 1) + 
  theme_void() + 
  labs(
    title = "All the possible magick pattern names"
  ) +
  coord_fixed(1)
```

![](16S-sequence-analysis_files/figure-gfm/patterns%20and%20fantaxtic-1.png)<!-- -->

Create manual pattern scale function:

``` r
#This assigns a pattern to each taxa and is used in place of 
#scale_pattern_type_manual; make sure all taxa from all data sets are included here
#so that this function can be used to create all graphs. The value 'gray100' leaves
#an assigned taxa blank.

scale_pattern_type_bacteria <- function(...){
  ggpattern:::manual_scale(
    'pattern_type',
    values = setNames(c('gray100', 'vertical', 
                        'bricks', 'gray100',
                        'gray100', 'horizontalsaw',
                        'hs_bdiagonal', 'verticalsaw',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 
                        'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'hs_horizontal',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100',
                        'gray100', 'gray100'),
                      
                      c("Oleispira Genus", "Pseudomonas Genus",
                        "Pseudomonas_anguilliseptica", "Aeromonas Genus", 
                        "Shewanella Genus", "Pseudoalteromonas Genus", 
                        "Vibrionaceae Family", "Rheinheimera Genus", 
                        "Comamonadaceae Family", "Polynucleobacter Genus",
                        "Clade_Ia Genus", "Clade_II Genus", 
                        "Caulobacter Genus", 
                        "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus",
                        "Sphingorhabdus Genus", "Sphingomonadaceae Family", 
                        "Sulfitobacter Genus", "Roseovarius Genus", 
                        "Rhodobacteraceae Family", "Tistrella Genus", 
                        "uncultured_bacterium", "Mesoflavibacter_sp.", 
                        "Antarctic_bacterium", "Flavobacterium Genus", 
                        "Flavobacterium_jumunjinense", "Flavobacteriaceae Family",
                        "Fluviicola Genus", "Arcicella Genus", 
                        "Flectobacillus Genus", "Pseudarcicella Genus",
                        "Other", "Other Alphaproteobacteria",
                        "Other Bacteroidia", "Other Gammaproteobacteria")),
    
    ...
  )
}


#'verticalsaw', 'crosshatch', 'verticalrightshingle', 
 #                       'horizontal3', 'hs_diagcross', 'horizontalsaw',
  #                      'rightshingle', 'hexagons', 'circles', 'verticalbricks',
   #                     'smallfishscales', 'left45', 'right45', 
    #                    'verticalsaw',   'gray100',
     #                   'right30', 'hs_fdiagonal', 'hs_cross', 'vertical2'
```

Add patterns to plot:

``` r
type_voyage_order <- c("port_uptake_1", "BWT_1", "ocean_uptake_1",
                       "BWT_BWE_1", "port_uptake_2", "BWT_2")

#Plot
plot_pattern <- plot_nested_bar(ps_obj = top_nested$ps_obj,
                top_level = "Class",
                nested_level = "Species",
                sample_order = type_voyage_order) +
  facet_wrap(~ Voyage, scale="free") +
  labs(x = "Sample type", y = "Relative abundance") +
  scale_x_discrete(labels=c("Port uptake V1", "BWT V1", "Ocean uptake V1",
                          "BWT+BWE V1", "Port uptake V2", "BWT V2")) +
  geom_bar_pattern(position = "fill",
                   stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_spacing = 0.05,
                   pattern = 'magick',
                   aes(pattern_type = Species)) +
  scale_pattern_type_bacteria(guide = "none") +
  
  theme_minimal()+theme(
    strip.text = element_text(size=12),
    panel.border = element_rect(fill=NA, colour = "black"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), 
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, 
                              colour="black"), 
    axis.text.y = element_text(size=10, colour="black"), 
    plot.title = element_text(hjust = 0.5), 
    axis.ticks.x = element_line(colour="#000000", linewidth=0.1), 
    axis.ticks.y = element_line(colour="#000000", linewidth=0.1),
    legend.text = element_text(size = 10))

print(plot_pattern)
```

![](16S-sequence-analysis_files/figure-gfm/add%20patterns%20to%20fantaxtic%20plot-1.png)<!-- -->

Save as .svg:

``` r
#Saving as an editable pdf
ggsave(filename="R output/BallastSeqR_fantaxtic_patterned.pdf", plot=plot_pattern, width=10, height=6, device=cairo_pdf)
```

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=9'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=8'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=5'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=2'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=3'
    ## are unlikely values in pixels
    ## Warning in png(png_file, width = width, height = height): 'width=18, height=3'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=4'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=3'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=2'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels
    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels
    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels
    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels
    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels
    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=5'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=10'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=3'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=8'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=5'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=11'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=1'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=3'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=2'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=13'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=7'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=11'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=4'
    ## are unlikely values in pixels

    ## Warning in png(png_file, width = width, height = height): 'width=18, height=3'
    ## are unlikely values in pixels

#### 7.6.1.2 Read abundance

Same plot as above, but here I’m using fantaxtic to create a plot
displaying the number of reads for each of the top 10 species in the top
3 classes, rather than the relative abundance, as is shown above.

``` r
#Order samples
order <- c("port_uptake_1", "BWT_1", "ocean_uptake_1",
                          "BWT_BWE_1", "port_uptake_2", "BWT_2")
order <- as.character(order)

#Get the top taxa
top_level <- "Class"
nested_level <- "Species"
sample_order <- order
top_asv <- nested_top_taxa(Ballast_physeq_merged,
                              top_tax_level = "Class",
                              nested_tax_level = "Species",
                              n_top_taxa = 3, 
                              n_nested_taxa = 10)


#Create names for NA taxa
ps_tmp <- top_asv$ps_obj %>%
  name_na_taxa()

#Add labels to taxa with the same names
ps_tmp <- ps_tmp %>%
  label_duplicate_taxa(tax_level = nested_level)

#Generate a palette basedon the phyloseq object
pal <- taxon_colours(ps_tmp,
                     tax_level = top_level)

#Convert physeq to df
psdf <- psmelt(ps_tmp)


#Move the merged labels to the appropriate positions in the plot:
#Top merged labels need to be at the top of the plot,
#nested merged labels at the bottom of each group
psdf <- move_label(psdf = psdf,
                   col_name = top_level,
                   label = "Other",
                   pos = 0)
psdf <- move_nested_labels(psdf,
                           top_level = top_level,
                           nested_level = nested_level,
                           top_merged_label = "Other",
                           nested_label = "Other",
                           pos = Inf)


# Reorder samples
if(!is.null(sample_order)){
  if(all(sample_order %in% unique(psdf$Sample))){
    psdf <- psdf %>%
      mutate(Sample = factor(Sample, levels = sample_order))
  } else {
    stop("Error: not all(sample_order %in% sample_names(ps_obj)).")
  }

}


#Generate a base plot
library(ggnested)
p <- ggnested(psdf,
              aes_string(main_group = top_level,
                         sub_group = nested_level,
                         x = "Sample",
                         y = "Abundance"),
              main_palette = pal) +
  facet_wrap(~ Voyage, scale="free") +
   scale_y_continuous(limits = c(0, 11000)) +
  labs(x = "Sample type", y = "Abundance") +
  scale_x_discrete(labels=c("Port uptake V1", "BWT V1", "Ocean uptake V1",
                          "BWT+BWE V1", "Port uptake V2", "BWT V2")) +
  theme_minimal()+theme(
    strip.text = element_text(size=12),
    panel.border = element_rect(fill=NA, colour = "black"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), 
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, 
                              colour="black"), 
    axis.text.y = element_text(size=10, colour="black"), 
    plot.title = element_text(hjust = 0.5), 
    axis.ticks.x = element_line(colour="#000000", linewidth=0.1), 
    axis.ticks.y = element_line(colour="#000000", linewidth=0.1),
    legend.text = element_text(size = 10))


#Add abundances
p <- p + geom_col()

#Use geom_col(position = position_fill()) if you want to plot the relative abundance

print(p)
```

![](16S-sequence-analysis_files/figure-gfm/total%20abundance%20fantaxtic%20class%20species-1.png)<!-- -->

``` r
#Make interactive plot
p_abundance <- ggplotly(p)

p_abundance
```

<div id="htmlwidget-f2a574f3146a1df57a01"
class="plotly html-widget html-fill-item"
style="width:672px;height:480px;">

</div>

<script type="application/json" data-for="htmlwidget-f2a574f3146a1df57a01">{"x":{"data":[{"orientation":"v","width":[0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[9771.6666666666661,4443,6324.833333333333,3511],"x":[3,1,2,4],"y":[274.66666666666788,113.33333333333303,5.5,0],"text":["Class: Other<br />Species: Other<br />Sample: ocean_uptake_1<br />Abundance:  274.6666667<br />group_subgroup: Other - Other<br />group_subgroup: Other - Other","Class: Other<br />Species: Other<br />Sample: port_uptake_1<br />Abundance:  113.3333333<br />group_subgroup: Other - Other<br />group_subgroup: Other - Other","Class: Other<br />Species: Other<br />Sample: BWT_1<br />Abundance:    5.5000000<br />group_subgroup: Other - Other<br />group_subgroup: Other - Other","Class: Other<br />Species: Other<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Other - Other<br />group_subgroup: Other - Other"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(229,229,229,1)","line":{"width":1.8897637795275593,"color":"rgba(229,229,229,1)"}},"name":"Other - Other","legendgroup":"Other - Other","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1795.6666666666667,5734],"x":[1,2],"y":[426.66666666666674,19.66666666666697],"text":["Class: Other<br />Species: Other<br />Sample: port_uptake_2<br />Abundance:  426.6666667<br />group_subgroup: Other - Other<br />group_subgroup: Other - Other","Class: Other<br />Species: Other<br />Sample: BWT_2<br />Abundance:   19.6666667<br />group_subgroup: Other - Other<br />group_subgroup: Other - Other"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(229,229,229,1)","line":{"width":1.8897637795275593,"color":"rgba(229,229,229,1)"}},"name":"Other - Other","legendgroup":"Other - Other","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991,0.90000000000000036,0.90000000000000036],"base":[6307,4438.1666666666661,9771.6666666666661,3511],"x":[2,1,3,4],"y":[17.83333333333303,4.8333333333339397,0,0],"text":["Class: Alphaproteobacteria<br />Species: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />Sample: BWT_1<br />Abundance:   17.8333333<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","Class: Alphaproteobacteria<br />Species: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />Sample: port_uptake_1<br />Abundance:    4.8333333<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","Class: Alphaproteobacteria<br />Species: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","Class: Alphaproteobacteria<br />Species: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,38,65,1)","line":{"width":1.8897637795275593,"color":"rgba(0,38,65,1)"}},"name":"Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","legendgroup":"Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1795.6666666666667,5734],"x":[1,2],"y":[0,0],"text":["Class: Alphaproteobacteria<br />Species: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","Class: Alphaproteobacteria<br />Species: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus<br />group_subgroup: Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,38,65,1)","line":{"width":1.8897637795275593,"color":"rgba(0,38,65,1)"}},"name":"Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","legendgroup":"Alphaproteobacteria - Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991,0.90000000000000036,0.90000000000000036],"base":[6297,4438.1666666666661,9771.6666666666661,3511],"x":[2,1,3,4],"y":[10,0,0,0],"text":["Class: Alphaproteobacteria<br />Species: Caulobacter Genus<br />Sample: BWT_1<br />Abundance:   10.0000000<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus","Class: Alphaproteobacteria<br />Species: Caulobacter Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus","Class: Alphaproteobacteria<br />Species: Caulobacter Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus","Class: Alphaproteobacteria<br />Species: Caulobacter Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(20,57,84,1)","line":{"width":1.8897637795275593,"color":"rgba(20,57,84,1)"}},"name":"Alphaproteobacteria - Caulobacter Genus","legendgroup":"Alphaproteobacteria - Caulobacter Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1775,5734],"x":[1,2],"y":[20.666666666666742,0],"text":["Class: Alphaproteobacteria<br />Species: Caulobacter Genus<br />Sample: port_uptake_2<br />Abundance:   20.6666667<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus","Class: Alphaproteobacteria<br />Species: Caulobacter Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus<br />group_subgroup: Alphaproteobacteria - Caulobacter Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(20,57,84,1)","line":{"width":1.8897637795275593,"color":"rgba(20,57,84,1)"}},"name":"Alphaproteobacteria - Caulobacter Genus","legendgroup":"Alphaproteobacteria - Caulobacter Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[9614.6666666666661,4438.1666666666661,6297,3511],"x":[3,1,2,4],"y":[157,0,0,0],"text":["Class: Alphaproteobacteria<br />Species: Clade_Ia Genus<br />Sample: ocean_uptake_1<br />Abundance:  157.0000000<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus","Class: Alphaproteobacteria<br />Species: Clade_Ia Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus","Class: Alphaproteobacteria<br />Species: Clade_Ia Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus","Class: Alphaproteobacteria<br />Species: Clade_Ia Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(41,77,103,1)","line":{"width":1.8897637795275593,"color":"rgba(41,77,103,1)"}},"name":"Alphaproteobacteria - Clade_Ia Genus","legendgroup":"Alphaproteobacteria - Clade_Ia Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1775,5734],"x":[1,2],"y":[0,0],"text":["Class: Alphaproteobacteria<br />Species: Clade_Ia Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus","Class: Alphaproteobacteria<br />Species: Clade_Ia Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus<br />group_subgroup: Alphaproteobacteria - Clade_Ia Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(41,77,103,1)","line":{"width":1.8897637795275593,"color":"rgba(41,77,103,1)"}},"name":"Alphaproteobacteria - Clade_Ia Genus","legendgroup":"Alphaproteobacteria - Clade_Ia Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[9553,4431.5,6297,3511],"x":[3,1,2,4],"y":[61.66666666666606,6.6666666666660603,0,0],"text":["Class: Alphaproteobacteria<br />Species: Clade_II Genus<br />Sample: ocean_uptake_1<br />Abundance:   61.6666667<br />group_subgroup: Alphaproteobacteria - Clade_II Genus<br />group_subgroup: Alphaproteobacteria - Clade_II Genus","Class: Alphaproteobacteria<br />Species: Clade_II Genus<br />Sample: port_uptake_1<br />Abundance:    6.6666667<br />group_subgroup: Alphaproteobacteria - Clade_II Genus<br />group_subgroup: Alphaproteobacteria - Clade_II Genus","Class: Alphaproteobacteria<br />Species: Clade_II Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_II Genus<br />group_subgroup: Alphaproteobacteria - Clade_II Genus","Class: Alphaproteobacteria<br />Species: Clade_II Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_II Genus<br />group_subgroup: Alphaproteobacteria - Clade_II Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(61,96,122,1)","line":{"width":1.8897637795275593,"color":"rgba(61,96,122,1)"}},"name":"Alphaproteobacteria - Clade_II Genus","legendgroup":"Alphaproteobacteria - Clade_II Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1775,5734],"x":[1,2],"y":[0,0],"text":["Class: Alphaproteobacteria<br />Species: Clade_II Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_II Genus<br />group_subgroup: Alphaproteobacteria - Clade_II Genus","Class: Alphaproteobacteria<br />Species: Clade_II Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Clade_II Genus<br />group_subgroup: Alphaproteobacteria - Clade_II Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(61,96,122,1)","line":{"width":1.8897637795275593,"color":"rgba(61,96,122,1)"}},"name":"Alphaproteobacteria - Clade_II Genus","legendgroup":"Alphaproteobacteria - Clade_II Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000036,0.90000000000000013,0.90000000000000036],"base":[4416.5,9545,6297,3511],"x":[1,3,2,4],"y":[15,8,0,0],"text":["Class: Alphaproteobacteria<br />Species: Rhodobacteraceae Family<br />Sample: port_uptake_1<br />Abundance:   15.0000000<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family","Class: Alphaproteobacteria<br />Species: Rhodobacteraceae Family<br />Sample: ocean_uptake_1<br />Abundance:    8.0000000<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family","Class: Alphaproteobacteria<br />Species: Rhodobacteraceae Family<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family","Class: Alphaproteobacteria<br />Species: Rhodobacteraceae Family<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(82,116,141,1)","line":{"width":1.8897637795275593,"color":"rgba(82,116,141,1)"}},"name":"Alphaproteobacteria - Rhodobacteraceae Family","legendgroup":"Alphaproteobacteria - Rhodobacteraceae Family","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1758.1666666666667,5734],"x":[1,2],"y":[16.833333333333258,0],"text":["Class: Alphaproteobacteria<br />Species: Rhodobacteraceae Family<br />Sample: port_uptake_2<br />Abundance:   16.8333333<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family","Class: Alphaproteobacteria<br />Species: Rhodobacteraceae Family<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family<br />group_subgroup: Alphaproteobacteria - Rhodobacteraceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(82,116,141,1)","line":{"width":1.8897637795275593,"color":"rgba(82,116,141,1)"}},"name":"Alphaproteobacteria - Rhodobacteraceae Family","legendgroup":"Alphaproteobacteria - Rhodobacteraceae Family","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000013,0.89999999999999991,0.90000000000000036],"base":[3478,6294.333333333333,4416.5,9545],"x":[4,2,1,3],"y":[33,2.6666666666669698,0,0],"text":["Class: Alphaproteobacteria<br />Species: Roseovarius Genus<br />Sample: BWT_BWE_1<br />Abundance:   33.0000000<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus","Class: Alphaproteobacteria<br />Species: Roseovarius Genus<br />Sample: BWT_1<br />Abundance:    2.6666667<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus","Class: Alphaproteobacteria<br />Species: Roseovarius Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus","Class: Alphaproteobacteria<br />Species: Roseovarius Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(102,136,160,1)","line":{"width":1.8897637795275593,"color":"rgba(102,136,160,1)"}},"name":"Alphaproteobacteria - Roseovarius Genus","legendgroup":"Alphaproteobacteria - Roseovarius Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1758.1666666666667,5734],"x":[1,2],"y":[0,0],"text":["Class: Alphaproteobacteria<br />Species: Roseovarius Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus","Class: Alphaproteobacteria<br />Species: Roseovarius Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus<br />group_subgroup: Alphaproteobacteria - Roseovarius Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(102,136,160,1)","line":{"width":1.8897637795275593,"color":"rgba(102,136,160,1)"}},"name":"Alphaproteobacteria - Roseovarius Genus","legendgroup":"Alphaproteobacteria - Roseovarius Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000013,0.89999999999999991,0.90000000000000036],"base":[3441,6266.833333333333,4416.5,9545],"x":[4,2,1,3],"y":[37,27.5,0,0],"text":["Class: Alphaproteobacteria<br />Species: Sphingomonadaceae Family<br />Sample: BWT_BWE_1<br />Abundance:   37.0000000<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family","Class: Alphaproteobacteria<br />Species: Sphingomonadaceae Family<br />Sample: BWT_1<br />Abundance:   27.5000000<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family","Class: Alphaproteobacteria<br />Species: Sphingomonadaceae Family<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family","Class: Alphaproteobacteria<br />Species: Sphingomonadaceae Family<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(123,155,179,1)","line":{"width":1.8897637795275593,"color":"rgba(123,155,179,1)"}},"name":"Alphaproteobacteria - Sphingomonadaceae Family","legendgroup":"Alphaproteobacteria - Sphingomonadaceae Family","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1755.1666666666667,5734],"x":[1,2],"y":[3,0],"text":["Class: Alphaproteobacteria<br />Species: Sphingomonadaceae Family<br />Sample: port_uptake_2<br />Abundance:    3.0000000<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family","Class: Alphaproteobacteria<br />Species: Sphingomonadaceae Family<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family<br />group_subgroup: Alphaproteobacteria - Sphingomonadaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(123,155,179,1)","line":{"width":1.8897637795275593,"color":"rgba(123,155,179,1)"}},"name":"Alphaproteobacteria - Sphingomonadaceae Family","legendgroup":"Alphaproteobacteria - Sphingomonadaceae Family","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[4351.6666666666661,6266.833333333333,9545,3441],"x":[1,2,3,4],"y":[64.83333333333394,0,0,0],"text":["Class: Alphaproteobacteria<br />Species: Sphingorhabdus Genus<br />Sample: port_uptake_1<br />Abundance:   64.8333333<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus","Class: Alphaproteobacteria<br />Species: Sphingorhabdus Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus","Class: Alphaproteobacteria<br />Species: Sphingorhabdus Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus","Class: Alphaproteobacteria<br />Species: Sphingorhabdus Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(143,175,198,1)","line":{"width":1.8897637795275593,"color":"rgba(143,175,198,1)"}},"name":"Alphaproteobacteria - Sphingorhabdus Genus","legendgroup":"Alphaproteobacteria - Sphingorhabdus Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[5698,1746.3333333333333],"x":[2,1],"y":[36,8.8333333333334849],"text":["Class: Alphaproteobacteria<br />Species: Sphingorhabdus Genus<br />Sample: BWT_2<br />Abundance:   36.0000000<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus","Class: Alphaproteobacteria<br />Species: Sphingorhabdus Genus<br />Sample: port_uptake_2<br />Abundance:    8.8333333<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus<br />group_subgroup: Alphaproteobacteria - Sphingorhabdus Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(143,175,198,1)","line":{"width":1.8897637795275593,"color":"rgba(143,175,198,1)"}},"name":"Alphaproteobacteria - Sphingorhabdus Genus","legendgroup":"Alphaproteobacteria - Sphingorhabdus Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000013,0.89999999999999991,0.90000000000000036],"base":[3409.6666666666665,6264.666666666667,4351.6666666666661,9545],"x":[4,2,1,3],"y":[31.333333333333485,2.1666666666660603,0,0],"text":["Class: Alphaproteobacteria<br />Species: Sulfitobacter Genus<br />Sample: BWT_BWE_1<br />Abundance:   31.3333333<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus","Class: Alphaproteobacteria<br />Species: Sulfitobacter Genus<br />Sample: BWT_1<br />Abundance:    2.1666667<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus","Class: Alphaproteobacteria<br />Species: Sulfitobacter Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus","Class: Alphaproteobacteria<br />Species: Sulfitobacter Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(164,194,217,1)","line":{"width":1.8897637795275593,"color":"rgba(164,194,217,1)"}},"name":"Alphaproteobacteria - Sulfitobacter Genus","legendgroup":"Alphaproteobacteria - Sulfitobacter Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1746.3333333333333,5698],"x":[1,2],"y":[0,0],"text":["Class: Alphaproteobacteria<br />Species: Sulfitobacter Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus","Class: Alphaproteobacteria<br />Species: Sulfitobacter Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus<br />group_subgroup: Alphaproteobacteria - Sulfitobacter Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(164,194,217,1)","line":{"width":1.8897637795275593,"color":"rgba(164,194,217,1)"}},"name":"Alphaproteobacteria - Sulfitobacter Genus","legendgroup":"Alphaproteobacteria - Sulfitobacter Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[4351.6666666666661,6264.666666666667,9545,3409.6666666666665],"x":[1,2,3,4],"y":[0,0,0,0],"text":["Class: Alphaproteobacteria<br />Species: Tistrella Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Tistrella Genus<br />group_subgroup: Alphaproteobacteria - Tistrella Genus","Class: Alphaproteobacteria<br />Species: Tistrella Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Tistrella Genus<br />group_subgroup: Alphaproteobacteria - Tistrella Genus","Class: Alphaproteobacteria<br />Species: Tistrella Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Tistrella Genus<br />group_subgroup: Alphaproteobacteria - Tistrella Genus","Class: Alphaproteobacteria<br />Species: Tistrella Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Tistrella Genus<br />group_subgroup: Alphaproteobacteria - Tistrella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(184,214,236,1)","line":{"width":1.8897637795275593,"color":"rgba(184,214,236,1)"}},"name":"Alphaproteobacteria - Tistrella Genus","legendgroup":"Alphaproteobacteria - Tistrella Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1711.1666666666667,5698],"x":[1,2],"y":[35.166666666666515,0],"text":["Class: Alphaproteobacteria<br />Species: Tistrella Genus<br />Sample: port_uptake_2<br />Abundance:   35.1666667<br />group_subgroup: Alphaproteobacteria - Tistrella Genus<br />group_subgroup: Alphaproteobacteria - Tistrella Genus","Class: Alphaproteobacteria<br />Species: Tistrella Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Tistrella Genus<br />group_subgroup: Alphaproteobacteria - Tistrella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(184,214,236,1)","line":{"width":1.8897637795275593,"color":"rgba(184,214,236,1)"}},"name":"Alphaproteobacteria - Tistrella Genus","legendgroup":"Alphaproteobacteria - Tistrella Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000013,0.89999999999999991,0.90000000000000036],"base":[9458.6666666666661,6229.5,4329,3409.6666666666665],"x":[3,2,1,4],"y":[86.33333333333394,35.16666666666697,22.66666666666606,0],"text":["Class: Alphaproteobacteria<br />Species: Other Alphaproteobacteria<br />Sample: ocean_uptake_1<br />Abundance:   86.3333333<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria","Class: Alphaproteobacteria<br />Species: Other Alphaproteobacteria<br />Sample: BWT_1<br />Abundance:   35.1666667<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria","Class: Alphaproteobacteria<br />Species: Other Alphaproteobacteria<br />Sample: port_uptake_1<br />Abundance:   22.6666667<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria","Class: Alphaproteobacteria<br />Species: Other Alphaproteobacteria<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(205,234,255,1)","line":{"width":1.8897637795275593,"color":"rgba(205,234,255,1)"}},"name":"Alphaproteobacteria - Other Alphaproteobacteria","legendgroup":"Alphaproteobacteria - Other Alphaproteobacteria","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1577.5,5697],"x":[1,2],"y":[133.66666666666674,1],"text":["Class: Alphaproteobacteria<br />Species: Other Alphaproteobacteria<br />Sample: port_uptake_2<br />Abundance:  133.6666667<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria","Class: Alphaproteobacteria<br />Species: Other Alphaproteobacteria<br />Sample: BWT_2<br />Abundance:    1.0000000<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria<br />group_subgroup: Alphaproteobacteria - Other Alphaproteobacteria"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(205,234,255,1)","line":{"width":1.8897637795275593,"color":"rgba(205,234,255,1)"}},"name":"Alphaproteobacteria - Other Alphaproteobacteria","legendgroup":"Alphaproteobacteria - Other Alphaproteobacteria","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[4325.1666666666661,6229.5,9458.6666666666661,3409.6666666666665],"x":[1,2,3,4],"y":[3.8333333333339397,0,0,0],"text":["Class: Bacteroidia<br />Species: Antarctic_bacterium<br />Sample: port_uptake_1<br />Abundance:    3.8333333<br />group_subgroup: Bacteroidia - Antarctic_bacterium<br />group_subgroup: Bacteroidia - Antarctic_bacterium","Class: Bacteroidia<br />Species: Antarctic_bacterium<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Antarctic_bacterium<br />group_subgroup: Bacteroidia - Antarctic_bacterium","Class: Bacteroidia<br />Species: Antarctic_bacterium<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Antarctic_bacterium<br />group_subgroup: Bacteroidia - Antarctic_bacterium","Class: Bacteroidia<br />Species: Antarctic_bacterium<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Antarctic_bacterium<br />group_subgroup: Bacteroidia - Antarctic_bacterium"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(38,65,0,1)","line":{"width":1.8897637795275593,"color":"rgba(38,65,0,1)"}},"name":"Bacteroidia - Antarctic_bacterium","legendgroup":"Bacteroidia - Antarctic_bacterium","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[5596.5,1577.5],"x":[2,1],"y":[100.5,0],"text":["Class: Bacteroidia<br />Species: Antarctic_bacterium<br />Sample: BWT_2<br />Abundance:  100.5000000<br />group_subgroup: Bacteroidia - Antarctic_bacterium<br />group_subgroup: Bacteroidia - Antarctic_bacterium","Class: Bacteroidia<br />Species: Antarctic_bacterium<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Antarctic_bacterium<br />group_subgroup: Bacteroidia - Antarctic_bacterium"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(38,65,0,1)","line":{"width":1.8897637795275593,"color":"rgba(38,65,0,1)"}},"name":"Bacteroidia - Antarctic_bacterium","legendgroup":"Bacteroidia - Antarctic_bacterium","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[4325.1666666666661,6229.5,9458.6666666666661,3409.6666666666665],"x":[1,2,3,4],"y":[0,0,0,0],"text":["Class: Bacteroidia<br />Species: Arcicella Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Arcicella Genus<br />group_subgroup: Bacteroidia - Arcicella Genus","Class: Bacteroidia<br />Species: Arcicella Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Arcicella Genus<br />group_subgroup: Bacteroidia - Arcicella Genus","Class: Bacteroidia<br />Species: Arcicella Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Arcicella Genus<br />group_subgroup: Bacteroidia - Arcicella Genus","Class: Bacteroidia<br />Species: Arcicella Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Arcicella Genus<br />group_subgroup: Bacteroidia - Arcicella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(57,84,20,1)","line":{"width":1.8897637795275593,"color":"rgba(57,84,20,1)"}},"name":"Bacteroidia - Arcicella Genus","legendgroup":"Bacteroidia - Arcicella Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[5466.833333333333,1577.5],"x":[2,1],"y":[129.66666666666697,0],"text":["Class: Bacteroidia<br />Species: Arcicella Genus<br />Sample: BWT_2<br />Abundance:  129.6666667<br />group_subgroup: Bacteroidia - Arcicella Genus<br />group_subgroup: Bacteroidia - Arcicella Genus","Class: Bacteroidia<br />Species: Arcicella Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Arcicella Genus<br />group_subgroup: Bacteroidia - Arcicella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(57,84,20,1)","line":{"width":1.8897637795275593,"color":"rgba(57,84,20,1)"}},"name":"Bacteroidia - Arcicella Genus","legendgroup":"Bacteroidia - Arcicella Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[9420,4305.333333333333,6229.5,3409.6666666666665],"x":[3,1,2,4],"y":[38.66666666666606,19.83333333333303,0,0],"text":["Class: Bacteroidia<br />Species: Flavobacteriaceae Family<br />Sample: ocean_uptake_1<br />Abundance:   38.6666667<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family","Class: Bacteroidia<br />Species: Flavobacteriaceae Family<br />Sample: port_uptake_1<br />Abundance:   19.8333333<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family","Class: Bacteroidia<br />Species: Flavobacteriaceae Family<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family","Class: Bacteroidia<br />Species: Flavobacteriaceae Family<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(77,103,41,1)","line":{"width":1.8897637795275593,"color":"rgba(77,103,41,1)"}},"name":"Bacteroidia - Flavobacteriaceae Family","legendgroup":"Bacteroidia - Flavobacteriaceae Family","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1577.5,5466.833333333333],"x":[1,2],"y":[0,0],"text":["Class: Bacteroidia<br />Species: Flavobacteriaceae Family<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family","Class: Bacteroidia<br />Species: Flavobacteriaceae Family<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family<br />group_subgroup: Bacteroidia - Flavobacteriaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(77,103,41,1)","line":{"width":1.8897637795275593,"color":"rgba(77,103,41,1)"}},"name":"Bacteroidia - Flavobacteriaceae Family","legendgroup":"Bacteroidia - Flavobacteriaceae Family","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[3739.333333333333,6206.666666666667,9420,3409.6666666666665],"x":[1,2,3,4],"y":[566,22.83333333333303,0,0],"text":["Class: Bacteroidia<br />Species: Flavobacterium Genus<br />Sample: port_uptake_1<br />Abundance:  566.0000000<br />group_subgroup: Bacteroidia - Flavobacterium Genus<br />group_subgroup: Bacteroidia - Flavobacterium Genus","Class: Bacteroidia<br />Species: Flavobacterium Genus<br />Sample: BWT_1<br />Abundance:   22.8333333<br />group_subgroup: Bacteroidia - Flavobacterium Genus<br />group_subgroup: Bacteroidia - Flavobacterium Genus","Class: Bacteroidia<br />Species: Flavobacterium Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacterium Genus<br />group_subgroup: Bacteroidia - Flavobacterium Genus","Class: Bacteroidia<br />Species: Flavobacterium Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacterium Genus<br />group_subgroup: Bacteroidia - Flavobacterium Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(96,122,61,1)","line":{"width":1.8897637795275593,"color":"rgba(96,122,61,1)"}},"name":"Bacteroidia - Flavobacterium Genus","legendgroup":"Bacteroidia - Flavobacterium Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[4599.333333333333,1452.6666666666667],"x":[2,1],"y":[867.5,124.83333333333326],"text":["Class: Bacteroidia<br />Species: Flavobacterium Genus<br />Sample: BWT_2<br />Abundance:  867.5000000<br />group_subgroup: Bacteroidia - Flavobacterium Genus<br />group_subgroup: Bacteroidia - Flavobacterium Genus","Class: Bacteroidia<br />Species: Flavobacterium Genus<br />Sample: port_uptake_2<br />Abundance:  124.8333333<br />group_subgroup: Bacteroidia - Flavobacterium Genus<br />group_subgroup: Bacteroidia - Flavobacterium Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(96,122,61,1)","line":{"width":1.8897637795275593,"color":"rgba(96,122,61,1)"}},"name":"Bacteroidia - Flavobacterium Genus","legendgroup":"Bacteroidia - Flavobacterium Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[9404,3739.333333333333,6206.666666666667,3409.6666666666665],"x":[3,1,2,4],"y":[16,0,0,0],"text":["Class: Bacteroidia<br />Species: Flavobacterium_jumunjinense<br />Sample: ocean_uptake_1<br />Abundance:   16.0000000<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense","Class: Bacteroidia<br />Species: Flavobacterium_jumunjinense<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense","Class: Bacteroidia<br />Species: Flavobacterium_jumunjinense<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense","Class: Bacteroidia<br />Species: Flavobacterium_jumunjinense<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(116,141,82,1)","line":{"width":1.8897637795275593,"color":"rgba(116,141,82,1)"}},"name":"Bacteroidia - Flavobacterium_jumunjinense","legendgroup":"Bacteroidia - Flavobacterium_jumunjinense","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1452.6666666666667,4599.333333333333],"x":[1,2],"y":[0,0],"text":["Class: Bacteroidia<br />Species: Flavobacterium_jumunjinense<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense","Class: Bacteroidia<br />Species: Flavobacterium_jumunjinense<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense<br />group_subgroup: Bacteroidia - Flavobacterium_jumunjinense"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(116,141,82,1)","line":{"width":1.8897637795275593,"color":"rgba(116,141,82,1)"}},"name":"Bacteroidia - Flavobacterium_jumunjinense","legendgroup":"Bacteroidia - Flavobacterium_jumunjinense","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991,0.90000000000000036,0.90000000000000036],"base":[6202.666666666667,3739.333333333333,9404,3409.6666666666665],"x":[2,1,3,4],"y":[4,0,0,0],"text":["Class: Bacteroidia<br />Species: Flectobacillus Genus<br />Sample: BWT_1<br />Abundance:    4.0000000<br />group_subgroup: Bacteroidia - Flectobacillus Genus<br />group_subgroup: Bacteroidia - Flectobacillus Genus","Class: Bacteroidia<br />Species: Flectobacillus Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flectobacillus Genus<br />group_subgroup: Bacteroidia - Flectobacillus Genus","Class: Bacteroidia<br />Species: Flectobacillus Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flectobacillus Genus<br />group_subgroup: Bacteroidia - Flectobacillus Genus","Class: Bacteroidia<br />Species: Flectobacillus Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flectobacillus Genus<br />group_subgroup: Bacteroidia - Flectobacillus Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(136,160,102,1)","line":{"width":1.8897637795275593,"color":"rgba(136,160,102,1)"}},"name":"Bacteroidia - Flectobacillus Genus","legendgroup":"Bacteroidia - Flectobacillus Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1452.6666666666667,4599.333333333333],"x":[1,2],"y":[0,0],"text":["Class: Bacteroidia<br />Species: Flectobacillus Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flectobacillus Genus<br />group_subgroup: Bacteroidia - Flectobacillus Genus","Class: Bacteroidia<br />Species: Flectobacillus Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Flectobacillus Genus<br />group_subgroup: Bacteroidia - Flectobacillus Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(136,160,102,1)","line":{"width":1.8897637795275593,"color":"rgba(136,160,102,1)"}},"name":"Bacteroidia - Flectobacillus Genus","legendgroup":"Bacteroidia - Flectobacillus Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[3729.333333333333,6202.666666666667,9404,3409.6666666666665],"x":[1,2,3,4],"y":[10,0,0,0],"text":["Class: Bacteroidia<br />Species: Fluviicola Genus<br />Sample: port_uptake_1<br />Abundance:   10.0000000<br />group_subgroup: Bacteroidia - Fluviicola Genus<br />group_subgroup: Bacteroidia - Fluviicola Genus","Class: Bacteroidia<br />Species: Fluviicola Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Fluviicola Genus<br />group_subgroup: Bacteroidia - Fluviicola Genus","Class: Bacteroidia<br />Species: Fluviicola Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Fluviicola Genus<br />group_subgroup: Bacteroidia - Fluviicola Genus","Class: Bacteroidia<br />Species: Fluviicola Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Fluviicola Genus<br />group_subgroup: Bacteroidia - Fluviicola Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(155,179,123,1)","line":{"width":1.8897637795275593,"color":"rgba(155,179,123,1)"}},"name":"Bacteroidia - Fluviicola Genus","legendgroup":"Bacteroidia - Fluviicola Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1381.6666666666667,4599.333333333333],"x":[1,2],"y":[71,0],"text":["Class: Bacteroidia<br />Species: Fluviicola Genus<br />Sample: port_uptake_2<br />Abundance:   71.0000000<br />group_subgroup: Bacteroidia - Fluviicola Genus<br />group_subgroup: Bacteroidia - Fluviicola Genus","Class: Bacteroidia<br />Species: Fluviicola Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Fluviicola Genus<br />group_subgroup: Bacteroidia - Fluviicola Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(155,179,123,1)","line":{"width":1.8897637795275593,"color":"rgba(155,179,123,1)"}},"name":"Bacteroidia - Fluviicola Genus","legendgroup":"Bacteroidia - Fluviicola Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[9384.3333333333339,3729.333333333333,6202.666666666667,3409.6666666666665],"x":[3,1,2,4],"y":[19.66666666666606,0,0,0],"text":["Class: Bacteroidia<br />Species: Mesoflavibacter_sp.<br />Sample: ocean_uptake_1<br />Abundance:   19.6666667<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.","Class: Bacteroidia<br />Species: Mesoflavibacter_sp.<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.","Class: Bacteroidia<br />Species: Mesoflavibacter_sp.<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.","Class: Bacteroidia<br />Species: Mesoflavibacter_sp.<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp."],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(175,198,143,1)","line":{"width":1.8897637795275593,"color":"rgba(175,198,143,1)"}},"name":"Bacteroidia - Mesoflavibacter_sp.","legendgroup":"Bacteroidia - Mesoflavibacter_sp.","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1381.6666666666667,4599.333333333333],"x":[1,2],"y":[0,0],"text":["Class: Bacteroidia<br />Species: Mesoflavibacter_sp.<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.","Class: Bacteroidia<br />Species: Mesoflavibacter_sp.<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp.<br />group_subgroup: Bacteroidia - Mesoflavibacter_sp."],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(175,198,143,1)","line":{"width":1.8897637795275593,"color":"rgba(175,198,143,1)"}},"name":"Bacteroidia - Mesoflavibacter_sp.","legendgroup":"Bacteroidia - Mesoflavibacter_sp.","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[3713.333333333333,6202.666666666667,9384.3333333333339,3409.6666666666665],"x":[1,2,3,4],"y":[16,0,0,0],"text":["Class: Bacteroidia<br />Species: Pseudarcicella Genus<br />Sample: port_uptake_1<br />Abundance:   16.0000000<br />group_subgroup: Bacteroidia - Pseudarcicella Genus<br />group_subgroup: Bacteroidia - Pseudarcicella Genus","Class: Bacteroidia<br />Species: Pseudarcicella Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Pseudarcicella Genus<br />group_subgroup: Bacteroidia - Pseudarcicella Genus","Class: Bacteroidia<br />Species: Pseudarcicella Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Pseudarcicella Genus<br />group_subgroup: Bacteroidia - Pseudarcicella Genus","Class: Bacteroidia<br />Species: Pseudarcicella Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Pseudarcicella Genus<br />group_subgroup: Bacteroidia - Pseudarcicella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(194,217,164,1)","line":{"width":1.8897637795275593,"color":"rgba(194,217,164,1)"}},"name":"Bacteroidia - Pseudarcicella Genus","legendgroup":"Bacteroidia - Pseudarcicella Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1259.5,4599.333333333333],"x":[1,2],"y":[122.16666666666674,0],"text":["Class: Bacteroidia<br />Species: Pseudarcicella Genus<br />Sample: port_uptake_2<br />Abundance:  122.1666667<br />group_subgroup: Bacteroidia - Pseudarcicella Genus<br />group_subgroup: Bacteroidia - Pseudarcicella Genus","Class: Bacteroidia<br />Species: Pseudarcicella Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Pseudarcicella Genus<br />group_subgroup: Bacteroidia - Pseudarcicella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(194,217,164,1)","line":{"width":1.8897637795275593,"color":"rgba(194,217,164,1)"}},"name":"Bacteroidia - Pseudarcicella Genus","legendgroup":"Bacteroidia - Pseudarcicella Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[3638.333333333333,6202.333333333333,9384.3333333333339,3409.6666666666665],"x":[1,2,3,4],"y":[75,0.33333333333393966,0,0],"text":["Class: Bacteroidia<br />Species: uncultured_bacterium<br />Sample: port_uptake_1<br />Abundance:   75.0000000<br />group_subgroup: Bacteroidia - uncultured_bacterium<br />group_subgroup: Bacteroidia - uncultured_bacterium","Class: Bacteroidia<br />Species: uncultured_bacterium<br />Sample: BWT_1<br />Abundance:    0.3333333<br />group_subgroup: Bacteroidia - uncultured_bacterium<br />group_subgroup: Bacteroidia - uncultured_bacterium","Class: Bacteroidia<br />Species: uncultured_bacterium<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - uncultured_bacterium<br />group_subgroup: Bacteroidia - uncultured_bacterium","Class: Bacteroidia<br />Species: uncultured_bacterium<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - uncultured_bacterium<br />group_subgroup: Bacteroidia - uncultured_bacterium"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(214,236,184,1)","line":{"width":1.8897637795275593,"color":"rgba(214,236,184,1)"}},"name":"Bacteroidia - uncultured_bacterium","legendgroup":"Bacteroidia - uncultured_bacterium","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[1259.5,4599.333333333333],"x":[1,2],"y":[0,0],"text":["Class: Bacteroidia<br />Species: uncultured_bacterium<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - uncultured_bacterium<br />group_subgroup: Bacteroidia - uncultured_bacterium","Class: Bacteroidia<br />Species: uncultured_bacterium<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - uncultured_bacterium<br />group_subgroup: Bacteroidia - uncultured_bacterium"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(214,236,184,1)","line":{"width":1.8897637795275593,"color":"rgba(214,236,184,1)"}},"name":"Bacteroidia - uncultured_bacterium","legendgroup":"Bacteroidia - uncultured_bacterium","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000036,0.90000000000000013,0.90000000000000036],"base":[3500,9288.6666666666661,6197.666666666667,3409.6666666666665],"x":[1,3,2,4],"y":[138.33333333333303,95.666666666667879,4.6666666666660603,0],"text":["Class: Bacteroidia<br />Species: Other Bacteroidia<br />Sample: port_uptake_1<br />Abundance:  138.3333333<br />group_subgroup: Bacteroidia - Other Bacteroidia<br />group_subgroup: Bacteroidia - Other Bacteroidia","Class: Bacteroidia<br />Species: Other Bacteroidia<br />Sample: ocean_uptake_1<br />Abundance:   95.6666667<br />group_subgroup: Bacteroidia - Other Bacteroidia<br />group_subgroup: Bacteroidia - Other Bacteroidia","Class: Bacteroidia<br />Species: Other Bacteroidia<br />Sample: BWT_1<br />Abundance:    4.6666667<br />group_subgroup: Bacteroidia - Other Bacteroidia<br />group_subgroup: Bacteroidia - Other Bacteroidia","Class: Bacteroidia<br />Species: Other Bacteroidia<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Bacteroidia - Other Bacteroidia<br />group_subgroup: Bacteroidia - Other Bacteroidia"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(234,255,205,1)","line":{"width":1.8897637795275593,"color":"rgba(234,255,205,1)"}},"name":"Bacteroidia - Other Bacteroidia","legendgroup":"Bacteroidia - Other Bacteroidia","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[821.16666666666663,4568.333333333333],"x":[1,2],"y":[438.33333333333337,31],"text":["Class: Bacteroidia<br />Species: Other Bacteroidia<br />Sample: port_uptake_2<br />Abundance:  438.3333333<br />group_subgroup: Bacteroidia - Other Bacteroidia<br />group_subgroup: Bacteroidia - Other Bacteroidia","Class: Bacteroidia<br />Species: Other Bacteroidia<br />Sample: BWT_2<br />Abundance:   31.0000000<br />group_subgroup: Bacteroidia - Other Bacteroidia<br />group_subgroup: Bacteroidia - Other Bacteroidia"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(234,255,205,1)","line":{"width":1.8897637795275593,"color":"rgba(234,255,205,1)"}},"name":"Bacteroidia - Other Bacteroidia","legendgroup":"Bacteroidia - Other Bacteroidia","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[911.33333333333326,6197.666666666667,9288.6666666666661,3409.6666666666665],"x":[1,2,3,4],"y":[2588.666666666667,0,0,0],"text":["Class: Gammaproteobacteria<br />Species: Aeromonas Genus<br />Sample: port_uptake_1<br />Abundance: 2588.6666667<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus","Class: Gammaproteobacteria<br />Species: Aeromonas Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus","Class: Gammaproteobacteria<br />Species: Aeromonas Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus","Class: Gammaproteobacteria<br />Species: Aeromonas Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(65,0,38,1)","line":{"width":1.8897637795275593,"color":"rgba(65,0,38,1)"}},"name":"Gammaproteobacteria - Aeromonas Genus","legendgroup":"Gammaproteobacteria - Aeromonas Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[821.16666666666663,4568.333333333333],"x":[1,2],"y":[0,0],"text":["Class: Gammaproteobacteria<br />Species: Aeromonas Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus","Class: Gammaproteobacteria<br />Species: Aeromonas Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus<br />group_subgroup: Gammaproteobacteria - Aeromonas Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(65,0,38,1)","line":{"width":1.8897637795275593,"color":"rgba(65,0,38,1)"}},"name":"Gammaproteobacteria - Aeromonas Genus","legendgroup":"Gammaproteobacteria - Aeromonas Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991,0.90000000000000036,0.90000000000000036],"base":[6047.166666666667,831.66666666666663,9288.6666666666661,3409.6666666666665],"x":[2,1,3,4],"y":[150.5,79.666666666666629,0,0],"text":["Class: Gammaproteobacteria<br />Species: Comamonadaceae Family<br />Sample: BWT_1<br />Abundance:  150.5000000<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family","Class: Gammaproteobacteria<br />Species: Comamonadaceae Family<br />Sample: port_uptake_1<br />Abundance:   79.6666667<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family","Class: Gammaproteobacteria<br />Species: Comamonadaceae Family<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family","Class: Gammaproteobacteria<br />Species: Comamonadaceae Family<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(84,20,57,1)","line":{"width":1.8897637795275593,"color":"rgba(84,20,57,1)"}},"name":"Gammaproteobacteria - Comamonadaceae Family","legendgroup":"Gammaproteobacteria - Comamonadaceae Family","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[495.66666666666669,4450.833333333333],"x":[1,2],"y":[325.49999999999994,117.5],"text":["Class: Gammaproteobacteria<br />Species: Comamonadaceae Family<br />Sample: port_uptake_2<br />Abundance:  325.5000000<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family","Class: Gammaproteobacteria<br />Species: Comamonadaceae Family<br />Sample: BWT_2<br />Abundance:  117.5000000<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family<br />group_subgroup: Gammaproteobacteria - Comamonadaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(84,20,57,1)","line":{"width":1.8897637795275593,"color":"rgba(84,20,57,1)"}},"name":"Gammaproteobacteria - Comamonadaceae Family","legendgroup":"Gammaproteobacteria - Comamonadaceae Family","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.89999999999999991,0.90000000000000013,0.90000000000000036],"base":[7681,831.66666666666663,6047.166666666667,3409.6666666666665],"x":[3,1,2,4],"y":[1607.6666666666661,0,0,0],"text":["Class: Gammaproteobacteria<br />Species: Oleispira Genus<br />Sample: ocean_uptake_1<br />Abundance: 1607.6666667<br />group_subgroup: Gammaproteobacteria - Oleispira Genus<br />group_subgroup: Gammaproteobacteria - Oleispira Genus","Class: Gammaproteobacteria<br />Species: Oleispira Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Oleispira Genus<br />group_subgroup: Gammaproteobacteria - Oleispira Genus","Class: Gammaproteobacteria<br />Species: Oleispira Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Oleispira Genus<br />group_subgroup: Gammaproteobacteria - Oleispira Genus","Class: Gammaproteobacteria<br />Species: Oleispira Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Oleispira Genus<br />group_subgroup: Gammaproteobacteria - Oleispira Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(103,41,77,1)","line":{"width":1.8897637795275593,"color":"rgba(103,41,77,1)"}},"name":"Gammaproteobacteria - Oleispira Genus","legendgroup":"Gammaproteobacteria - Oleispira Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[495.66666666666669,4450.833333333333],"x":[1,2],"y":[0,0],"text":["Class: Gammaproteobacteria<br />Species: Oleispira Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Oleispira Genus<br />group_subgroup: Gammaproteobacteria - Oleispira Genus","Class: Gammaproteobacteria<br />Species: Oleispira Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Oleispira Genus<br />group_subgroup: Gammaproteobacteria - Oleispira Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(103,41,77,1)","line":{"width":1.8897637795275593,"color":"rgba(103,41,77,1)"}},"name":"Gammaproteobacteria - Oleispira Genus","legendgroup":"Gammaproteobacteria - Oleispira Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013,0.90000000000000036,0.90000000000000036],"base":[831.66666666666663,6047.166666666667,7681,3409.6666666666665],"x":[1,2,3,4],"y":[0,0,0,0],"text":["Class: Gammaproteobacteria<br />Species: Polynucleobacter Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus","Class: Gammaproteobacteria<br />Species: Polynucleobacter Genus<br />Sample: BWT_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus","Class: Gammaproteobacteria<br />Species: Polynucleobacter Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus","Class: Gammaproteobacteria<br />Species: Polynucleobacter Genus<br />Sample: BWT_BWE_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(122,61,96,1)","line":{"width":1.8897637795275593,"color":"rgba(122,61,96,1)"}},"name":"Gammaproteobacteria - Polynucleobacter Genus","legendgroup":"Gammaproteobacteria - Polynucleobacter Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[397,4450.833333333333],"x":[1,2],"y":[98.666666666666686,0],"text":["Class: Gammaproteobacteria<br />Species: Polynucleobacter Genus<br />Sample: port_uptake_2<br />Abundance:   98.6666667<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus","Class: Gammaproteobacteria<br />Species: Polynucleobacter Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus<br />group_subgroup: Gammaproteobacteria - Polynucleobacter Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(122,61,96,1)","line":{"width":1.8897637795275593,"color":"rgba(122,61,96,1)"}},"name":"Gammaproteobacteria - Polynucleobacter Genus","legendgroup":"Gammaproteobacteria - Polynucleobacter Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036,0.90000000000000013,0.89999999999999991],"base":[1922.6666666666665,6463.6666666666661,6042.833333333333,831.66666666666663],"x":[4,3,2,1],"y":[1487,1217.3333333333339,4.3333333333339397,0],"text":["Class: Gammaproteobacteria<br />Species: Pseudoalteromonas Genus<br />Sample: BWT_BWE_1<br />Abundance: 1487.0000000<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus","Class: Gammaproteobacteria<br />Species: Pseudoalteromonas Genus<br />Sample: ocean_uptake_1<br />Abundance: 1217.3333333<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus","Class: Gammaproteobacteria<br />Species: Pseudoalteromonas Genus<br />Sample: BWT_1<br />Abundance:    4.3333333<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus","Class: Gammaproteobacteria<br />Species: Pseudoalteromonas Genus<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(141,82,116,1)","line":{"width":1.8897637795275593,"color":"rgba(141,82,116,1)"}},"name":"Gammaproteobacteria - Pseudoalteromonas Genus","legendgroup":"Gammaproteobacteria - Pseudoalteromonas Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[397,4450.833333333333],"x":[1,2],"y":[0,0],"text":["Class: Gammaproteobacteria<br />Species: Pseudoalteromonas Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus","Class: Gammaproteobacteria<br />Species: Pseudoalteromonas Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudoalteromonas Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(141,82,116,1)","line":{"width":1.8897637795275593,"color":"rgba(141,82,116,1)"}},"name":"Gammaproteobacteria - Pseudoalteromonas Genus","legendgroup":"Gammaproteobacteria - Pseudoalteromonas Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.89999999999999991,0.90000000000000036],"base":[3875.5,706.33333333333337,673.66666666666663,6431],"x":[2,4,1,3],"y":[2167.333333333333,1216.333333333333,158,32.66666666666606],"text":["Class: Gammaproteobacteria<br />Species: Pseudomonas Genus<br />Sample: BWT_1<br />Abundance: 2167.3333333<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus","Class: Gammaproteobacteria<br />Species: Pseudomonas Genus<br />Sample: BWT_BWE_1<br />Abundance: 1216.3333333<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus","Class: Gammaproteobacteria<br />Species: Pseudomonas Genus<br />Sample: port_uptake_1<br />Abundance:  158.0000000<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus","Class: Gammaproteobacteria<br />Species: Pseudomonas Genus<br />Sample: ocean_uptake_1<br />Abundance:   32.6666667<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(160,102,136,1)","line":{"width":1.8897637795275593,"color":"rgba(160,102,136,1)"}},"name":"Gammaproteobacteria - Pseudomonas Genus","legendgroup":"Gammaproteobacteria - Pseudomonas Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[1164.8333333333333,397],"x":[2,1],"y":[3286,0],"text":["Class: Gammaproteobacteria<br />Species: Pseudomonas Genus<br />Sample: BWT_2<br />Abundance: 3286.0000000<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus","Class: Gammaproteobacteria<br />Species: Pseudomonas Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus<br />group_subgroup: Gammaproteobacteria - Pseudomonas Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(160,102,136,1)","line":{"width":1.8897637795275593,"color":"rgba(160,102,136,1)"}},"name":"Gammaproteobacteria - Pseudomonas Genus","legendgroup":"Gammaproteobacteria - Pseudomonas Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991,0.90000000000000036,0.90000000000000036],"base":[2064,395.33333333333331,675,6428.6666666666661],"x":[2,1,4,3],"y":[1811.5,278.33333333333331,31.333333333333371,2.3333333333339397],"text":["Class: Gammaproteobacteria<br />Species: Pseudomonas_anguilliseptica<br />Sample: BWT_1<br />Abundance: 1811.5000000<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica","Class: Gammaproteobacteria<br />Species: Pseudomonas_anguilliseptica<br />Sample: port_uptake_1<br />Abundance:  278.3333333<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica","Class: Gammaproteobacteria<br />Species: Pseudomonas_anguilliseptica<br />Sample: BWT_BWE_1<br />Abundance:   31.3333333<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica","Class: Gammaproteobacteria<br />Species: Pseudomonas_anguilliseptica<br />Sample: ocean_uptake_1<br />Abundance:    2.3333333<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(179,123,155,1)","line":{"width":1.8897637795275593,"color":"rgba(179,123,155,1)"}},"name":"Gammaproteobacteria - Pseudomonas_anguilliseptica","legendgroup":"Gammaproteobacteria - Pseudomonas_anguilliseptica","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[1122.5,397],"x":[2,1],"y":[42.333333333333258,0],"text":["Class: Gammaproteobacteria<br />Species: Pseudomonas_anguilliseptica<br />Sample: BWT_2<br />Abundance:   42.3333333<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica","Class: Gammaproteobacteria<br />Species: Pseudomonas_anguilliseptica<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica<br />group_subgroup: Gammaproteobacteria - Pseudomonas_anguilliseptica"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(179,123,155,1)","line":{"width":1.8897637795275593,"color":"rgba(179,123,155,1)"}},"name":"Gammaproteobacteria - Pseudomonas_anguilliseptica","legendgroup":"Gammaproteobacteria - Pseudomonas_anguilliseptica","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.89999999999999991,0.90000000000000036],"base":[877.66666666666652,454,363.33333333333331,6428.6666666666661],"x":[2,4,1,3],"y":[1186.3333333333335,221,32,0],"text":["Class: Gammaproteobacteria<br />Species: Rheinheimera Genus<br />Sample: BWT_1<br />Abundance: 1186.3333333<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus","Class: Gammaproteobacteria<br />Species: Rheinheimera Genus<br />Sample: BWT_BWE_1<br />Abundance:  221.0000000<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus","Class: Gammaproteobacteria<br />Species: Rheinheimera Genus<br />Sample: port_uptake_1<br />Abundance:   32.0000000<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus","Class: Gammaproteobacteria<br />Species: Rheinheimera Genus<br />Sample: ocean_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(198,143,175,1)","line":{"width":1.8897637795275593,"color":"rgba(198,143,175,1)"}},"name":"Gammaproteobacteria - Rheinheimera Genus","legendgroup":"Gammaproteobacteria - Rheinheimera Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[1099.3333333333333,397],"x":[2,1],"y":[23.166666666666742,0],"text":["Class: Gammaproteobacteria<br />Species: Rheinheimera Genus<br />Sample: BWT_2<br />Abundance:   23.1666667<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus","Class: Gammaproteobacteria<br />Species: Rheinheimera Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus<br />group_subgroup: Gammaproteobacteria - Rheinheimera Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(198,143,175,1)","line":{"width":1.8897637795275593,"color":"rgba(198,143,175,1)"}},"name":"Gammaproteobacteria - Rheinheimera Genus","legendgroup":"Gammaproteobacteria - Rheinheimera Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.90000000000000036,0.90000000000000036,0.89999999999999991],"base":[617.33333333333326,6200.6666666666661,330.33333333333337,299.83333333333331],"x":[2,3,4,1],"y":[260.33333333333326,228,123.66666666666663,63.5],"text":["Class: Gammaproteobacteria<br />Species: Shewanella Genus<br />Sample: BWT_1<br />Abundance:  260.3333333<br />group_subgroup: Gammaproteobacteria - Shewanella Genus<br />group_subgroup: Gammaproteobacteria - Shewanella Genus","Class: Gammaproteobacteria<br />Species: Shewanella Genus<br />Sample: ocean_uptake_1<br />Abundance:  228.0000000<br />group_subgroup: Gammaproteobacteria - Shewanella Genus<br />group_subgroup: Gammaproteobacteria - Shewanella Genus","Class: Gammaproteobacteria<br />Species: Shewanella Genus<br />Sample: BWT_BWE_1<br />Abundance:  123.6666667<br />group_subgroup: Gammaproteobacteria - Shewanella Genus<br />group_subgroup: Gammaproteobacteria - Shewanella Genus","Class: Gammaproteobacteria<br />Species: Shewanella Genus<br />Sample: port_uptake_1<br />Abundance:   63.5000000<br />group_subgroup: Gammaproteobacteria - Shewanella Genus<br />group_subgroup: Gammaproteobacteria - Shewanella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(217,164,194,1)","line":{"width":1.8897637795275593,"color":"rgba(217,164,194,1)"}},"name":"Gammaproteobacteria - Shewanella Genus","legendgroup":"Gammaproteobacteria - Shewanella Genus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[397,1099.3333333333333],"x":[1,2],"y":[0,0],"text":["Class: Gammaproteobacteria<br />Species: Shewanella Genus<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Shewanella Genus<br />group_subgroup: Gammaproteobacteria - Shewanella Genus","Class: Gammaproteobacteria<br />Species: Shewanella Genus<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Shewanella Genus<br />group_subgroup: Gammaproteobacteria - Shewanella Genus"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(217,164,194,1)","line":{"width":1.8897637795275593,"color":"rgba(217,164,194,1)"}},"name":"Gammaproteobacteria - Shewanella Genus","legendgroup":"Gammaproteobacteria - Shewanella Genus","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000036,0.90000000000000013,0.89999999999999991],"base":[1922.3333333333333,298.66666666666669,607.83333333333326,299.83333333333331],"x":[3,4,2,1],"y":[4278.333333333333,31.666666666666686,9.5,0],"text":["Class: Gammaproteobacteria<br />Species: Vibrionaceae Family<br />Sample: ocean_uptake_1<br />Abundance: 4278.3333333<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family","Class: Gammaproteobacteria<br />Species: Vibrionaceae Family<br />Sample: BWT_BWE_1<br />Abundance:   31.6666667<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family","Class: Gammaproteobacteria<br />Species: Vibrionaceae Family<br />Sample: BWT_1<br />Abundance:    9.5000000<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family","Class: Gammaproteobacteria<br />Species: Vibrionaceae Family<br />Sample: port_uptake_1<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(236,184,214,1)","line":{"width":1.8897637795275593,"color":"rgba(236,184,214,1)"}},"name":"Gammaproteobacteria - Vibrionaceae Family","legendgroup":"Gammaproteobacteria - Vibrionaceae Family","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.89999999999999991,0.90000000000000013],"base":[397,1099.3333333333333],"x":[1,2],"y":[0,0],"text":["Class: Gammaproteobacteria<br />Species: Vibrionaceae Family<br />Sample: port_uptake_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family","Class: Gammaproteobacteria<br />Species: Vibrionaceae Family<br />Sample: BWT_2<br />Abundance:    0.0000000<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family<br />group_subgroup: Gammaproteobacteria - Vibrionaceae Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(236,184,214,1)","line":{"width":1.8897637795275593,"color":"rgba(236,184,214,1)"}},"name":"Gammaproteobacteria - Vibrionaceae Family","legendgroup":"Gammaproteobacteria - Vibrionaceae Family","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000036,0.90000000000000013,0.89999999999999991,0.90000000000000036],"base":[0,0,0,0],"x":[3,2,1,4],"y":[1922.3333333333333,607.83333333333326,299.83333333333331,298.66666666666669],"text":["Class: Gammaproteobacteria<br />Species: Other Gammaproteobacteria<br />Sample: ocean_uptake_1<br />Abundance: 1922.3333333<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria","Class: Gammaproteobacteria<br />Species: Other Gammaproteobacteria<br />Sample: BWT_1<br />Abundance:  607.8333333<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria","Class: Gammaproteobacteria<br />Species: Other Gammaproteobacteria<br />Sample: port_uptake_1<br />Abundance:  299.8333333<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria","Class: Gammaproteobacteria<br />Species: Other Gammaproteobacteria<br />Sample: BWT_BWE_1<br />Abundance:  298.6666667<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(255,205,234,1)","line":{"width":1.8897637795275593,"color":"rgba(255,205,234,1)"}},"name":"Gammaproteobacteria - Other Gammaproteobacteria","legendgroup":"Gammaproteobacteria - Other Gammaproteobacteria","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.90000000000000013,0.89999999999999991],"base":[0,0],"x":[2,1],"y":[1099.3333333333333,397],"text":["Class: Gammaproteobacteria<br />Species: Other Gammaproteobacteria<br />Sample: BWT_2<br />Abundance: 1099.3333333<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria","Class: Gammaproteobacteria<br />Species: Other Gammaproteobacteria<br />Sample: port_uptake_2<br />Abundance:  397.0000000<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria<br />group_subgroup: Gammaproteobacteria - Other Gammaproteobacteria"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(255,205,234,1)","line":{"width":1.8897637795275593,"color":"rgba(255,205,234,1)"}},"name":"Gammaproteobacteria - Other Gammaproteobacteria","legendgroup":"Gammaproteobacteria - Other Gammaproteobacteria","showlegend":false,"xaxis":"x2","yaxis":"y2","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":39.246160232461605,"r":7.3059360730593621,"b":126.52552926525527,"l":53.466168534661698},"paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,0.46204708533475658],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,4.5999999999999996],"tickmode":"array","ticktext":["Port uptake V1","BWT V1","Ocean uptake V1","BWT+BWE V1"],"tickvals":[1,2,3,4],"categoryorder":"array","categoryarray":["Port uptake V1","BWT V1","Ocean uptake V1","BWT+BWE V1"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-90,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"annotations":[{"text":"Sample type","x":0.5,"y":0,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"top","annotationType":"axis","yshift":-109.85471149854712},{"text":"Abundance","x":0,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis","xshift":-36.795350767953508},{"text":"Voyage 1","x":0.23102354266737829,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"},{"text":"Voyage 2","x":0.76897645733262165,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(26,26,26,1)","family":"","size":15.940224159402243},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-550,11550],"tickmode":"array","ticktext":["0","3000","6000","9000"],"tickvals":[0,3000.0000000000005,6000,9000],"categoryorder":"array","categoryarray":["0","3000","6000","9000"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.13283520132835203,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835205},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.46204708533475658,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":0.46204708533475658,"y0":0,"y1":27.629721876297225,"yanchor":1,"ysizemode":"pixel"},{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0.53795291466524342,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","layer":"below","x0":0.53795291466524342,"x1":1,"y0":0,"y1":27.629721876297225,"yanchor":1,"ysizemode":"pixel"}],"xaxis2":{"type":"linear","autorange":false,"range":[0.40000000000000002,2.6000000000000001],"tickmode":"array","ticktext":["Port uptake V1","BWT V1"],"tickvals":[1,2],"categoryorder":"array","categoryarray":["Port uptake V1","BWT V1"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"tickangle":-90,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0.53795291466524342,1],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y2","title":"","hoverformat":".2f"},"yaxis2":{"type":"linear","autorange":false,"range":[-550,11550],"tickmode":"array","ticktext":["0","3000","6000","9000"],"tickvals":[0,3000.0000000000005,6000,9000],"categoryorder":"array","categoryarray":["0","3000","6000","9000"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.13283520132835203,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835205},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"domain":[0,1],"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x2","title":"","hoverformat":".2f"},"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.283520132835198},"title":{"text":"Species","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"1d385bfa1ecf":{"main_group":{},"sub_group":{},"x":{},"y":{},"fill":{},"colour":{},"type":"bar"}},"cur_data":"1d385bfa1ecf","visdat":{"1d385bfa1ecf":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>

Save as svg:

``` r
ggsave(filename="R output/BallastSeqR_voyage_species_bar_interactive_fantaxtic_reads.svg", 
       plot=p, width=10, height=6, device=svg)
```

## 7.7 Species level, grouped by voyage, sample type, and tank

### 7.7.1 Create plot with fantaxtic

#### 7.7.1.1 Relative abundance

``` r
#Merge by sample_type and voyage; there is a separate category in the metadata that combines these categories together
Ballast_physeq_merged_group <- speedyseq::merge_samples2(Ballast_physeq, "Sample_group", fun_otu = mean)
SD <- speedyseq::merge_samples2(sample_data(Ballast_physeq), "Sample_group")
print(SD[, "Sample_group"])
```

    ## Sample Data:        [ 9 samples by 1 sample variables ]:
    ##                 Sample_group 
    ##                 <chr>        
    ## 1 V2_UP_T6      V2_UP_T6     
    ## 2 V2_UP_T2      V2_UP_T2     
    ## 3 V2_BWT_T6     V2_BWT_T6    
    ## 4 V1_UP_T5      V1_UP_T5     
    ## 5 V1_UP_T6      V1_UP_T6     
    ## 6 V1_BWT_T5     V1_BWT_T5    
    ## 7 V1_Sea_T5     V1_Sea_T5    
    ## 8 V1_BWT_T6     V1_BWT_T6    
    ## 9 V1_BWT_BWE_T5 V1_BWT_BWE_T5

``` r
#Now samples will be displayed by sample_group rather than individially when graphed using fantaxtic

#Load fantaxtic and create plot
library("fantaxtic")
                        
top_nested_group <- nested_top_taxa(Ballast_physeq_merged_group,
                              #grouping = "Type_voyage",
                              top_tax_level = "Class",
                              nested_tax_level = "Species",
                              n_top_taxa = 3, 
                              n_nested_taxa = 10)

#Order samples
order <- sample_data(top_nested_group$ps_obj) %>%
  data.frame() %>%
  arrange(Sample_group) %>%
  pull(Sample_group) %>%
  as.character()

#Order samples
order <- c("V1_UP_T5", "V1_UP_T6", "V1_BWT_T5", "V1_BWT_T6", "V1_Sea_T5",
           "V1_BWT_BWE_T5", "V2_UP_T2", "V2_UP_T6", "V2_BWT_T6")

#Plot
plot2 <- plot_nested_bar(ps_obj = top_nested_group$ps_obj,
                top_level = "Class",
                nested_level = "Species",
                sample_order = order) +
  facet_nested(~ Voyage + Sample_type, scale="free") +
  labs(x = "Sample group", y = "Relative abundance") +
  theme_minimal()+theme(
    strip.text = element_text(size=12),
    panel.border = element_rect(fill=NA, colour = "black"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(), 
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=10, 
                              colour="black"), 
    axis.text.y = element_text(size=10, colour="black"), 
    plot.title = element_text(hjust = 0.5), 
    axis.ticks.x = element_line(colour="#000000", linewidth=0.1), 
    axis.ticks.y = element_line(colour="#000000", linewidth=0.1),
    legend.text = element_text(size = 10))
```

# 8 Bacterial community composition dot plots

## 8.1 Genus level, grouped by sample type and voyage

``` r
#Calculate relative abundance at the genus level for each sample type and voyage
BallastSeqR_genus_voyage <- Ballast_physeq %>%
  tax_glom(taxrank = "Genus") %>%                     #Agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)}) %>% #Transform to rel. abundance
  psmelt() %>%                                        #Melt to long format
  group_by(Sample_type, Voyage, 
           Kingdom, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarize(Mean = 
                     mean(Abundance, na.rm=TRUE),
                    .groups = "keep") %>%             #Calculate average
  filter(Mean > 0.01) %>%                             #Filter 
  arrange(Class)

#Create color palette
genus_col = colorRampPalette(brewer.pal(9,"Dark2"))(12)
```

    ## Warning in brewer.pal(9, "Dark2"): n too large, allowed maximum for palette Dark2 is 8
    ## Returning the palette you asked for with that many colors

``` r
#Create plot
genus_dot <- ggplot(BallastSeqR_genus_voyage, aes(x=Voyage, y=Genus, 
                                                  color=Class)) +  
  geom_point(aes(size=Mean), alpha = 1/2) +
  scale_color_manual(values = genus_col) +
  facet_nested( ~ Sample_type, scales = "free", 
                nest_line = element_line(colour = "black")) +
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA, colour = "black"), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, 
                                    size=8, colour="black"), 
        axis.text.y = element_text(size=8, colour="black"), 
        plot.title = element_text(hjust = 0.5), 
        axis.ticks.x = element_line(colour="#000000", linewidth = 0.1), 
        axis.ticks.y = element_line(colour="#000000", linewidth = 0.1),
        legend.text = element_text(size = 7),
        ggh4x.facet.nestline = element_line(linetype = "solid"))

print(genus_dot)
```

![](16S-sequence-analysis_files/figure-gfm/genus%20dot%20plot-1.png)<!-- -->

## 8.2 Genus level, grouped by voyage and sample type

Same plot as above, but voyage is grouped first and sample type is
subset under it.

``` r
#Create plot
genus_dot_2 <- ggplot(BallastSeqR_genus_voyage, aes(x=Sample_type, y=Genus, 
                                                  color=Class)) +  
  geom_point(aes(size=Mean), alpha = 1/2) +
  scale_color_manual(values = genus_col) +
  facet_nested( ~ Voyage, scales = "free", 
                nest_line = element_line(colour = "black")) +
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA, colour = "black"), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, 
                                    size=8, colour="black"), 
        axis.text.y = element_text(size=8, colour="black"), 
        plot.title = element_text(hjust = 0.5), 
        axis.ticks.x = element_line(colour="#000000", linewidth = 0.1), 
        axis.ticks.y = element_line(colour="#000000", linewidth = 0.1),
        legend.text = element_text(size = 7),
        ggh4x.facet.nestline = element_line(linetype = "solid"))

print(genus_dot_2)
```

![](16S-sequence-analysis_files/figure-gfm/genus%20dot%20plot%202-1.png)<!-- -->

## 8.3 Species level, grouped by sample type and voyage

Same as the first genus level plot, but displaying the abundance of
species present in abundance \>1%.

``` r
#Calculate relative abundance at the genus level for each sample type and voyage
BallastSeqR_species_voyage <- Ballast_physeq %>%
  tax_glom(taxrank = "Species") %>%                   #Agglomerate at species level
  transform_sample_counts(function(x) {x/sum(x)}) %>% #Transform to rel. abundance
  psmelt() %>%                                        #Melt to long format
  group_by(Sample_type, Voyage, 
           Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::summarize(Mean = 
                     mean(Abundance, na.rm=TRUE),
                    .groups = "keep") %>%             #Calculate average
  filter(Mean > 0.01) %>%                             #Filter out <1% 
  arrange(Class)

#Create color palette
species_col = colorRampPalette(brewer.pal(9,"Dark2"))(12)
```

    ## Warning in brewer.pal(9, "Dark2"): n too large, allowed maximum for palette Dark2 is 8
    ## Returning the palette you asked for with that many colors

``` r
#Create plot
species_dot <- ggplot(BallastSeqR_species_voyage, aes(x=Voyage, y=Species, 
                                                  color=Class)) +  
  geom_point(aes(size=Mean), alpha = 1/2) +
  scale_color_manual(values = species_col) +
  facet_nested( ~ Sample_type, scales = "free", 
                nest_line = element_line(colour = "black")) +
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA, colour = "black"), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, 
                                    size=8, colour="black"), 
        axis.text.y = element_text(size=8, colour="black"), 
        plot.title = element_text(hjust = 0.5), 
        axis.ticks.x = element_line(colour="#000000", linewidth = 0.1), 
        axis.ticks.y = element_line(colour="#000000", linewidth = 0.1),
        legend.text = element_text(size = 7),
        ggh4x.facet.nestline = element_line(linetype = "solid"))

print(species_dot)
```

![](16S-sequence-analysis_files/figure-gfm/species%20dot%20plot-1.png)<!-- -->

## 8.4 Species level, grouped by voyage and sample type

Same as the second genus level plot, but displaying the abundance of
species present in abundance \>1%.

``` r
#Create plot
species_dot_2 <- ggplot(BallastSeqR_species_voyage, aes(x=Sample_type, y=Species, 
                                                  color=Class)) +  
  geom_point(aes(size=Mean), alpha = 1/2) +
  scale_color_manual(values = species_col) +
  facet_nested( ~ Voyage, scales = "free", 
                nest_line = element_line(colour = "black")) +
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA, colour = "black"), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, 
                                    size=8, colour="black"), 
        axis.text.y = element_text(size=8, colour="black"), 
        plot.title = element_text(hjust = 0.5), 
        axis.ticks.x = element_line(colour="#000000", linewidth = 0.1), 
        axis.ticks.y = element_line(colour="#000000", linewidth = 0.1),
        legend.text = element_text(size = 7),
        ggh4x.facet.nestline = element_line(linetype = "solid"))

print(species_dot_2)
```

![](16S-sequence-analysis_files/figure-gfm/species%20dot%20plot%202-1.png)<!-- -->

# 9 Filter ASVs to determine how many were given species-level assignments

Recreate a fresh phyloseq object from the original .qza files so that
the taxonomy is examined here without the relabeling (e.g. `tax_fix`)
applied to `Ballast_physeq` earlier in the analysis.

``` r
#Create a phyloseq object from the .qza files exported from qiime2 using
#the qiime2R package
Ballast_physeq2 <- qza_to_phyloseq(
  features="Qiime output/ballast-dada2-table.qza",
  tree="Qiime output/ballast-rooted-tree.qza",
  taxonomy="Qiime output/ballast-taxonomy.qza",
  metadata = "Qiime output/Metadata.tsv"
)

#Remove chloroplast, mitochondrial and archaeal sequences
Ballast_physeq2 <- subset_taxa(Ballast_physeq2, Kingdom != "d__Archaea")
Ballast_physeq2 <- subset_taxa(Ballast_physeq2, Order != "Chloroplast")
Ballast_physeq2 <- subset_taxa(Ballast_physeq2, Family != "Mitochondria")

#Remove blanks and positive controls
Ballast_physeq2 = subset_samples(Ballast_physeq2, Sample_type != "control")
Ballast_physeq2 = subset_samples(Ballast_physeq2, Sample_type != "negative")

#Remove samples with few reads
Ballast_physeq2 <- subset_samples(Ballast_physeq2,
                                  !(Sample_number %in% samples_to_exclude))

#Filter out low abundance ASVs
Ballast_physeq2 = prune_taxa(taxa_sums(Ballast_physeq2) > 2, Ballast_physeq2)

#Order factors
sample_data(Ballast_physeq2)$Sample_type <- factor(
  sample_data(Ballast_physeq2)$Sample_type,
  levels = sample_type_levels,
  labels = sample_type_labels)

sample_data(Ballast_physeq2)$Sample_number <- factor(
  sample_data(Ballast_physeq2)$Sample_number,
  levels = sample_numbers)

sample_data(Ballast_physeq2)$Voyage <- factor(
  sample_data(Ballast_physeq2)$Voyage,
  levels = c("1", "2"),
  labels = c("Voyage 1", "Voyage 2"))
```

``` r
#Use tax_mutate to edit the species level names
tax_table(Ballast_physeq2) <- tax_table(Ballast_physeq2) %>%
  tax_mutate(Species=gsub("_", " ", Species))

taxtab2 <- as.data.frame(tax_table(Ballast_physeq2))

#Create new dataframe that excludes any NA's in the species column
taxtabfilt <- taxtab2 %>% filter(!is.na(Species))

#Remove unassigned, uncultured and similar
taxtabfilt <- dplyr::filter(taxtabfilt, !grepl('unidentified|uncultured', Species))
taxtabfilt <- dplyr::filter(taxtabfilt, !grepl('metagenome|sp.', Species))

#Remove non-species 
taxtabfilt <- taxtabfilt %>%
  filter(!Species %in% 
           c("Antarctic bacterium", "bacterium enrichment", "microbial mat",
             "Rhodobacteraceae bacterium", "alpha proteobacterium", 
             "Actinobacteria bacterium", "Candidatus Planktophila")) 
```

``` r
#Use tax_mutate to edit the species level names
tax_table(Ballast_physeq2) <- tax_table(Ballast_physeq2) %>%
  tax_mutate(Genus=gsub("_", " ", Genus))

taxtab3 <- as.data.frame(tax_table(Ballast_physeq2))

#Create new dataframe that excludes any NA's in the species column
taxtabfilt_genus <- taxtab3 %>% filter(!is.na(Genus))

#Remove unassigned, uncultured and similar
taxtabfilt_genus <- dplyr::filter(taxtabfilt_genus, 
                                  !grepl('unidentified|uncultured', Genus))
taxtabfilt_genus <- dplyr::filter(taxtabfilt_genus, 
                                  !grepl('clade|group', Genus))

#Remove non-genera 
taxtabfilt_genus <- taxtabfilt_genus %>%
  filter(!Genus %in% 
           c("env.OPS 17", "AKYH767",
             "Bacteroidetes vadinHA17", "Clade IV", "TRA3-20", "R7C24", 
             "JGI 0000069-P22", "DEV007", "UBA12409", " 67-14", "PeM15",
             "IMCC26256", "KD4-96", "S25-593", "T34", "Clade Ia", "Clade II", 
             "67-14", "C39", "Clade III", "GWD2-49-16", "MWH-Ta3", "OLB12", 
             "BSV13", "SH3-11", "SM1A02", "SH-PL14", "HTCC5015", "BIyi10",
             "Clostridium sensu stricto 1", "SUP05 cluster")) 
```

# 10 Session Info

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.2 (2025-10-31 ucrt)
    ##  os       Windows 11 x64 (build 26100)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       America/New_York
    ##  date     2026-07-30
    ##  pandoc   3.6.3 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   1.8.25 @ C:\\PROGRA~1\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package              * version    date (UTC) lib source
    ##    abind                  1.4-8      2024-09-12 [1] CRAN (R 4.5.2)
    ##    ade4                   1.7-24     2026-03-21 [1] CRAN (R 4.5.3)
    ##    ape                  * 5.8-1      2024-12-16 [1] CRAN (R 4.5.3)
    ##    backports              1.5.1      2026-04-03 [1] CRAN (R 4.5.3)
    ##    base64enc              0.1-6      2026-02-02 [1] CRAN (R 4.5.2)
    ##    Biobase              * 2.70.0     2025-10-29 [1] repository (https://github.com/bioc/Biobase@9964e15)
    ##    BiocGenerics         * 0.56.0     2025-10-29 [1] repository (https://github.com/bioc/BiocGenerics@16cf16d)
    ##    BiocManager            1.30.27    2025-11-14 [1] CRAN (R 4.5.3)
    ##    BiocParallel           1.44.0     2025-10-29 [1] repository (https://github.com/bioc/BiocParallel@3d6f2f6)
    ##    biomformat             1.38.3     2026-03-15 [1] repository (https://github.com/bioc/biomformat@0747f37)
    ##    Biostrings             2.78.0     2025-10-29 [1] repository (https://github.com/bioc/Biostrings@eda5d66)
    ##    cachem                 1.1.0      2024-05-16 [1] CRAN (R 4.5.3)
    ##    checkmate              2.3.4      2026-02-03 [1] CRAN (R 4.5.3)
    ##    circlize               0.4.18     2026-04-04 [1] CRAN (R 4.5.3)
    ##    cli                    3.6.6      2026-04-09 [1] CRAN (R 4.5.3)
    ##    clue                   0.3-68     2026-03-26 [1] CRAN (R 4.5.3)
    ##  P cluster              * 2.1.8.1    2025-03-12 [?] CRAN (R 4.5.2)
    ##  P codetools              0.2-20     2024-03-31 [?] CRAN (R 4.5.2)
    ##    colorspace           * 2.1-2      2025-09-22 [1] CRAN (R 4.5.3)
    ##    ComplexHeatmap       * 2.26.1     2026-01-30 [1] repository (https://github.com/bioc/ComplexHeatmap@880562b)
    ##    crayon                 1.5.3      2024-06-20 [1] CRAN (R 4.5.3)
    ##    crosstalk              1.2.2      2025-08-26 [1] CRAN (R 4.5.3)
    ##    data.table             1.18.4     2026-05-06 [1] CRAN (R 4.5.3)
    ##    DelayedArray           0.36.1     2026-03-31 [1] repository (https://github.com/bioc/DelayedArray@1f8a9cb)
    ##    DESeq2               * 1.50.2     2025-11-12 [1] repository (https://github.com/bioc/DESeq2@d90821a)
    ##    devtools               2.5.2      2026-04-30 [1] CRAN (R 4.5.3)
    ##    digest                 0.6.39     2025-11-19 [1] CRAN (R 4.5.3)
    ##    doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.5.3)
    ##    dplyr                * 1.2.1      2026-04-03 [1] CRAN (R 4.5.3)
    ##    DT                     0.34.0     2025-09-02 [1] CRAN (R 4.5.3)
    ##    ellipsis               0.3.3      2026-04-04 [1] CRAN (R 4.5.3)
    ##    evaluate               1.0.5      2025-08-27 [1] CRAN (R 4.5.3)
    ##    fantaxtic            * 0.2.1      2026-06-23 [1] Github (gmteunisse/fantaxtic@b822d7f)
    ##    farver                 2.1.2      2024-05-13 [1] CRAN (R 4.5.3)
    ##    fastmap                1.2.0      2024-05-15 [1] CRAN (R 4.5.3)
    ##    forcats              * 1.0.1      2025-09-25 [1] CRAN (R 4.5.3)
    ##    foreach                1.5.2      2022-02-02 [1] CRAN (R 4.5.3)
    ##  P foreign                0.8-90     2025-03-31 [?] CRAN (R 4.5.2)
    ##    Formula                1.2-5      2023-02-24 [1] CRAN (R 4.5.2)
    ##    fs                     2.1.0      2026-04-18 [1] CRAN (R 4.5.3)
    ##    generics             * 0.1.4      2025-05-09 [1] CRAN (R 4.5.3)
    ##    GenomicRanges        * 1.62.1     2025-12-08 [1] repository (https://github.com/bioc/GenomicRanges@ce11a45)
    ##    GetoptLong             1.1.1      2026-04-08 [1] CRAN (R 4.5.3)
    ##    ggh4x                * 0.3.1      2025-05-30 [1] CRAN (R 4.5.3)
    ##    ggnested             * 0.1.1      2026-06-23 [1] Github (gmteunisse/ggnested@e72a35c)
    ##    ggpattern            * 1.3.2-1    2026-06-23 [1] Github (coolbutuseless/ggpattern@d38fea4)
    ##    ggplot2              * 4.0.3      2026-04-22 [1] CRAN (R 4.5.3)
    ##    ggtext                 0.1.2      2022-09-16 [1] CRAN (R 4.5.3)
    ##    GlobalOptions          0.1.4      2026-04-08 [1] CRAN (R 4.5.3)
    ##    glue                   1.8.1      2026-04-17 [1] CRAN (R 4.5.3)
    ##    gridExtra              2.3        2017-09-09 [1] CRAN (R 4.5.3)
    ##    gridpattern            1.4.2      2026-06-23 [1] CRAN (R 4.5.2)
    ##    gridtext               0.1.6      2026-02-19 [1] CRAN (R 4.5.3)
    ##    gtable                 0.3.6      2024-10-25 [1] CRAN (R 4.5.3)
    ##    Hmisc                  5.2-6      2026-06-19 [1] CRAN (R 4.5.3)
    ##    hms                    1.1.4      2025-10-17 [1] CRAN (R 4.5.3)
    ##    htmlTable              2.5.0      2026-04-22 [1] CRAN (R 4.5.3)
    ##    htmltools              0.5.9      2025-12-04 [1] CRAN (R 4.5.3)
    ##    htmlwidgets            1.6.4      2023-12-06 [1] CRAN (R 4.5.3)
    ##    httr                   1.4.8      2026-02-13 [1] CRAN (R 4.5.3)
    ##    igraph                 2.3.2      2026-05-29 [1] CRAN (R 4.5.3)
    ##    IRanges              * 2.44.0     2025-10-29 [1] repository (https://github.com/bioc/IRanges@964a290)
    ##    iterators              1.0.14     2022-02-05 [1] CRAN (R 4.5.3)
    ##    jsonlite               2.0.0      2025-03-27 [1] CRAN (R 4.5.3)
    ##    kableExtra           * 1.4.1      2026-07-08 [1] CRAN (R 4.5.2)
    ##    knitr                  1.51       2025-12-20 [1] CRAN (R 4.5.3)
    ##    labeling               0.4.3      2023-08-29 [1] CRAN (R 4.5.2)
    ##  P lattice                0.22-7     2025-04-02 [?] CRAN (R 4.5.2)
    ##    lazyeval               0.2.3      2026-04-04 [1] CRAN (R 4.5.3)
    ##    lifecycle              1.0.5      2026-01-08 [1] CRAN (R 4.5.3)
    ##    locfit                 1.5-9.12   2025-03-05 [1] CRAN (R 4.5.3)
    ##    lubridate            * 1.9.5      2026-02-04 [1] CRAN (R 4.5.3)
    ##    magick               * 2.9.1      2026-02-28 [1] CRAN (R 4.5.3)
    ##    magrittr               2.0.5      2026-04-04 [1] CRAN (R 4.5.3)
    ##  P MASS                   7.3-65     2025-02-28 [?] CRAN (R 4.5.2)
    ##  P Matrix                 1.7-4      2025-08-28 [?] CRAN (R 4.5.2)
    ##    MatrixGenerics       * 1.22.0     2025-10-29 [1] repository (https://github.com/bioc/MatrixGenerics@75d9a54)
    ##    matrixStats          * 1.5.0      2025-01-07 [1] CRAN (R 4.5.3)
    ##    memoise                2.0.1      2021-11-26 [1] CRAN (R 4.5.3)
    ##  P mgcv                   1.9-3      2025-04-04 [?] CRAN (R 4.5.2)
    ##    microViz             * 0.13.1     2026-06-23 [1] Github (david-barnett/microViz@0422bf8)
    ##    multtest               2.66.0     2025-10-29 [1] repository (https://github.com/bioc/multtest@2722ca2)
    ##    NADA                   1.6-1.2    2025-08-29 [1] CRAN (R 4.5.3)
    ##  P nlme                 * 3.1-168    2025-03-31 [?] CRAN (R 4.5.2)
    ##  P nnet                   7.3-20     2025-01-01 [?] CRAN (R 4.5.2)
    ##    otel                   0.2.0      2025-08-29 [1] CRAN (R 4.5.3)
    ##    pairwiseAdonis       * 0.4.1      2026-06-23 [1] Github (pmartinezarbizu/pairwiseAdonis@cb190f7)
    ##    patchwork            * 1.3.2      2025-08-25 [1] CRAN (R 4.5.3)
    ##    permute              * 0.9-10     2026-02-06 [1] CRAN (R 4.5.3)
    ##    phyloseq             * 1.54.2     2026-02-27 [1] repository (https://github.com/bioc/phyloseq@2d15ed3)
    ##    phyloseq.extended    * 0.0.0.9000 2026-06-23 [1] Github (mahendra-mariadassou/phyloseq-extended@093aac5)
    ##    picante              * 1.8.2      2020-06-10 [1] CRAN (R 4.5.3)
    ##    pillar                 1.11.1     2025-09-17 [1] CRAN (R 4.5.3)
    ##    pkgbuild               1.4.8      2025-05-26 [1] CRAN (R 4.5.3)
    ##    pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.5.3)
    ##    pkgload                1.5.3      2026-06-15 [1] CRAN (R 4.5.3)
    ##    plotly               * 4.12.0     2026-01-24 [1] CRAN (R 4.5.3)
    ##    plotwidgets            0.5.1      2022-05-10 [1] CRAN (R 4.5.3)
    ##    plyr                 * 1.8.9      2023-10-02 [1] CRAN (R 4.5.3)
    ##    png                    0.1-9      2026-03-15 [1] CRAN (R 4.5.3)
    ##    purrr                * 1.2.2      2026-04-10 [1] CRAN (R 4.5.3)
    ##    qiime2R              * 0.99.6     2026-06-23 [1] Github (jbisanz/qiime2R@fb7a907)
    ##    R6                     2.6.1      2025-02-15 [1] CRAN (R 4.5.3)
    ##    ragg                   1.5.2      2026-03-23 [1] CRAN (R 4.5.3)
    ##    RColorBrewer         * 1.1-3      2022-04-03 [1] CRAN (R 4.5.2)
    ##    Rcpp                   1.1.1-1.1  2026-04-24 [1] CRAN (R 4.5.3)
    ##    readr                * 2.2.0      2026-02-19 [1] CRAN (R 4.5.3)
    ##    renv                   1.2.3      2026-05-16 [1] CRAN (R 4.5.3)
    ##    reshape2             * 1.4.5      2025-11-12 [1] CRAN (R 4.5.3)
    ##    rhdf5                  2.54.1     2025-12-02 [1] repository (https://github.com/bioc/rhdf5@7f691e4)
    ##  D rhdf5filters           1.22.0     2025-10-29 [1] repository (https://github.com/bioc/rhdf5filters@3465c24)
    ##    Rhdf5lib               1.32.0     2025-10-29 [1] repository (https://github.com/bioc/Rhdf5lib@f62ae28)
    ##    rjson                  0.2.23     2024-09-16 [1] CRAN (R 4.5.2)
    ##    rlang                  1.2.0      2026-04-06 [1] CRAN (R 4.5.3)
    ##    rmarkdown              2.31       2026-03-26 [1] CRAN (R 4.5.3)
    ##  P rpart                  4.1.24     2025-01-07 [?] CRAN (R 4.5.2)
    ##    rstudioapi             0.19.0     2026-06-11 [1] CRAN (R 4.5.3)
    ##    S4Arrays               1.10.1     2025-12-01 [1] repository (https://github.com/bioc/S4Arrays@a4cccba)
    ##    S4Vectors            * 0.48.1     2026-04-04 [1] repository (https://github.com/bioc/S4Vectors@ae25d08)
    ##    S7                     0.2.2      2026-04-22 [1] CRAN (R 4.5.3)
    ##    scales               * 1.4.0      2025-04-24 [1] CRAN (R 4.5.3)
    ##    seecolor             * 0.2.0      2023-02-24 [1] CRAN (R 4.5.3)
    ##    Seqinfo              * 1.0.0      2025-10-29 [1] repository (https://github.com/bioc/Seqinfo@9fc5a61)
    ##    sessioninfo            1.2.4      2026-06-04 [1] CRAN (R 4.5.3)
    ##    shape                  1.4.6.1    2024-02-23 [1] CRAN (R 4.5.2)
    ##    SparseArray            1.10.10    2026-03-30 [1] repository (https://github.com/bioc/SparseArray@ae957c5)
    ##    speedyseq            * 0.5.3.9021 2026-06-23 [1] Github (mikemc/speedyseq@0057652)
    ##    stringi                1.8.7      2025-03-27 [1] CRAN (R 4.5.2)
    ##    stringr              * 1.6.0      2025-11-04 [1] CRAN (R 4.5.3)
    ##    SummarizedExperiment * 1.40.0     2025-10-29 [1] repository (https://github.com/bioc/SummarizedExperiment@469a2de)
    ##  P survival               3.8-3      2024-12-17 [?] CRAN (R 4.5.2)
    ##    svglite                2.2.2      2025-10-21 [1] CRAN (R 4.5.3)
    ##    systemfonts            1.3.2      2026-03-05 [1] CRAN (R 4.5.3)
    ##    textshaping            1.0.5      2026-03-06 [1] CRAN (R 4.5.3)
    ##    tibble               * 3.3.1      2026-01-11 [1] CRAN (R 4.5.3)
    ##    tidyr                * 1.3.2      2025-12-19 [1] CRAN (R 4.5.3)
    ##    tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.5.3)
    ##    tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.5.3)
    ##    timechange             0.4.0      2026-01-29 [1] CRAN (R 4.5.3)
    ##    truncnorm              1.0-9      2023-03-20 [1] CRAN (R 4.5.3)
    ##    tzdb                   0.5.0      2025-03-15 [1] CRAN (R 4.5.3)
    ##    usethis                3.2.1      2025-09-06 [1] CRAN (R 4.5.3)
    ##    utf8                   1.2.6      2025-06-08 [1] CRAN (R 4.5.3)
    ##    vctrs                  0.7.3      2026-04-11 [1] CRAN (R 4.5.3)
    ##    vegan                * 2.7-5      2026-05-25 [1] CRAN (R 4.5.3)
    ##    viridis              * 0.6.5      2024-01-29 [1] CRAN (R 4.5.3)
    ##    viridisLite          * 0.4.3      2026-02-04 [1] CRAN (R 4.5.3)
    ##    withr                  3.0.3      2026-06-19 [1] CRAN (R 4.5.3)
    ##    xfun                   0.59       2026-06-19 [1] CRAN (R 4.5.2)
    ##    xml2                   1.6.0      2026-06-22 [1] CRAN (R 4.5.2)
    ##    XVector                0.50.0     2025-10-29 [1] repository (https://github.com/bioc/XVector@6b7e2a1)
    ##    yaml                   2.3.12     2025-12-10 [1] CRAN (R 4.5.3)
    ##    zCompositions          1.6.2      2026-06-23 [1] CRAN (R 4.5.2)
    ## 
    ##  [1] C:/Users/sab236/Box/Home Folder sab236/Private/Ballast-water-bacteria/renv/library/windows/R-4.5/x86_64-w64-mingw32
    ##  [2] C:/Users/sab236/AppData/Local/R/cache/R/renv/sandbox/windows/R-4.5/x86_64-w64-mingw32/ebc25411
    ## 
    ##  * ── Packages attached to the search path.
    ##  P ── Loaded and on-disk path mismatch.
    ##  D ── DLL MD5 mismatch, broken installation.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
