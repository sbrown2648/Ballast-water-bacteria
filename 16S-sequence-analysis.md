Ballast water bacteria 16S rRNA analysis
================
Sarah Brown
June 9, 2023

## Load Packages and Import Data

Here, we load the required packages and then import the .qza files
exported from QIIME2 using the qiime2R package. The final object is a
phyloseq object called Ballast_physeq.

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.2.3

``` r
library(vegan)
```

    ## Warning: package 'vegan' was built under R version 4.2.3

    ## Loading required package: permute

    ## Warning: package 'permute' was built under R version 4.2.3

    ## Loading required package: lattice

    ## Warning: package 'lattice' was built under R version 4.2.3

    ## This is vegan 2.6-4

``` r
library(plyr)
```

    ## Warning: package 'plyr' was built under R version 4.2.3

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.2.3

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
```

    ## Warning: package 'scales' was built under R version 4.2.3

``` r
library(grid)
library(reshape2)
```

    ## Warning: package 'reshape2' was built under R version 4.2.3

``` r
library(phyloseq)
library(picante)
```

    ## Warning: package 'picante' was built under R version 4.2.3

    ## Loading required package: ape

    ## Warning: package 'ape' was built under R version 4.2.3

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

    ## Loading required package: nlme

    ## Warning: package 'nlme' was built under R version 4.2.3

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
library(tidyr)
```

    ## Warning: package 'tidyr' was built under R version 4.2.3

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:reshape2':
    ## 
    ##     smiths

``` r
library(viridis)
```

    ## Warning: package 'viridis' was built under R version 4.2.3

    ## Loading required package: viridisLite

    ## Warning: package 'viridisLite' was built under R version 4.2.3

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:scales':
    ## 
    ##     viridis_pal

``` r
library(qiime2R)
library(DESeq2)
```

    ## Warning: package 'DESeq2' was built under R version 4.2.2

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.2.2

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

    ## Warning: package 'GenomicRanges' was built under R version 4.2.2

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 4.2.2

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## Warning: package 'matrixStats' was built under R version 4.2.3

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

## Preparing the Data

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

### Examine the total number of reads and ASV’s

``` r
readsumsdf = data.frame(nreads = sort(taxa_sums(Ballast_physeq), TRUE), 
                        sorted = 1:ntaxa(Ballast_physeq), type = "ASVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(Ballast_physeq), TRUE), sorted = 1:nsamples(Ballast_physeq), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
```

![](16S-sequence-analysis_files/figure-gfm/asv-1.png)<!-- -->

### Create a rarefaction curve

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
    ##    0.11    0.03    0.14

### Create table showing the number of reads per sample

``` r
#Sum reads and sort by samples with least to greatest number of reads
samples_reads <- sort(sample_sums(Ballast_physeq))

#Create table 
knitr::kable(samples_reads) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

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

### Filter samples

Based on the table above, there are some samples with very few to no
reads, so we’ll filter these out here. In this case, we keep all samples
with more than 1,000 reads.

Removing these samples with low sequencing depths is important because
samples with low reads typically are lower quality and have a greater
probability of containing contaminant sequences ([Weiss et al.,
2017](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y)).

``` r
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
```

### Label and order variables

Now we can relabel the variables; if you type in Ballast_physeq you will
see that sample_data is the matrix that holds the information that we
want to change, so we need to include sample_data in our code here.

``` r
#Order factors
sample_data(Ballast_physeq)$Sample_type <- factor(sample_data(Ballast_physeq)$Sample_type, 
               levels = c("port_uptake", "ocean_uptake", "BWT", "BWT_BWE"),
               labels = c("port_uptake", "ocean_uptake", "BWT", "BWT_BWE"))

sample_data(Ballast_physeq)$Sample_number <- factor(sample_data(Ballast_physeq)$Sample_number, 
               levels = c("80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "100", "101", "102", "103", "110", "111", "112", "113", "114", "115"),
               labels = c("80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "100", "101", "102", "103", "110", "111", "112", "113", "114", "115"))

sample_data(Ballast_physeq)$Voyage <- factor(sample_data(Ballast_physeq)$Voyage, 
               levels = c("1", "2"),
               labels = c("1", "2"))
```

### Abundance transformation

We have to transform the count data to account for differences in
library size between samples. This is often done using rarefaction, but
rarefaction results in the elimination of valid data and can also
increase the rate of false positives when testing for ASV’s that are
diferentially abundant in different sample categories.

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

## Alpha Diversity

Note that while these data are not rarefied, they are normalized. We’ll
use the Shannon diversity index (based on richness AND evenness;
examines how many different taxa are present and how evenly they’re
distributed within a sample) to analyze alpha diversity between variable
types. This means it considers both the number of species and the
inequality between species abundances.

### Calculate Shannon Diversity per Sample

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
#richness <- data.frame(estimate_richness(Ballast_physeqRA, measures = c("Shannon")))
#richness <- setNames(cbind(rownames(richness), richness, row.names = NULL), 
#                     c("sample-id", "Shannon"))

#Add the sample metadata to the dataframe
#s <- data.frame(sample_data(Ballast_physeqRA))
#s <- setNames(cbind(rownames(s), s, row.names = NULL), 
#              c("sample-id", "effluent", "week", "polymer_type", 
#                "bead_diameter", "Channel", "sample_type", "particle_type"))

#alphadiv <- merge(s, richness, by = "sample-id")

#Order factors
#alphadiv$polymer_type <- factor(alphadiv$polymer_type, 
#                        levels = c("Glass", "HDPE", "LDPE", "PP", "PS", "Water"),
#                       labels = c("Glass", "HDPE", "LDPE", "PP", "PS", "Water")) 

#Shows the calculated indices
#knitr::kable(head(alphadiv)) %>% 
  #kableExtra::kable_styling("striped", 
  #                          latex_options="scale_down") %>% 
  #kableExtra::scroll_box(width = "100%")
```

### Create Alpha Diversity Boxplots

#### Comparison: particles vs. water

``` r
#Plot
#ggplot(data=alphadiv, aes(x=sample_type, y=Shannon), alpha=0.1) + 
#  geom_boxplot(aes(fill=sample_type)) +
#  geom_point(position=position_dodge(width=0.75),aes(group=sample_type)) +
#  facet_wrap(~effluent, scale="free") +
#  theme(legend.position="right",
#        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
#        panel.background = element_blank())
```

#### Comparison: particle types (microplastic, glass, or water) between effluent sources

``` r
#Plot
#ggplot(data=alphadiv, aes(x=effluent, y=Shannon), alpha=0.1) + 
#  geom_boxplot(aes(fill=effluent)) +
#  geom_point(position=position_dodge(width=0.75),aes(group=effluent)) +
#  facet_wrap(~particle_type, scale="free") +
#  theme(legend.position="right",
#        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
#        panel.background = element_blank())
```

### Alpha Diversity statistics

## Ordination plots

### Find plot color hex codes

``` r
brewer.pal(4, "Dark2")
```

    ## [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A"

### Create NMDS with the Bray-Curtis dissimilarity matrix

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

    ## Run 0 stress 0.1420338 
    ## Run 1 stress 0.07649612 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1432104  max resid 0.338809 
    ## Run 2 stress 0.07650044 
    ## ... Procrustes: rmse 0.002260751  max resid 0.008134204 
    ## ... Similar to previous best
    ## Run 3 stress 0.07649569 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000242526  max resid 0.0005346942 
    ## ... Similar to previous best
    ## Run 4 stress 0.07649482 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001099109  max resid 0.004510155 
    ## ... Similar to previous best
    ## Run 5 stress 0.07650056 
    ## ... Procrustes: rmse 0.002175056  max resid 0.007691398 
    ## ... Similar to previous best
    ## Run 6 stress 0.07650031 
    ## ... Procrustes: rmse 0.001922759  max resid 0.007579412 
    ## ... Similar to previous best
    ## Run 7 stress 0.07649621 
    ## ... Procrustes: rmse 0.0009348649  max resid 0.004271182 
    ## ... Similar to previous best
    ## Run 8 stress 0.07650048 
    ## ... Procrustes: rmse 0.002154482  max resid 0.007662833 
    ## ... Similar to previous best
    ## Run 9 stress 0.07650136 
    ## ... Procrustes: rmse 0.002338164  max resid 0.007678458 
    ## ... Similar to previous best
    ## Run 10 stress 0.07650086 
    ## ... Procrustes: rmse 0.001916288  max resid 0.007723962 
    ## ... Similar to previous best
    ## Run 11 stress 0.07649574 
    ## ... Procrustes: rmse 0.001118349  max resid 0.004579129 
    ## ... Similar to previous best
    ## Run 12 stress 0.07649592 
    ## ... Procrustes: rmse 0.001179995  max resid 0.004786112 
    ## ... Similar to previous best
    ## Run 13 stress 0.07649556 
    ## ... Procrustes: rmse 0.0009011258  max resid 0.003437049 
    ## ... Similar to previous best
    ## Run 14 stress 0.07650051 
    ## ... Procrustes: rmse 0.002165229  max resid 0.007677577 
    ## ... Similar to previous best
    ## Run 15 stress 0.129975 
    ## Run 16 stress 0.0764961 
    ## ... Procrustes: rmse 0.0009225967  max resid 0.004179171 
    ## ... Similar to previous best
    ## Run 17 stress 0.1299779 
    ## Run 18 stress 0.129975 
    ## Run 19 stress 0.1339254 
    ## Run 20 stress 0.1483185 
    ## *** Best solution repeated 12 times

``` r
#Plot, color coding by Sample_type
all.nmds.sample.type <- plot_ordination(
  physeq = Ballast_physeqRA,
  ordination = all.nmds.source.ord) + 
  scale_colour_manual(values = c("port_uptake" = "#E7298A", "ocean_uptake" = "#7570B3", "BWT" = "#1B9E77", "BWT_BWE" = "#D95F02"), "Sample Type") +
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

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(all.nmds.sample.type)
```

![](16S-sequence-analysis_files/figure-gfm/NMDS%20Bray-1.png)<!-- -->

### Create PCoA with the Bray-Curtis dissimilarity matrix

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
  scale_colour_manual(values = c("port_uptake" = "#E7298A", "ocean_uptake" = "#7570B3", "BWT" = "#1B9E77", "BWT_BWE" = "#D95F02"), "Sample Type") +
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
        legend.justification = c("right", "top")) 

print(all.pcoa.sample.type)
```

![](16S-sequence-analysis_files/figure-gfm/PCoA%20Bray-1.png)<!-- -->

### 

## Preliminary taxonomy fix

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

## Bacterial community composition bar charts

### Family level, grouped by sample type

## Session Info

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
    ##  date     2023-06-12
    ##  pandoc   2.18 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package              * version    date (UTC) lib source
    ##    ade4                   1.7-22     2023-02-06 [1] CRAN (R 4.2.3)
    ##    annotate               1.76.0     2022-11-03 [1] Bioconductor
    ##    AnnotationDbi          1.60.2     2023-03-10 [1] Bioconductor
    ##    ape                  * 5.7-1      2023-03-13 [1] CRAN (R 4.2.3)
    ##    backports              1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
    ##    base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
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
    ##    codetools              0.2-19     2023-02-01 [1] CRAN (R 4.2.2)
    ##    colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.2.3)
    ##    ComplexHeatmap       * 2.14.0     2022-11-01 [1] Bioconductor
    ##    crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.2.3)
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
    ##    evaluate               0.21       2023-05-05 [1] CRAN (R 4.2.3)
    ##    fansi                  1.0.4      2023-01-22 [1] CRAN (R 4.2.3)
    ##    farver                 2.1.1      2022-07-06 [1] CRAN (R 4.2.3)
    ##    fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.2.3)
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
    ##    gtable                 0.3.3      2023-03-21 [1] CRAN (R 4.2.3)
    ##    highr                  0.10       2022-12-22 [1] CRAN (R 4.2.3)
    ##    Hmisc                  5.1-0      2023-05-08 [1] CRAN (R 4.2.3)
    ##    htmlTable              2.4.1      2022-07-07 [1] CRAN (R 4.2.3)
    ##    htmltools              0.5.5      2023-03-23 [1] CRAN (R 4.2.3)
    ##    htmlwidgets            1.6.2      2023-03-17 [1] CRAN (R 4.2.3)
    ##    httpuv                 1.6.11     2023-05-11 [1] CRAN (R 4.2.3)
    ##    httr                   1.4.6      2023-05-08 [1] CRAN (R 4.2.3)
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
    ##    lifecycle              1.0.3      2022-10-07 [1] CRAN (R 4.2.3)
    ##    locfit                 1.5-9.7    2023-01-02 [1] CRAN (R 4.2.3)
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
    ##    NADA                   1.6-1.1    2020-03-22 [1] CRAN (R 4.2.3)
    ##    nlme                 * 3.1-162    2023-01-31 [1] CRAN (R 4.2.3)
    ##    nnet                   7.3-19     2023-05-03 [1] CRAN (R 4.2.3)
    ##    patchwork            * 1.1.2.9000 2023-04-24 [1] Github (thomasp85/patchwork@c14c960)
    ##    permute              * 0.9-7      2022-01-27 [1] CRAN (R 4.2.3)
    ##    phyloseq             * 1.42.0     2022-11-01 [1] Bioconductor
    ##    picante              * 1.8.2      2020-06-10 [1] CRAN (R 4.2.3)
    ##    pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.2.3)
    ##    pkgbuild               1.4.0      2022-11-27 [1] CRAN (R 4.2.3)
    ##    pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.2.3)
    ##    pkgload                1.3.2      2022-11-16 [1] CRAN (R 4.2.3)
    ##    plyr                 * 1.8.8      2022-11-11 [1] CRAN (R 4.2.3)
    ##    png                    0.1-8      2022-11-29 [1] CRAN (R 4.2.2)
    ##    prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.2.3)
    ##    processx               3.8.1      2023-04-18 [1] CRAN (R 4.2.1)
    ##    profvis                0.3.8      2023-05-02 [1] CRAN (R 4.2.3)
    ##    promises               1.2.0.1    2021-02-11 [1] CRAN (R 4.2.3)
    ##    ps                     1.7.5      2023-04-18 [1] CRAN (R 4.2.3)
    ##    purrr                  1.0.1      2023-01-10 [1] CRAN (R 4.2.3)
    ##    qiime2R              * 0.99.6     2023-04-07 [1] Github (jbisanz/qiime2R@2a3cee1)
    ##    R6                     2.5.1      2021-08-19 [1] CRAN (R 4.2.3)
    ##    RColorBrewer         * 1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
    ##    Rcpp                   1.0.10     2023-01-22 [1] CRAN (R 4.2.3)
    ##    RCurl                  1.98-1.12  2023-03-27 [1] CRAN (R 4.2.3)
    ##    remotes                2.4.2      2021-11-30 [1] CRAN (R 4.2.3)
    ##    reshape2             * 1.4.4      2020-04-09 [1] CRAN (R 4.2.3)
    ##    rhdf5                  2.42.0     2022-11-01 [1] Bioconductor
    ##  D rhdf5filters           1.10.1     2023-03-24 [1] Bioconductor
    ##    Rhdf5lib               1.20.0     2022-11-01 [1] Bioconductor
    ##    rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.2.0)
    ##    rlang                  1.1.1      2023-04-28 [1] CRAN (R 4.2.3)
    ##    rmarkdown              2.22       2023-06-01 [1] CRAN (R 4.2.3)
    ##    rpart                  4.1.19     2022-10-21 [1] CRAN (R 4.2.3)
    ##    RSQLite                2.3.1      2023-04-03 [1] CRAN (R 4.2.3)
    ##    rstudioapi             0.14       2022-08-22 [1] CRAN (R 4.2.3)
    ##    rvest                  1.0.3      2022-08-19 [1] CRAN (R 4.2.3)
    ##    S4Vectors            * 0.36.2     2023-02-26 [1] Bioconductor
    ##    scales               * 1.2.1      2022-08-20 [1] CRAN (R 4.2.3)
    ##    sessioninfo            1.2.2      2021-12-06 [1] CRAN (R 4.2.3)
    ##    shape                  1.4.6      2021-05-19 [1] CRAN (R 4.2.0)
    ##    shiny                  1.7.4      2022-12-15 [1] CRAN (R 4.2.3)
    ##    speedyseq            * 0.5.3.9018 2023-05-01 [1] Github (mikemc/speedyseq@ceb941f)
    ##    stringi                1.7.12     2023-01-11 [1] CRAN (R 4.2.2)
    ##    stringr                1.5.0      2022-12-02 [1] CRAN (R 4.2.3)
    ##    SummarizedExperiment * 1.28.0     2022-11-01 [1] Bioconductor
    ##    survival               3.5-5      2023-03-12 [1] CRAN (R 4.2.3)
    ##    svglite                2.1.1      2023-01-10 [1] CRAN (R 4.2.3)
    ##    systemfonts            1.0.4      2022-02-11 [1] CRAN (R 4.2.3)
    ##    tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.2.3)
    ##    tidyr                * 1.3.0      2023-01-24 [1] CRAN (R 4.2.3)
    ##    tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.2.3)
    ##    truncnorm              1.0-9      2023-03-20 [1] CRAN (R 4.2.3)
    ##    urlchecker             1.0.1      2021-11-30 [1] CRAN (R 4.2.3)
    ##    usethis                2.2.0      2023-06-06 [1] CRAN (R 4.2.3)
    ##    utf8                   1.2.3      2023-01-31 [1] CRAN (R 4.2.3)
    ##    vctrs                  0.6.2      2023-04-19 [1] CRAN (R 4.2.3)
    ##    vegan                * 2.6-4      2022-10-11 [1] CRAN (R 4.2.3)
    ##    viridis              * 0.6.3      2023-05-03 [1] CRAN (R 4.2.3)
    ##    viridisLite          * 0.4.2      2023-05-02 [1] CRAN (R 4.2.3)
    ##    webshot                0.5.4      2022-09-26 [1] CRAN (R 4.2.3)
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
