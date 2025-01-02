Preprocess RNA-seq
================
Kevin Thomas
6/20/2023

- [Preprocessing cell type-specific RNA-seq
  data](#preprocessing-cell-type-specific-rna-seq-data)
  - [Setup](#setup)
  - [Transcript import](#transcript-import)
  - [Filtering out small rRNA reads](#filtering-out-small-rrna-reads)
  - [Metadata import](#metadata-import)
  - [Make DGEList](#make-dgelist)
    - [Identify and remove outliers](#identify-and-remove-outliers)

# Preprocessing cell type-specific RNA-seq data

## Setup

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, table,
    ##     tapply, union, unique, unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

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

``` r
library(tximport)
library(tximeta)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::count()        masks matrixStats::count()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(openxlsx)
library(edgeR)
```

    ## Loading required package: limma
    ## 
    ## Attaching package: 'limma'
    ## 
    ## The following object is masked from 'package:DESeq2':
    ## 
    ##     plotMA
    ## 
    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
# Parallel backend
BPPARAM =
  BiocParallel::SnowParam(
    workers       = parallel::detectCores()-2,
    exportglobals = FALSE,
    progressbar   = TRUE,
    type = "SOCK"
  )
BiocParallel::register(BPPARAM)
```

## Transcript import

After pseudoalignment with salmon, transcript counts are imported into R
via `tximport`.

``` r
# Get quant files from output of files
data_dir <- "~/workspace/datasets/rnaseq/lra-bms-advantaseq/data"
count_files <- dir(
  data_dir, 
  pattern = "*.sf.gz",
  full.names = TRUE,
  recursive = TRUE
)
count_files <- count_files[-length(count_files)] ## Remove the last quant file, which contains the undetermined barcodes

# Annotations for tximport
annotations <- read_csv(
  file = "../references/gencode_v32_virus_tx2gene_v1.2.csv"
)
```

    ## Rows: 228646 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): transcript, gene_name
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Transcript import
tx_counts <- tximport(
  files = count_files,
  type = "salmon",
  txIn = TRUE,
  txOut = FALSE,
  tx2gene = annotations,
  importer = data.table::fread
)
```

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 
    ## removing duplicated transcript rows from tx2gene
    ## transcripts missing from tx2gene: 37
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

## Filtering out small rRNA reads

``` r
# Filter out small rRNA reads (RNA5 prefix)
filtered_counts <- purrr::map(
  .x = c("abundance", "counts", "length"),
  .f = \(x){
    cleaned_counts <-
      tx_counts[[x]] |>
      tibble::as_tibble(rownames = "gene_symbol") |>
      dplyr::filter(
        stringr::str_detect(
          string = gene_symbol,
          pattern = "^RNA5",
          negate = TRUE
        )
      )
    
    cleaned_counts <-
      cleaned_counts |>
      dplyr::mutate(across(-matches("gene_symbol"), as.numeric)) |>
      tibble::column_to_rownames("gene_symbol") |>
      as.matrix()
  }) |> 
  set_names(c("abundance", "counts", "length"))
```

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
    ## `.name_repair` is omitted as of tibble 2.0.0.
    ## ℹ Using compatibility `.name_repair`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
filtered_counts[["countsFromAbundance"]] <- "no"
```

## Metadata import

``` r
advantaseq_md <- readRDS("../data/advantaseq_md.RDS")
```

## Make DGEList

``` r
# Make DeSeq Dataset object
dds_import <- DESeqDataSetFromTximport(
  txi = filtered_counts,
  colData = advantaseq_md,
  design = ~ test_group
)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

Datasets must have at leaset 5 million counts for further use.

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](00_preprocess_rna_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
\### Filter out libraries and genes with low counts

``` r
# Remove libraries with low counts
low_count_idx <- counts(dds_import) %>%
  colSums() %>%
  (function(x) which(x < 5e6))
dds_import <- dds_import[,-low_count_idx]

# Count matrix, keep transcripts with at least 10 counts represented in a group of interest
ct_mtx <- counts(dds_import)
keep <- which(rowSums(counts(dds_import) >= 10) >= min(table(dds_import$test_group)))
ct_mtx <- ct_mtx[keep,]

# DGEList
d0 <- DGEList(ct_mtx, samples = as.data.frame(colData(dds_import)))

# Normalize and filter out transcripts with CPM < 1
d0 <- calcNormFactors(d0)
drop <- which(apply(cpm(d0), 1, max) < 1)
d <- d0[-drop,] 
dim(d) # number of genes and samples left
```

    ## [1] 19274   128

### Identify and remove outliers

``` r
d <- calcNormFactors(d)
# Perform principal component analysis to identify outliers
cpm(d)%>%
  t() %>%
  prcomp(center = TRUE, scale. = TRUE) %>%
  .$x %>%
  bind_cols(d$samples) %>%
  ggplot(aes(x = PC1, y = PC2, col = SampleType)) +
  geom_text(aes(label = SampleName)) +
  theme_bw()
```

![](00_preprocess_rna_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Remove outlier AD055-22
d <- d[,!d$samples$SampleName %in% c("AD055-22")]
d <- calcNormFactors(d)
```

``` r
# Save
saveRDS(object = d, file = "../data/rna_dgelist.RDS")
```

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] edgeR_4.2.1                 limma_3.60.4               
    ##  [3] openxlsx_4.2.6.1            lubridate_1.9.3            
    ##  [5] forcats_1.0.0               stringr_1.5.1              
    ##  [7] dplyr_1.1.4                 purrr_1.0.2                
    ##  [9] readr_2.1.5                 tidyr_1.3.1                
    ## [11] tibble_3.2.1                ggplot2_3.5.1              
    ## [13] tidyverse_2.0.0             tximeta_1.22.1             
    ## [15] tximport_1.32.0             DESeq2_1.44.0              
    ## [17] SummarizedExperiment_1.34.0 Biobase_2.64.0             
    ## [19] MatrixGenerics_1.16.0       matrixStats_1.3.0          
    ## [21] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
    ## [23] IRanges_2.38.1              S4Vectors_0.42.1           
    ## [25] BiocGenerics_0.50.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rstudioapi_0.16.0        jsonlite_1.8.8           magrittr_2.0.3          
    ##   [4] GenomicFeatures_1.56.0   farver_2.1.2             rmarkdown_2.27          
    ##   [7] BiocIO_1.14.0            zlibbioc_1.50.0          vctrs_0.6.5             
    ##  [10] memoise_2.0.1            Rsamtools_2.20.0         RCurl_1.98-1.16         
    ##  [13] htmltools_0.5.8.1        S4Arrays_1.4.1           progress_1.2.3          
    ##  [16] AnnotationHub_3.12.0     curl_5.2.1               SparseArray_1.4.8       
    ##  [19] httr2_1.0.2              cachem_1.1.0             GenomicAlignments_1.40.0
    ##  [22] lifecycle_1.0.4          pkgconfig_2.0.3          Matrix_1.7-0            
    ##  [25] R6_2.5.1                 fastmap_1.2.0            GenomeInfoDbData_1.2.12 
    ##  [28] digest_0.6.36            colorspace_2.1-1         AnnotationDbi_1.66.0    
    ##  [31] RSQLite_2.3.7            labeling_0.4.3           filelock_1.0.3          
    ##  [34] fansi_1.0.6              timechange_0.3.0         httr_1.4.7              
    ##  [37] abind_1.4-5              compiler_4.4.1           bit64_4.0.5             
    ##  [40] withr_3.0.1              BiocParallel_1.38.0      DBI_1.2.3               
    ##  [43] highr_0.11               R.utils_2.12.3           biomaRt_2.60.1          
    ##  [46] rappdirs_0.3.3           DelayedArray_0.30.1      rjson_0.2.21            
    ##  [49] tools_4.4.1              zip_2.3.1                R.oo_1.26.0             
    ##  [52] glue_1.7.0               restfulr_0.0.15          grid_4.4.1              
    ##  [55] generics_0.1.3           gtable_0.3.5             tzdb_0.4.0              
    ##  [58] R.methodsS3_1.8.2        ensembldb_2.28.0         data.table_1.15.4       
    ##  [61] hms_1.1.3                xml2_1.3.6               utf8_1.2.4              
    ##  [64] XVector_0.44.0           BiocVersion_3.19.1       pillar_1.9.0            
    ##  [67] vroom_1.6.5              BiocFileCache_2.12.0     lattice_0.22-6          
    ##  [70] rtracklayer_1.64.0       bit_4.0.5                tidyselect_1.2.1        
    ##  [73] locfit_1.5-9.10          Biostrings_2.72.1        knitr_1.48              
    ##  [76] ProtGenerics_1.36.0      xfun_0.46                statmod_1.5.0           
    ##  [79] stringi_1.8.4            UCSC.utils_1.0.0         lazyeval_0.2.2          
    ##  [82] yaml_2.3.10              evaluate_0.24.0          codetools_0.2-20        
    ##  [85] BiocManager_1.30.23      cli_3.6.3                munsell_0.5.1           
    ##  [88] Rcpp_1.0.13              dbplyr_2.5.0             png_0.1-8               
    ##  [91] XML_3.99-0.17            parallel_4.4.1           blob_1.2.4              
    ##  [94] prettyunits_1.2.0        AnnotationFilter_1.28.0  bitops_1.0-8            
    ##  [97] txdbmaker_1.0.1          scales_1.3.0             crayon_1.5.3            
    ## [100] rlang_1.1.4              KEGGREST_1.44.1
