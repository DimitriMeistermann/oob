oob: a toolbox for analyzing OMIC data Out Of Bounds
================

This package provide function for analyzing OMIC data, with a focus on
plots, functional enrichment and bulk/single-cell RNA-Seq.

## Installation

In a R console:

    install.packages("devtools")
    devtools::install_github("https://github.com/DimitriMeistermann/oob", build_vignettes = FALSE)

For a manual installation of dependencies:

    install.packages("BiocManager")
    BiocManager::install(
      c("ggplot2","basilisk","batchelor","BiocParallel","BiocGenerics","circlize","ComplexHeatmap","data.table","dbscan",
      "ggbeeswarm","ggrepel","glmnet","grDevices","grid","irlba","labeling","lsa","MASS","Matrix","pvclust","reshape2",
      "reticulate","rgl","Rcpp","RcppArmadillo","rlang","scales","scater","scattermore","scran","SingleCellExperiment",
      "SummarizedExperiment","stringr","cli","uwot","graphics","methods","utils","knitr","rmarkdown","qualpalr","WGCNA")
    )

*oob* needs also several python dependencies accessible in the
environment used by *reticulate*, the R package for calling python
functions. See
<https://rstudio.github.io/reticulate/articles/python_packages.html> for
more information on installation of python package via reticulate.

For a quick installation:

    reticulate::py_install(c("numpy","leidenalg","trimap"))

The package is ready to load.

``` r
  library(oob)
```

## Data formatting

Fast read and write of text files containing a dataframe or matrix with
`fastRead` and `fastWrite`. Default format are tabulated separated
value. `chead` is then useful for a quick view of first columns/rows.
`rn`, `cn`, `len`, `inter` are respectively aliases for `rownames`,
`colnames`, `length` and `intersect`.

``` r
data("bulkLogCounts") #Bulk RNA-Seq log counts from from Kilens, Meistermann & al 2018

dir.create("test",showWarnings = FALSE)

fastWrite(bulkLogCounts,"test/bulkLogCounts.tsv")
rm(bulkLogCounts)
bulkLogCounts<-fastRead("test/bulkLogCounts.tsv")
chead(bulkLogCounts, n = 10) #quick view of bulkLogCounts
#>        P3E07LSP32 R3E08MIP01 A3E09MIP20 A3E10MIP25 E3E11LEJ07 R3E12BJP01 R3F01BJP03 R3F02MIP03 K3F03BJP03 K3F04MIP04
#> A1BG       0.7031     0.7982     0.0442     0.0442     2.6915     1.5739     0.7296     1.3773     0.9448     1.3465
#> A2M        0.0256     0.9958     0.0000     0.0000     0.2450     1.2847     0.1163     0.0157     0.9306     0.8830
#> A2ML1      1.0645     1.1010     0.3719     1.2047     0.3760     1.6005     1.3253     0.3663     1.2232     0.3834
#> A4GALT     1.1246     1.5890     2.8936     2.3702     1.2442     1.7533     1.6792     1.5236     0.4104     0.4421
#> AAAS       3.3294     3.0880     2.6811     3.0237     2.3443     2.8623     3.0138     3.1165     3.6051     3.6174
#> AACS       2.0621     1.3823     1.0639     1.6621     1.0525     1.8674     2.0116     1.2994     1.4468     0.0000
#> AADAT      0.1166     1.2819     1.0213     0.1895     2.1091     1.7016     1.4465     1.2142     1.8129     1.1494
#> AAED1      1.2258     1.0843     0.0000     0.0000     1.6439     0.7909     0.7430     0.0181     0.6644     0.9484
#> AAGAB      4.0413     3.7013     3.5889     3.6993     3.8790     4.0598     3.7633     3.3143     3.8737     3.9129
#> AAK1       1.1232     1.4225     1.5251     0.0622     1.6696     1.5750     1.6063     0.0825     0.8750     0.8958

rn(bulkLogCounts)[1:10]
#>  [1] "A1BG"   "A2M"    "A2ML1"  "A4GALT" "AAAS"   "AACS"   "AADAT"  "AAED1"  "AAGAB"  "AAK1"
cn(bulkLogCounts)[1:10]
#>  [1] "P3E07LSP32" "R3E08MIP01" "A3E09MIP20" "A3E10MIP25" "E3E11LEJ07" "R3E12BJP01" "R3F01BJP03" "R3F02MIP03" "K3F03BJP03" "K3F04MIP04"
len(rn(bulkLogCounts))
#> [1] 16959
inter(rn(bulkLogCounts)[1:10],rn(bulkLogCounts)[5:15])
#> [1] "AAAS"  "AACS"  "AADAT" "AAED1" "AAGAB" "AAK1"
```

List of vectors can be processed to factor via `VectorListToFactor` or
reciprocally, `factorToVectorList` converts A factor to a list of
vectors, where each level is an element and contains the observation
names having this level. `read.vectorList` and `write.vectorList` are
used to read or write list of vectors in text files.

``` r
data("DEgenesPrime_Naive") #Differential expression (prime vs naive cell lines) from Kilens, Meistermann & al 2018
genesPerSet<-factorToVectorList(DEgenesPrime_Naive$isDE,factorNames = rn(DEgenesPrime_Naive))
lapply(genesPerSet,head)
#> $DOWNREG
#> [1] "A1BG"    "AADAT"   "AASS"    "ABAT"    "ABHD17B" "ABI2"   
#> 
#> $NONE
#> [1] "A2M"   "A2ML1" "AAAS"  "AACS"  "AAED1" "AAGAB"
#> 
#> $UPREG
#> [1] "A4GALT" "AARS2"  "ABCA1"  "ABCC13" "ABCC2"  "ABCF2"

write.vectorList(genesPerSet,"test/genesDE.tsv")
read.vectorList("test/genesDE.tsv") |> lapply(head)
#> $DOWNREG
#> [1] "A1BG"    "AADAT"   "AASS"    "ABAT"    "ABHD17B" "ABI2"   
#> 
#> $NONE
#> [1] "A2M"   "A2ML1" "AAAS"  "AACS"  "AAED1" "AAGAB"
#> 
#> $UPREG
#> [1] "A4GALT" "AARS2"  "ABCA1"  "ABCC13" "ABCC2"  "ABCF2"

top10FC<-rn(DEgenesPrime_Naive)[ whichTop(DEgenesPrime_Naive$log2FoldChange,top = 10) ] #whichTop used here for retrieving top 10 genes in term of log2(Fold-Change)
copyReadyVector(top10FC) #Ready to be copied in a console !
#> [1] "c('KHDC1L','DNMT3L','NLRP7','SPIC','OLAH','MAGEB2','SUN3','FAM151A','DPPA5','SYCP3')"
```

When set `set.seed` is used the Random Number Generator give a different
result each time it is used. `getRandState`, `setRandState` are used for
retrieving/setting this precise seed state.

``` r
set.seed(666)
rnorm(10)
#>  [1]  0.75331105  2.01435467 -0.35513446  2.02816784 -2.21687445  0.75839618 -1.30618526 -0.80251957 -1.79224083 -0.04203245
rnorm(10)
#>  [1]  2.15004262 -1.77023084  0.86465359 -1.72015590  0.13412567 -0.07582656  0.85830054  0.34490035 -0.58245269  0.78617038
set.seed(666)
rnorm(10)
#>  [1]  0.75331105  2.01435467 -0.35513446  2.02816784 -2.21687445  0.75839618 -1.30618526 -0.80251957 -1.79224083 -0.04203245
randState<-getRandState()
rnorm(10)
#>  [1]  2.15004262 -1.77023084  0.86465359 -1.72015590  0.13412567 -0.07582656  0.85830054  0.34490035 -0.58245269  0.78617038
setRandState(randState)
rnorm(10)
#>  [1]  2.15004262 -1.77023084  0.86465359 -1.72015590  0.13412567 -0.07582656  0.85830054  0.34490035 -0.58245269  0.78617038
```

## Statistics utilities

The following functions can be applied to a vector of numeric values:

- `Mode`: (the mode of a distribution)

- `gmean`: geometrical mean $\sqrt[n]{\prod^n_ix_i}$ - `cv`: coefficient
  of variation $\sigma/\mu$ - `cv2`: coefficient of variation (squared
  standard deviation and mean) $\sigma^2/\mu^2$ - `Mode` (most
  represented value in distribution) - `se`: Standard mean error
  $\sigma/\sqrt{n}$

- Distance of correlation $2-cor(x)$ for returning a distance object
  usable by function of clustering as `hclust`, two can be used:

- `corrDist`: Method used for computing correlation can be specified.
  See `method` argument from `cor`.

- `linearScale` can returned a function able to map one set of value to
  another

``` r
minutesTofreqScale<-linearScale(c(0,60),c(0,1))
minutesTofreqScale(30)
#> [1] 0.5
```

`reScale` extend the concept of scaling ranges to matrices. A typical
use is to put the range of a batch corrected count table from RNA-Seq on
the same range than the count table before correction.

``` r
library(ggplot2)
m1<-matrix(rnorm(100),ncol=50)
p1<-qplotDensity(m1[1,],returnGraph = TRUE)+ggtitle("Distribution of matrix to adjust")
m2<-matrix(rnorm(100,20),ncol=50)
p2<-qplotDensity(m2[1,],returnGraph = TRUE)+ggtitle("Distribution of matrix with target ranges")
m1Rescale<-reScale(m1,m2)
p3<-qplotDensity(m1Rescale[1,],returnGraph = TRUE)+ggtitle("Readjusted matrix")

multiplot(p1,p2,p3,cols = 2) #multiplot for quick display of several plot at once
```

<img src="ReadMe_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

`qplotDensity` used here is a quick plotting function for visualizing
distribution. *oob* contains several quick plot wrapper from
*ggplot2*: - `qplotAutoX`: Quick plot of a numeric vector with index as
x-axis. - `qplotBarplot`: Same but with bars instead of points.

`uncenter` is similar to `reScale` but just shift the values of each
feature so the new minimum is 0. `rowScale` scales a matrix in the same
way than `scale`, but with rows as features, the common format for
RNA-Seq.

`aggregMatPerVector` is

``` r
data("sampleAnnot")
aggregMatPerVector(bulkLogCounts[c(
  'KHDC1L','DNMT3L','NLRP7','SPIC','OLAH','MAGEB2'
),], sampleAnnot$culture_media)
#>               E7    Fibro   KSR+FGF2  mTeSR1     RSeT   T2iLGO
#> KHDC1L 2.2209000 1.484425 1.49697778 1.60990 5.553344 8.647927
#> DNMT3L 2.5441500 1.417825 1.86412222 1.64406 5.764344 8.024754
#> NLRP7  1.6627375 0.956650 0.83468889 0.97890 4.465220 6.317400
#> SPIC   0.2433500 0.044725 0.08795556 0.36446 1.471476 3.400135
#> OLAH   1.3360875 0.248050 0.31224444 0.75298 2.090316 3.762662
#> MAGEB2 0.5168875 0.165325 0.26643333 0.62160 2.083200 3.695100
```

## Normalization of RNA-Seq

The package can perform several quick normalization methods. Let’s
simulate a fake count table and attribute real gene symbols as row
names.

``` r
data("geneLengthGRCh38") # vector of gene length for GRCh38 Human reference, named by gene symbols
library(MASS)
countMat<-t(sapply(vector("numeric",length = length(geneLengthGRCh38)),function(x){ #generate a fake RNA-Seq dataset
  rnegbin(10,theta = abs(rnorm(1,mean = 10,sd = 20)),mu = abs(rnorm(1,mean = 10,sd = 20)))
}))
colnames(countMat)<-letters[1:ncol(countMat)];rownames(countMat)<-names(geneLengthGRCh38)
chead(countMat)
#>           a  b  c  d  e
#> TSPAN6   20  9 23 22 10
#> TNMD     13 57 11 27 13
#> DPM1      0  0  0  0  3
#> SCYL3    25 10  6 11 10
#> C1orf112 14 30 15 11 15

#Quality control
computeQCmetricSamples(countMat)
#>       mean       sd        CV TotalGenEx TotalCount
#> a 17.84183 18.53505 1.0388539      23804     446349
#> b 17.83439 17.89979 1.0036672      23775     446163
#> c 17.82716 17.89028 1.0035406      23764     445982
#> d 17.85502 17.54812 0.9828118      23800     446679
#> e 17.88184 18.96104 1.0603517      23775     447350
#> f 17.81964 20.50809 1.1508696      23773     445794
#> g 17.74733 17.75208 1.0002677      23794     443985
#> h 17.73734 16.95955 0.9561496      23776     443735
#> i 17.78247 17.15649 0.9647978      23818     444864
#> j 17.80625 16.79565 0.9432443      23806     445459

#Different kind of normalization
CPM(countMat) |> chead() #Count Per Million normalization
#>                 a         b        c        d         e
#> TSPAN6   44.80799  20.17200 51.57159 49.25237 22.353862
#> TNMD     29.12519 127.75600 24.66467 60.44609 29.060020
#> DPM1      0.00000   0.00000  0.00000  0.00000  6.706158
#> SCYL3    56.00998  22.41333 13.45346 24.62619 22.353862
#> C1orf112 31.36559  67.24000 33.63364 24.62619 33.530792
normDeseq(countMat) |> chead() #DESeq2 normalization
#>                 a         b         c        d         e
#> TSPAN6   19.40455  8.739744 22.259615 21.26538  9.684225
#> TNMD     12.61296 55.351712 10.645903 26.09842 12.589492
#> DPM1      0.00000  0.000000  0.000000  0.00000  2.905267
#> SCYL3    24.25569  9.710827  5.806856 10.63269  9.684225
#> C1orf112 13.58319 29.132480 14.517140 10.63269 14.526337
quickSCnorm(countMat) |> chead() #Single-cell (scran) normalization (return log counts by default)
#>                 a        b        c        d        e
#> TSPAN6   4.390121 3.320394 4.583889 4.520336 3.454398
#> TNMD     3.805213 5.856305 3.583936 4.804103 3.802213
#> DPM1     0.000000 0.000000 0.000000 0.000000 1.995848
#> SCYL3    4.698222 3.457882 2.806395 3.581871 3.454398
#> C1orf112 3.904738 4.952546 3.998950 3.581871 3.994809
RPKM(countMat,gene.length = geneLengthGRCh38) |> chead() #Reads Per Kilobase per Million (RPKM) normalization for full-length transcript, short read sequencing.
#>                  a         b         c         d         e
#> TSPAN6    9.880482  4.448071 11.371905 10.860501  4.929187
#> TNMD     18.090181 79.351552 15.319672 37.544157 18.049702
#> DPM1      0.000000  0.000000  0.000000  0.000000  5.556055
#> SCYL3     8.137438  3.256332  1.954592  3.577827  3.247692
#> C1orf112  5.256509 11.268644  5.636609  4.127063  5.619372
TPMfullLength(countMat,gene.length = geneLengthGRCh38) |> chead() #Transcript per milion (TPM) normalization for full-length transcript, short read sequencing.
#>                  a         b         c         d         e
#> TSPAN6    7.623907  3.391072  8.785572  8.196325  3.864663
#> TNMD     13.958616 60.495183 11.835491 28.334246 14.151627
#> DPM1      0.000000  0.000000  0.000000  0.000000  4.356151
#> SCYL3     6.278952  2.482527  1.510056  2.700155  2.546309
#> C1orf112  4.055990  8.590868  4.354664  3.114658  4.405793
```

## Principal Component Analysis tools

Principal Component Analysis can be performed with `PCA` or approximated
with `fastPCA` with parametrs optimezed by default for RNA-Seq. Result
can be then plotted with `pca2d`.

``` r
data("sampleAnnot") #colData of the count table

PCAres<-PCA(bulkLogCounts)
PCAresQuick<-fastPCA(bulkLogCounts,nPC = 2) # /!\ for fastPCA percentage of explained var is on the total of computed PC (less PC -> increase percentage)
multiplot(
  pca2d(PCAres,returnGraph = TRUE,main = "PCA on bulkLogCounts",colorBy = sampleAnnot$culture_media),
  pca2d(PCAresQuick,returnGraph = TRUE,main = "estimation with fastPCA",colorBy = sampleAnnot$culture_media),
  pca2d(PCAres,plotVars = TRUE,outlierLabel = TRUE,returnGraph = TRUE,main = "Contribution of genes for PC 1 & 2"),
  barplotPercentVar(PCAres,returnGraph = TRUE,nPC = 30)+ggtitle("Scree plot of explained variance per PC"),
  cols = 2
)
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-9-1.png" width="100%" style="display: block; margin: auto;" />

In the gene contribution plot, `pointdensity.nrd` is used for estimating
the 2 density of points. `pcaAddSamples` is another PCA function than
can be used to add new samples to the PCA.

For linking experimental variable to principal components we can perform
a Principal Component analysis of variance.

``` r
library(ggbeeswarm)
resPC_aov<-PCaov(PCAres,colData = sampleAnnot[,c("culture_media","line","passage")])

ggBorderedFactors(
  ggplot(resPC_aov,aes(x=PC,y=sumSqPercent,fill= feature))+
    geom_beeswarm(pch=21,size=4,cex = 3)+
    xlab("Principal component")+ylab("% Sum of Squares")+
    scale_fill_manual(values = oobColors(nlevels(resPC_aov$feature)))+
    theme(
      panel.grid.major.y = element_line(colour = "grey75"),
      panel.grid.minor.y = element_line(colour = "grey75"),
      panel.background = element_rect(fill = NA,colour="black")
    )
)
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

The culture media seems to be linked here with PC 1 & 2. Note in the
code the use of `mostDistantColor` for generating a color scale with
maximum contrast and `ggBorderedFactors` for moving vertical line break
between factors in a *ggplot2* function.

## Characterize features of RNA-Seq

We can first filter the count table by computing an overdispersion score
for each gene with `getMostVariableGenes`. This function must be
executed on a normalized, not log transformed count table. We hence have
to “unlog” it.

``` r
normCount<-2^(bulkLogCounts-1)
dispData<-getMostVariableGenes(normCount,minCount=1)
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

``` r

overdispersedGenes<-rn(dispData)[dispData$residuals>0]
length(overdispersedGenes)
#> [1] 4028
```

We can then compute a PCA on these genes that will result in a reduced
dimension space that will be used for downstream analyses.

``` r
reducedSpace<-t(fastPCA(bulkLogCounts[overdispersedGenes,],nPC = 30))
```

For example, we can perform non linear dimension reductions. Several are
provided with the package: `NMDS`, `UMAP` and `TRIMAP`. The result is
plotted with `proj2d`.

``` r
reducedSpace<-t(fastPCA(bulkLogCounts[overdispersedGenes,],nPC = 30)$x)
chead(reducedSpace)
#>       P3E07LSP32 R3E08MIP01 A3E09MIP20  A3E10MIP25  E3E11LEJ07
#> PC1   0.08257998  -8.905635  18.272929   6.0898056 -89.2020967
#> PC2  42.50375008   2.073339 -16.236969  -8.3071475 -13.3994248
#> PC3   1.07668171 -15.241565 -14.630569 -13.8937491  -0.3845952
#> PC4  -6.53768721  -3.562575   3.176719  -0.4743893  -4.0716714
#> PC5 -15.38921203   1.870331  -3.234557  -1.6699761  -5.8401730
coordNMDS<-NMDS(reducedSpace)
#> initial  value 14.274929 
#> final  value 14.261433 
#> converged
coordUMAP<-UMAP(reducedSpace,n_neighbors = ncol(reducedSpace))
coordTRIMAP<-TRIMAP(reducedSpace)

multiplot(cols = 2,
          pca2d(PCAres,colorBy = sampleAnnot$culture_media,returnGraph = TRUE, main = "PCA",ratioPerPercentVar = TRUE),
          proj2d(coordNMDS,colorBy = sampleAnnot$culture_media,returnGraph = TRUE, main = "NMDS",fixedCoord = TRUE),
          proj2d(coordUMAP,colorBy = sampleAnnot$culture_media,returnGraph = TRUE, main = "UMAP",fixedCoord = TRUE),
          proj2d(coordTRIMAP,colorBy = sampleAnnot$culture_media,returnGraph = TRUE, main = "TRIMAP",fixedCoord = TRUE)
)
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-13-1.png" width="100%" style="display: block; margin: auto;" />

The knn graph used can also be saved when computing the UMAP, and used
for a community detection (Leiden) clustering.

``` r
UMAPwithNN<-UMAP(reducedSpace,n_neighbors = ncol(reducedSpace),ret_nn = TRUE)
sampleAnnot$leidenClusters<-leidenFromUMAP(UMAPwithNN)

proj2d(UMAPwithNN$embedding,colorBy = sampleAnnot$leidenClusters,
       plotFactorsCentroids = TRUE,plotLabelRepel = TRUE,fixedCoord = TRUE, legendTitle = "Leiden cluster")
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

We now pick the best marker per cluster with `getMarkers`, which use the
`auroc` function for computing marker scores, along with logFC and
p-value.

``` r
markerPerCluster<-getMarkers(bulkLogCounts[overdispersedGenes,],group = sampleAnnot$leidenClusters)
chead(markerPerCluster)
#>            k1.lfc  k1.auroc  k1.score      k1.pval      k1.padj
#> A1BG   -1.0768310 0.1554422 -1.389625 5.855828e-08 2.227316e-07
#> A4GALT  1.1012386 0.9578231  1.969951 6.869410e-16 9.013024e-15
#> AAK1   -0.1796876 0.4027211 -0.233799 1.857275e-01 2.166552e-01
#> AARS2   0.8521467 0.8782313  1.428687 8.959301e-10 4.516654e-09
#> AASS   -0.8434000 0.2734694 -0.827172 1.090721e-03 1.883973e-03
```

Now we pick up the best marker per cluster, and plot them. Note the
log10(n+1) scale accessible as a scale for any ggplot with `log10plus1`.
`extractFeatureMarkerData` is used to extract one value per cluster per
gene (by default the marker score).

``` r
markerScorePerCluster<-extractFeatureMarkerData(markerPerCluster)
head(markerScorePerCluster) #matrix of marker scores
#>                k1         k2         k3
#> A1BG   -1.3896249  0.0000000  2.6999342
#> A4GALT  1.9699515 -1.3579995 -0.5425076
#> AAK1   -0.2337990 -0.1989542  0.9002681
#> AARS2   1.4286870 -0.6045051 -1.0212258
#> AASS   -0.8271720  2.3471458 -1.6273016
#> ABCA1   0.7786247 -1.1834598  0.3623288

top20MarkerPerCluster<-apply(markerScorePerCluster, 2, function(x) rownames(markerPerCluster)[whichTop(x, top = 20)] ) |> data.frame() |> as.list()

topGenePerCluster<-sapply(top20MarkerPerCluster,function(x) x[1])

plotExpr(normCount[topGenePerCluster,],group = sampleAnnot$leidenClusters, legendTitle = "Leiden cluster")
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

We then used those marker for plotting an heatmap with the dedicated
function.

``` r
htMarker(bulkLogCounts,markerData = markerScorePerCluster,topn = 20 ,group = sampleAnnot$leidenClusters, colData = sampleAnnot[c("culture_media","line")])
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-17-1.png" width="100%" style="display: block; margin: auto;" />
Note that we can obtain the same result with more struggles using the
more general function `heatmap.DM`:

``` r
#convert the list to a factor
top20MarkerFactor<-VectorListToFactor(top20MarkerPerCluster)

#We subsampled the matrix so each cluster give the same number of sample in the heatmap
markEXpressionSubSampled<-subSampleColumnPerGroup(bulkLogCounts[names(top20MarkerFactor),],groupVector = sampleAnnot$leidenClusters)

heatmap.DM(markEXpressionSubSampled, colData = sampleAnnot[cn(markEXpressionSubSampled),c("culture_media","line")],
           column_split=sampleAnnot[cn(markEXpressionSubSampled),"leidenClusters"],row_split=top20MarkerFactor,
           cluster_row_slices = FALSE, cluster_column_slices =FALSE)
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-18-1.png" width="100%" style="display: block; margin: auto;" />

``` r
#> 
```

`autoGparFontSizeMatrix` is used by `heatmap.DM` to adjust row and
column size to the number of samples / features. `heatmap.DM` is a
wrapper for `ComplexHeatmap`. It uses `genColorsForAnnots` and
`genTopAnnot` are for generating top annotations. The proper heatmap
color scale is generated using `computeColorScaleFun` which is able to
map colors to percentiles.

Example:

``` r
values=sort(rnorm(100))
computeColorScaleFun(colors = c("blue","white","red"),values = values,returnColorFun = FALSE,useProb = TRUE, probs = c(.25,.5,.75)) |> 
  plotPalette()
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

`computeColorScaleFun` can also return a ggplot scale if
`returnGGscale`.

*oob* features other functions for colors: - `oobColors`: Generate
pretty and distinct values for a qualitative scale.

- `mostDistantColor`: Most distant theoretical color palette for a
  qualitative scale (wrapper for *qualpal*).

- `ggplotColours`: generate the default colors of *ggplot2* for a
  qualitative scale.

- `convertColorAdd2Sub`: Convert color from additive to subtracting
  mixing.

Finally, `customUpsetPlot` is used to compare intersection between gene
set with an additional value of enrichment. Here we compare the marker
of unsupervised cluster and differential expression results between
naive and primed pluripotency.

``` r
sets = list(k1 = rn(markerPerCluster)[markerPerCluster$k1.padj < 0.05],
            k2 = rn(markerPerCluster)[markerPerCluster$k2.padj < 0.05],
            k3= rn(markerPerCluster)[markerPerCluster$k3.padj < 0.05],
            primeVSnaiveUP=rn(DEgenesPrime_Naive)[DEgenesPrime_Naive$isDE=="UPREG"],
            primeVSnaiveDOWN=rn(DEgenesPrime_Naive)[DEgenesPrime_Naive$isDE=="DOWNREG"])

richUpset(sets, universe = rn(markerPerCluster))
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-20-1.png" width="100%" style="display: block; margin: auto;" />

`volcanoPlot.DESeq2` is designed to show the result of differential
expression analysis in the most meaningful way.

``` r
data("DEgenesPrime_Naive")
volcanoPlot.DESeq2(DEgenesPrime_Naive,formula = "~culture_media+Run",condColumn = "culture_media",downLevel = "KSR+FGF2",upLevel = "T2iLGO")
```

<img src="D:/PostdocUnsync/myLibraries/oob/ReadMe_files/figure-gfm/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

## Other utilities functions:

- `corGeneToOthers` : Compute the correlation of a gene to all other
  genes in a matrix.

- `best.cutree`: Find the best cutree for a hierarchical clustering.

- `hierarchicalClustering`: Hierarchical clustering from a matrix.

- `supprNAnames`: Delete row/column in a matrix/df with NA row/colnames

- `takefirst`: Similar to unique but conserve vector names or return
  index where you can find each first value of multiple element.

- `formatAnnotFromMeta`: Format a data frame following metadata. Also
  add an attribute “colorScale” for mapping colors to the features of
  the data frame.

- `formatNumber2Character`: Convert numeric to string, add 0 to the
  number to respect lexicographical order.

- `strsplitNth`: Character split with chosen returned element, return a
  vector.

- `make.unique2`: Similar to make.unique, but also add a sequence member
  for the first encountered duplicated element.

- `matrixCoord1D_2D`: Return the row and column index (2D coordinate)
  from a 1D coordinate in a matrix.

- `matrixFromDimnames`: Return the row and column index (2D coordinate)
  from a 1D coordinate in a matrix.

## Session info

``` r
sessionInfo()
#> R version 4.3.2 (2023-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 11 x64 (build 22631)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=French_France.utf8  LC_CTYPE=French_France.utf8    LC_MONETARY=French_France.utf8 LC_NUMERIC=C                   LC_TIME=French_France.utf8    
#> 
#> time zone: Europe/Helsinki
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] rmarkdown_2.25   ggbeeswarm_0.7.2 MASS_7.3-60.0.1  oob_0.9.0        rlang_1.1.3      ggplot2_3.4.4   
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.3.2               later_1.3.2                 batchelor_1.18.1            bitops_1.0-7                filelock_1.0.3              tibble_3.2.1                basilisk.utils_1.14.1       lifecycle_1.0.4            
#>   [9] edgeR_4.0.9                 doParallel_1.0.17           rprojroot_2.0.4             processx_3.8.3              lattice_0.22-5              qualpalr_0.4.4              SnowballC_0.7.1             magrittr_2.0.3             
#>  [17] sass_0.4.8                  limma_3.58.1                jquerylib_0.1.4             yaml_2.3.8                  remotes_2.4.2.1             metapod_1.10.1              httpuv_1.6.13               sessioninfo_1.2.2          
#>  [25] pkgbuild_1.4.3              reticulate_1.34.0           RColorBrewer_1.1-3          ResidualMatrix_1.12.0       abind_1.4-5                 pkgload_1.3.4               zlibbioc_1.48.0             GenomicRanges_1.54.1       
#>  [33] purrr_1.0.2                 BiocGenerics_0.48.1         RCurl_1.98-1.14             rgl_1.2.8                   circlize_0.4.15             GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2           
#>  [41] ggrepel_0.9.5               irlba_2.3.5.1               dqrng_0.3.2                 DelayedMatrixStats_1.24.0   codetools_0.2-19            DelayedArray_0.28.0         scuttle_1.12.0              tidyselect_1.2.0           
#>  [49] shape_1.4.6                 farver_2.1.1                randtoolbox_2.0.4           ScaledMatrix_1.10.0         viridis_0.6.4               matrixStats_1.2.0           stats4_4.3.2                base64enc_0.1-3            
#>  [57] jsonlite_1.8.8              GetoptLong_1.0.5            BiocNeighbors_1.20.2        ellipsis_0.3.2              survival_3.5-7              scater_1.30.1               iterators_1.0.14            foreach_1.5.2              
#>  [65] dbscan_1.1-12               tools_4.3.2                 snow_0.4-4                  Rcpp_1.0.12                 glue_1.7.0                  gridExtra_2.3               SparseArray_1.2.3           xfun_0.41                  
#>  [73] MatrixGenerics_1.14.0       usethis_2.2.2               GenomeInfoDb_1.38.5         dplyr_1.1.4                 withr_3.0.0                 fastmap_1.1.1               basilisk_1.14.2             bluster_1.12.0             
#>  [81] fansi_1.0.6                 callr_3.7.3                 digest_0.6.34               rsvd_1.0.5                  R6_2.5.1                    mime_0.12                   colorspace_2.1-0            Cairo_1.6-2                
#>  [89] scattermore_1.2             utf8_1.2.4                  generics_0.1.3              data.table_1.14.10          FNN_1.1.4                   htmlwidgets_1.6.4           S4Arrays_1.2.0              rngWELL_0.10-9             
#>  [97] uwot_0.1.16                 pkgconfig_2.0.3             gtable_0.3.4                ComplexHeatmap_2.18.0       SingleCellExperiment_1.24.0 XVector_0.42.0              htmltools_0.5.7             profvis_0.3.8              
#> [105] clue_0.3-65                 scales_1.3.0                Biobase_2.62.0              png_0.1-8                   scran_1.30.1                knitr_1.45                  rstudioapi_0.15.0           reshape2_1.4.4             
#> [113] rjson_0.2.21                cachem_1.0.8                GlobalOptions_0.1.2         stringr_1.5.1               parallel_4.3.2              miniUI_0.1.1.1              vipor_0.4.7                 desc_1.4.3                 
#> [121] pillar_1.9.0                grid_4.3.2                  vctrs_0.6.5                 urlchecker_1.0.1            lsa_0.73.3                  promises_1.2.1              BiocSingular_1.18.0         beachmat_2.18.0            
#> [129] xtable_1.8-4                cluster_2.1.6               beeswarm_0.4.0              evaluate_0.23               magick_2.8.2                cli_3.6.2                   locfit_1.5-9.8              compiler_4.3.2             
#> [137] crayon_1.5.2                labeling_0.4.3              RcppArmadillo_0.12.6.6.1    ps_1.7.6                    plyr_1.8.9                  fs_1.6.3                    stringi_1.8.3               viridisLite_0.4.2          
#> [145] BiocParallel_1.36.0         assertthat_0.2.1            munsell_0.5.0               devtools_2.4.5              glmnet_4.1-8                Matrix_1.6-5                dir.expiry_1.10.0           sparseMatrixStats_1.14.0   
#> [153] statmod_1.5.0               shiny_1.8.0                 SummarizedExperiment_1.32.0 highr_0.10                  igraph_1.6.0                memoise_2.0.1               bslib_0.6.1
```

