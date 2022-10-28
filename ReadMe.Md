# oob: a toolbox for analyzing OMIC data, Out Of Bounds of black boxes.

This package provide function for analyzing OMIC data, with a focus on plots, functional enrichment and bulk/single-cell RNA-Seq.

## Installation

In a R console: `devtools::install("https://github.com/DimitriMeistermann/oob")`

*oob* needs also several python dependencies accessible in the environment used by *reticulate*, the R package for calling python functions. See [https://rstudio.github.io/reticulate/articles/python_packages.html](https://rstudio.github.io/reticulate/articles/python_packages.html) for more information on installation of python package via reticulate.

For a quick installation:
```{r}
reticulate::py_install(c("igraph","numpy","leidenalg","trimap"))
```