#' Compute p-value and other metric for enrichment of set intersection.
#'
#' @param intersectionSize Numeric. Number of common element between sets.
#' @param setSizes Vector of numeric. Number of element in each set.
#' @param universeSize Total number of element in the universe.
#'
#' @details If two sets, perform an hypergeometric test, if more the p-value is
#' computed from binomial distribution.
#'
#' @return A list of values:
#' - observed: observed number of elements in intersection
#'  same value as `intersectionSize`)
#' - expected: expected number of elements in intersection.
#' - OEdeviation: Effect size of the association
#'  `(observed - expected) / sqrt(universeSize)`
#' - pval: pvalue (upper tail) where NULL hypotheses is that observed
#'  intersection is from the distribution of intersection size
#'  with p = nExpected / nUniverse.
#'
#' @export
#'
#' @examples
#' setSizes<-c(150,200,250)
#' universeSize<-10000
#' intersectionSize<-5
#' enrichSetIntersection(intersectionSize, setSizes, universeSize)
enrichSetIntersection <-
    function(intersectionSize, setSizes, universeSize) {
        nSet <- length(setSizes)
        expected <- prod(setSizes) / (universeSize ^ (nSet - 1))
        if (nSet == 2) {
            pval <-
                phyper(
                    q = intersectionSize - 0.5,
                    m = setSizes[1],
                    n = universeSize - setSizes[1],
                    k = setSizes[2],
                    lower.tail = FALSE
                )
        } else if (nSet > 2) {
            pval <-
                pbinom(intersectionSize - 1,
                    universeSize,
                    expected / universeSize,
                    lower.tail = FALSE)
        } else{
            stop("Number of set should be superior to 1")
        }
        OEdeviation <- (intersectionSize - expected) / sqrt(universeSize)
        return(
            list(
                "observed" = intersectionSize,
                "expected" = expected,
                "OEdeviation" = OEdeviation,
                "pval" = pval
            )
        )
    }


#' Compute enrichment statistics for set intersection.
#'
#' @param featurePerGroupList A list of sets (vector of character)
#' @param universe NULL or vector of Character. The entire list of features
#'   (Universal set)
#'
#' @return A list that contains:
#'
#' 1. A data frame (*enrichStats*) containing the enrichment statistics for each
#' intersection of sets. Row names are corresponding to the intersection code
#' (e.g. 101 for the intersection of the first and third set). The columns are:
#'    - `observed`: observed number of elements in intersection
#'    - `expected`: expected number of elements in intersection
#'    - `OEdeviation`: Effect size of the association
#'      `(observed - expected) / sqrt(universeSize)`
#'    - `pval`: pvalue (upper tail) where NULL hypotheses is that observed
#'   intersection is from the distribution of intersection size
#'   with *p = nExpected / nUniverse*.
#'    - `padj`: p-value adjusted for multiple testing
#'    - `log10padj`: -log10(padj), with a floor at 384
#'    - `combDegree`: Number of sets in the intersection
#'    - `combsize`: Number of elements in the intersection
#'
#' 2. A combination matrix (*upsetMatrix*), can be used as input to
#' [ComplexHeatmap::UpSet]
#'
#' 3. *universeSize*: Total number of element in the universe.
#'
#' 4. *setsize*: Number of elements in each set.
#'
#' @export
#'
#' @examples
#' lt <- list(set1 = sample(letters, 5),
#'                     set2 = sample(letters, 10),
#'                     set3 = sample(letters, 15)
#'                     )
#' lt$set4 <- unique(c(lt$set1,lt$set2))
#'
#' richUpsetStats(lt, universe = letters)
#'
#' @seealso [richUpset()]
richUpsetStats<-function(featurePerGroupList, universe=NULL){
    if (is.null(universe))
        universe <- unique(unlist(featurePerGroupList))
    isInGroupMatrix <-
        list_to_matrix(featurePerGroupList, universal_set = universe)
    upsetMatrix <-
        make_comb_mat(isInGroupMatrix, mode = "intersect")
    upsetMatrix <-
        upsetMatrix[comb_degree(upsetMatrix) > 1]
    # retain only intersections of sets

    combsize <- comb_size(upsetMatrix)
    setsize <- set_size(upsetMatrix)

    #Are the intersections sets (or Venn diagram region) enriched or not ?
    regionEnrich <-
        lapply(comb_name(upsetMatrix), function(region) {
            colOfcomp <- which(strsplit(region, split = "")[[1]] == "1")
            enrichSetIntersection(combsize[region], setsize[colOfcomp],
                length(universe))
        })
    regionEnrichRes <-
        data.frame(row.names = comb_name(upsetMatrix))
    for (el in names(regionEnrich[[1]]))
        regionEnrichRes[[el]] <- vapply(regionEnrich, function(x)
            x[[el]], numeric(1))
    rm(regionEnrich)
    regionEnrichRes$padj <-
        p.adjust(regionEnrichRes$pval, method = "BH")
    regionEnrichRes$log10padj <- -log10(regionEnrichRes$padj)
    regionEnrichRes$log10padj[regionEnrichRes$log10padj == Inf] <-
        384

    regionEnrichRes$combDegree <- comb_degree(upsetMatrix)
    regionEnrichRes$combsize <- combsize
    return(list(
        "enrichStats" = regionEnrichRes,
        "upsetMatrix" = upsetMatrix,
        "universeSize" = length(universe),
        "setsize" = setsize
    ))
}


#' Upset plot with additional enrichment values.
#'
#' @inheritParams richUpsetStats
#' @param pvalThreshold Numeric. The threshold for the pvalue, represented by a
#'   horizontal red bar. Default is 0.01.
#'
#' @description In addition to Upset plot, this method computes and represents
#' additional values useful for understanding the relationship between sets. The
#' main ones are a p-value for each overlap, and the effect size of the
#' association as the OEdev: `(observed - expected) / sqrt(universeSize)`
#'
#'
#' @return Plot in the current graphical device. In the pval bar graph, the red
#' line indicates an adjusted pval of 0.01 (-log10 = 2). U indicates the
#' universe size (total number of elements).
#' @export
#'
#' @examples
#' lt <- list(set1 = sample(letters, 5),
#'                     set2 = sample(letters, 10),
#'                     set3 = sample(letters, 15)
#'                     )
#' lt$set4 <- unique(c(lt$set1,lt$set2))
#'
#' richUpset(lt, universe = letters)
#'
#' @seealso [richUpsetStats()]
richUpset <-
    function(featurePerGroupList,
            universe = NULL,
            pvalThreshold = 0.01) {

        richUpsetStatsRes <- richUpsetStats(featurePerGroupList, universe)
        regionEnrichRes<-richUpsetStatsRes$enrichStats
        universeSize <- richUpsetStatsRes$universeSize
        setsize <- richUpsetStatsRes$setsize
        combsize <- regionEnrichRes$combsize
        upsetMatrix <- richUpsetStatsRes$upsetMatrix

        combDegree <- paste0(regionEnrichRes$combDegree , "\u00b0")

        OEdevMatrix <- regionEnrichRes[, "OEdeviation", drop = FALSE] |> t()
        rownames(OEdevMatrix) <- "Observed /\nExpected\ndeviation"

        enrich_ha <- HeatmapAnnotation(
            "pval" = anno_barplot(
                regionEnrichRes$log10padj,
                gp = gpar(fill = "black"),
                height = unit(3, "cm"),
                axis_param = list(side = "left"),
                ylim = c(0, max(
                    max(regionEnrichRes$log10padj) * 1.1, 2
                ))
            ),
            annotation_name_side = "left",
            annotation_name_rot = 0,
            annotation_name_gp = gpar(fontface = "bold"),
            annotation_label = c("-log10\nadj.\np-val")
        )
        intersect_ha <- HeatmapAnnotation(
            "intersection_size" = anno_barplot(
                combsize,
                gp = gpar(fill = "black"),
                height = unit(3, "cm"),
                axis_param = list(side = "left"),
                ylim = c(0, max(combsize) * 1.1)
            ),
            annotation_name_side = "left",
            annotation_name_rot = 0,
            annotation_name_gp = gpar(fontface = "bold"),
            annotation_label = "Observed \nintersection\nsize\n(expected)"
        )
        set_size_ha <- rowAnnotation(
            "set_size" = anno_barplot(
                setsize,
                gp = gpar(fill = "black"),
                width = unit(2, "cm"),
                ylim = c(0, max(setsize) * 1.3)
            ),
            annotation_name_side = "bottom",
            annotation_name_rot = 0,
            annotation_name_gp = gpar(fontface = "bold"),
            annotation_label = paste0("Set\nsize\n(U=", universeSize, ")")
        )

        ht <- draw(
            UpSet(
                upsetMatrix,
                top_annotation = intersect_ha,
                bottom_annotation = enrich_ha,
                right_annotation = set_size_ha,
                border = TRUE,
                column_split = combDegree,
                row_names_gp = gpar(fontsize = min(1 / max(
                    nchar(rownames(upsetMatrix))
                ) * 260, 20))#automatic fontsize to avoid out of bound text
            )
            %v% Heatmap(
                OEdevMatrix,
                show_column_names = FALSE,
                cell_fun = function(j, i, x, y, w, h, col) {
                    if (colSums(col2rgb(col)) < 382.5)
                        col <- "white"
                    else
                        col <- "black"
                    grid.text(
                        as.character(
                            signif(
                                regionEnrichRes[j, "OEdeviation"], 2
                            )
                        ), x, y,
                        gp = gpar(col = col, fontface = 2)
                    )
                },
                show_heatmap_legend = FALSE,
                col = colorRamp2(
                    colors=c("darkblue", "white", "red2"),
                    breaks = c(-1,0,1)
                ),
                rect_gp = gpar(col = "black"),
                row_names_side = "left",
                row_names_gp = gpar(fontface = 2)
            )
        )

        #Offset to counterbalance column split space
        colPerSplit <- vapply(column_order(ht), length, 1L)
        offsetPerSplit <- seq(0, length(colPerSplit) - 1)
        offsets <-
            unlist(lapply(seq_along(colPerSplit), function(i)
                rep(offsetPerSplit[i], colPerSplit[i])), use.names = FALSE)

        rowOrder <- rev(row_order(ht)[[1]])
        columnOrder <- unlist(column_order(ht))
        xCoordinatesColText <-
            unit(seq_along(columnOrder), "native") + unit(offsets, "mm")

        decorate_annotation("intersection_size", {
            grid.text(
                paste0(
                    combsize[columnOrder],
                    " (",
                    round(regionEnrichRes$expected[columnOrder], 1),
                    ")"
                ),
                x = xCoordinatesColText,
                y = unit(combsize[columnOrder], "native") + unit(6, "pt"),
                default.units = "native",
                just = "center",
                gp = gpar(fontsize = 8)
            )
        })
        decorate_annotation("pval", {
            grid.text(
                round(regionEnrichRes$log10padj[columnOrder], 3),
                x = xCoordinatesColText,
                y = unit(regionEnrichRes$log10padj[columnOrder], "native") +
                    unit(6, "pt"),
                default.units = "native",
                just = "center",
                gp = gpar(fontsize = 8)
            )
        })


        bestPval <- min(regionEnrichRes$padj)

        decorate_annotation("pval",
                            addHbarUpset(
                                -log10(pvalThreshold),
                                offsetPerSplit,
                                colPerSplit,
                                gp = gpar(alpha = .5, col = "red")
                            ))

        decorate_annotation("set_size", {
            grid.text(
                round(setsize[rowOrder], 2),
                y = seq_along(setsize),
                x = unit(setsize[rowOrder], "native") + unit(7, "pt"),
                default.units = "native",
                just = "center",
                gp = gpar(fontsize = 10),
                rot = -90
            )
        })
    }

addHbarUpset <- function(y, offsetPerSplit, colPerSplit, gp = NULL) {
    return({
        y <- unit(y, "native")

        x1 <- cumsum(colPerSplit) + 0.5
        x0 <- c(0.5, x1[seq_len(length(x1) - 1)])

        x0 <- unit(x0, "native") + unit(offsetPerSplit, "mm")
        x1 <- unit(x1, "native") + unit(offsetPerSplit, "mm")

        grid.segments(
            x0 = x0,
            x1 = x1,
            y0 = y,
            y1 = y,
            gp = gp
        )
    })
}
