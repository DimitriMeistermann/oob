#' Load quickly a text file to a dataframe/matrix
#'
#' @param fileName Path to the file.
#' @param sep Separator ("\\t"=tab-separated values).
#' @param row.names Column that contains the row names,set to NULL for not
#'   attributing row names.
#' @param as.matrix Return a dataframe if FALSE, a matrix if TRUE.
#' @param stringsAsFactors Convert automatically strings to factors if TRUE.
#' @param ... Other parameters passed to data.table::fread.
#'
#' @return A dataframe or a matrix.
#' @export
#' @seealso fastWrite
#' @examples
#' data("bulkLogCounts")
#' tf <- tempfile(fileext = ".tsv")
#' fastWrite(bulkLogCounts,tf)
#' fastRead(tf)
fastRead <-
    function(fileName,
                    sep = '\t',
                    row.names = 1,
                    as.matrix = FALSE,
                    stringsAsFactors = FALSE,
                    ...) {
        dat <-
            as.data.frame(data.table::fread(fileName,
                stringsAsFactors = stringsAsFactors,
                sep = sep,
                ...))
        if (!is.null(row.names)) {
            if (row.names != FALSE) {
                rownames(dat) <- dat[, row.names]
                dat <- dat[, -row.names, drop = FALSE]
            }
        }
        if (as.matrix)
            dat <- as.matrix(dat)
        return(dat)
    }

#' Write quickly a text file from dataframe/matrix
#'
#' @param x matrix/dataframe to save.
#' @param fileName Path to the file to write.
#' @param headRow Name in the header of the rownames column.
#' @param row.names Whether or not to save rownames.
#' @param col.names Whether or not to save colnames.
#' @param dec Character that serve as a decimal separator.
#' @param sep Separator ("\\t"=tab-separated values).
#' @param ... Other parameters passed to data.table::fwrite.
#' @export
#' @seealso fastRead
#' @return NULL, write a file.
#' @examples
#' data("bulkLogCounts")
#' tf <- tempfile(fileext = ".tsv")
#' fastWrite(bulkLogCounts,tf)
fastWrite <-
    function (x,
            fileName = "default.tsv",
            headRow = "Name",
            row.names = TRUE,
            col.names = TRUE,
            dec = ".",
            sep = "\t",
            ...)
    {
        if (is.null(rownames(x)))
            row.names <- FALSE
        if (is.null(colnames(x)))
            col.names <- FALSE
        if (row.names) {
            x <- cbind(rownames(x), x)
            if (col.names)
                colnames(x)[1] <- headRow
        }
        data.table::fwrite(
            x = data.frame(x),
            file = fileName,
            sep = sep,
            row.names = FALSE,
            col.names = col.names,
            quote = FALSE,
            dec = dec,
            ...
        )
    }

#' Write a list in a text file.
#'
#' @description Write a list in a text file with a specific separator. If
#' list.names and vector.names are TRUE each element of the list is saved on a
#' cycle of 3 rows. First row = element name, second = vector names, third =
#' vector.
#'
#' @param list The list object to save.
#' @param filename Path to the file to write.
#' @param sep Separator ("\\t"=tab-separated values).
#' @param list.names Save list names in the row before the values.
#' @param vector.names Save value name vector in the row before the values.
#' @return NULL, write a file.
#' @export
#' @seealso read.vectorList
#' @examples
#' tf <- tempfile(fileext = ".tsv")
#' write.vectorList(list(1:10,letters[1:10]),tf)
write.vectorList <-
    function(list,
            filename,
            sep = "\t",
            list.names = TRUE,
            vector.names = FALSE) {
        if ((!is.list(list)))
            stop("list must be a list")
        sink(filename)
        for (i in seq_along(list)) {
            if (!(is.null(names(list)) |
                (!list.names)))
                cat(names(list[i]), "\n", sep = "")
            element <- list[[i]]
            if (!(is.vector(element) |
                    is.factor(element)))
                stop("each element of the list should be a vector")
            if (!(is.null(names(element)) |
                    (!vector.names)))
                cat(paste0(names(element), collapse = sep), "\n", sep = "")
            cat(paste0(as.character(element), collapse = sep), "\n", sep = "")
        }
        sink()
    }

#' Load a list from a text file.
#'
#' @description Currently works for file saved with parameters list.names=TRUE
#' and vector.names=FALSE from write.vectorList.
#'
#' @param fileName Path to the file to load
#' @param sep Separator ("\\t"=tab-separated values).
#'
#' @return A list.
#' @export
#' @seealso write.vectorList
#' @examples
#' tf <- tempfile(fileext = ".tsv")
#' write.vectorList(list(1:10,letters[1:10]),tf)
#' read.vectorList(tf)
read.vectorList <- function(fileName, sep = "\t") {
    con <- file(fileName)
    txt <- readLines(con)
    elNames <- str_remove_all(txt[seq(1, length(txt), 2)], sep)
    res <- strsplit(txt[seq(2, length(txt), 2)], sep)
    names(res) <- elNames
    close(con)
    res
}
