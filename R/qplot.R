#' GGplot2 quick plot
#'
#' @param x,y,... Aesthetics passed into each layer
#' @param data Data frame to use (optional).  If not specified, will create
#'   one, extracting vectors from the current environment.
#' @param facets faceting formula to use. Picks [facet_wrap()] or
#'   [facet_grid()] depending on whether the formula is one-
#'   or two-sided
#' @param margins See `facet_grid()`: display marginal facets?
#' @param geom Character vector specifying geom(s) to draw. Defaults to
#'  "point" if x and y are specified, and "histogram" if only x is specified.
#' @param xlim,ylim X and y axis limits
#' @param log Which variables to log transform ("x", "y", or "xy")
#' @param main,xlab,ylab Character vector (or expression) giving plot title,
#'   x axis label, and y axis label respectively.
#' @param asp The y/x aspect ratio
#' @export
#' @details This function is a copy of ggplot2::qplot, as it is deprecated since ggplot2 >=3.4 but is still useful.
#' @examples
#' # Use data from data.frame
#' oobqplot(mpg, wt, data = mtcars)
#' oobqplot(mpg, wt, data = mtcars, colour = cyl)
#' oobqplot(mpg, wt, data = mtcars, size = cyl)
#' oobqplot(mpg, wt, data = mtcars, facets = vs ~ am)

oobqplot<-function (x, y, ..., data, facets = NULL, margins = FALSE, geom = "auto",
										xlim = c(NA, NA), ylim = c(NA, NA), log = "", main = NULL,
										xlab = NULL, ylab = NULL, asp = NA)
{
	caller_env <- parent.frame()
	check_character(geom)
	exprs <- enquos(x = x, y = y, ...)
	is_missing <- vapply(exprs, quo_is_missing, logical(1))
	is_constant <- (!names(exprs) %in% ggplot_global$all_aesthetics) |
		vapply(exprs, quo_is_call, logical(1), name = "I")

	mapping <- new_aes(exprs[!is_missing & !is_constant], env = parent.frame())
	consts <- exprs[is_constant]
	aes_names <- names(mapping)
	mapping <- rename_aes(mapping)
	if (is.null(xlab)) {
		if (quo_is_missing(exprs$x)) {
			xlab <- ""
		}
		else {
			xlab <- as_label(exprs$x)
		}
	}
	if (is.null(ylab)) {
		if (quo_is_missing(exprs$y)) {
			ylab <- ""
		}
		else {
			ylab <- as_label(exprs$y)
		}
	}
	if (missing(data)) {
		data <- data.frame()
		facetvars <- all.vars(facets)
		facetvars <- facetvars[facetvars != "."]
		names(facetvars) <- facetvars
		facetsdf <- as.data.frame(mget(facetvars, envir = caller_env))
		if (nrow(facetsdf))
			data <- facetsdf
	}
	if ("auto" %in% geom) {
		if ("sample" %in% aes_names) {
			geom[geom == "auto"] <- "qq"
		}
		else if (missing(y)) {
			x <- eval_tidy(mapping$x, data, caller_env)
			if (is.discrete(x)) {
				geom[geom == "auto"] <- "bar"
			}
			else {
				geom[geom == "auto"] <- "histogram"
			}
			if (is.null(ylab))
				ylab <- "count"
		}
		else {
			if (missing(x)) {
				mapping$x <- quo(seq_along(!!mapping$y))
			}
			geom[geom == "auto"] <- "point"
		}
	}
	p <- ggplot(data, mapping, environment = caller_env)
	if (is.null(facets)) {
		p <- p + facet_null()
	}
	else if (is.formula(facets) && length(facets) == 2) {
		p <- p + facet_wrap(facets)
	}
	else {
		p <- p + facet_grid(rows = deparse(facets), margins = margins)
	}
	if (!is.null(main))
		p <- p + ggtitle(main)
	for (g in geom) {
		params <- lapply(consts, eval_tidy)
		p <- p + do.call(paste0("geom_", g), params)
	}
	logv <- function(var) var %in% strsplit(log, "")[[1]]
	if (logv("x"))
		p <- p + scale_x_log10()
	if (logv("y"))
		p <- p + scale_y_log10()
	if (!is.na(asp))
		p <- p + theme(aspect.ratio = asp)
	if (!missing(xlab))
		p <- p + xlab(xlab)
	if (!missing(ylab))
		p <- p + ylab(ylab)
	if (!missing(xlim) && !all(is.na(xlim)))
		p <- p + xlim(xlim)
	if (!missing(ylim) && !all(is.na(ylim)))
		p <- p + ylim(ylim)
	p
}

rename_aes<-function (x) {
	names(x) <- standardise_aes_names(names(x))
	duplicated_names <- names(x)[duplicated(names(x))]
	if (length(duplicated_names) > 0L) {
		cli::cli_warn("Duplicated aesthetics after name standardisation: {.field {unique0(duplicated_names)}}")
	}
	x
}

standardise_aes_names<-function (x) {
	x <- sub("color", "colour", x, fixed = TRUE)
	revalue(x, ggplot_global$base_to_ggplot)
}

revalue<-function (x, replace) {
	if (is.character(x)) {
		replace <- replace[names(replace) %in% x]
		if (length(replace) == 0)
			return(x)
		x[match(names(replace), x)] <- replace
	}
	else if (is.factor(x)) {
		lev <- levels(x)
		replace <- replace[names(replace) %in% lev]
		if (length(replace) == 0)
			return(x)
		lev[match(names(replace), lev)] <- replace
		levels(x) <- lev
	}
	else if (!is.null(x)) {
		stop_input_type(x, "a factor or character vector")
	}
	x
}


new_aes<-function (x, env = globalenv()) {
	check_object(x, is.list, "a {.cls list}")
	x <- lapply(x, new_aesthetic, env = env)
	structure(x, class = "uneval")
}

new_aesthetic<-function (x, env = globalenv()) {
	if (is_quosure(x)) {
		if (!quo_is_symbolic(x)) {
			x <- quo_get_expr(x)
		}
		return(x)
	}
	if (is_symbolic(x)) {
		x <- new_quosure(x, env = env)
		return(x)
	}
}

check_object<-function (x, check_fun, what, ..., allow_null = FALSE, arg = caller_arg(x),
					call = caller_env()){
	if (!missing(x)) {
		if (check_fun(x)) {
			return(invisible(NULL))
		}
		if (allow_null && is_null(x)) {
			return(invisible(NULL))
		}
	}
	stop_input_type(x, as_cli(what), ..., allow_null = allow_null,
									arg = arg, call = call)
}

as_cli<-function (..., env = caller_env())
{
	cli::cli_fmt(cli::cli_text(..., .envir = env))
}

stop_input_type<-function (x, what, ..., allow_na = FALSE, allow_null = FALSE,
					show_value = TRUE, arg = caller_arg(x), call = caller_env()) {
	cli <- env_get_list(nms = c("format_arg", "format_code"),
											last = topenv(), default = function(x) sprintf("`%s`",
																																		 x), inherit = TRUE)
	if (allow_na) {
		what <- c(what, cli$format_code("NA"))
	}
	if (allow_null) {
		what <- c(what, cli$format_code("NULL"))
	}
	if (length(what)) {
		what <- oxford_comma(what)
	}
	message <- sprintf("%s must be %s, not %s.", cli$format_arg(arg),
										 what, obj_type_friendly(x, value = show_value))
	abort(message, ..., call = call, arg = arg)
}

oxford_comma<-function (chr, sep = ", ", final = "or"){
	n <- length(chr)
	if (n < 2) {
		return(chr)
	}
	head <- chr[seq_len(n - 1)]
	last <- chr[n]
	head <- paste(head, collapse = sep)
	if (n > 2) {
		paste0(head, sep, final, " ", last)
	}
	else {
		paste0(head, " ", final, " ", last)
	}
}

check_character<-function (x, ..., allow_null = FALSE, arg = caller_arg(x), call = caller_env()){
	if (!missing(x)) {
		if (is_character(x)) {
			return(invisible(NULL))
		}
		if (allow_null && is_null(x)) {
			return(invisible(NULL))
		}
	}
	stop_input_type(x, "a character vector", ..., allow_na = FALSE,
									allow_null = allow_null, arg = arg, call = call)
}

is.discrete<-function (x){
	is.factor(x) || is.character(x) || is.logical(x)
}

is.formula<-function (x)
	inherits(x, "formula")

obj_type_friendly<-function (x, value = TRUE){
	if (is_missing(x)) {
		return("absent")
	}
	if (is.object(x)) {
		if (inherits(x, "quosure")) {
			type <- "quosure"
		}
		else {
			type <- paste(class(x), collapse = "/")
		}
		return(sprintf("a <%s> object", type))
	}
	if (!is_vector(x)) {
		return(.rlang_as_friendly_type(typeof(x)))
	}
	n_dim <- length(dim(x))
	if (!n_dim) {
		if (!is_list(x) && length(x) == 1) {
			if (is_na(x)) {
				return(switch(typeof(x), logical = "`NA`", integer = "an integer `NA`",
											double = if (is.nan(x)) {
												"`NaN`"
											} else {
												"a numeric `NA`"
											}, complex = "a complex `NA`", character = "a character `NA`",
											.rlang_stop_unexpected_typeof(x)))
			}
			show_infinites <- function(x) {
				if (x > 0) {
					"`Inf`"
				}
				else {
					"`-Inf`"
				}
			}
			str_encode <- function(x, width = 30, ...) {
				if (nchar(x) > width) {
					x <- substr(x, 1, width - 3)
					x <- paste0(x, "...")
				}
				encodeString(x, ...)
			}
			if (value) {
				if (is.numeric(x) && is.infinite(x)) {
					return(show_infinites(x))
				}
				if (is.numeric(x) || is.complex(x)) {
					number <- as.character(round(x, 2))
					what <- if (is.complex(x))
						"the complex number"
					else "the number"
					return(paste(what, number))
				}
				return(switch(typeof(x), logical = if (x) "`TRUE`" else "`FALSE`",
											character = {
												what <- if (nzchar(x)) "the string" else "the empty string"
												paste(what, str_encode(x, quote = "\""))
											}, raw = paste("the raw value", as.character(x)),
											.rlang_stop_unexpected_typeof(x)))
			}
			return(switch(typeof(x), logical = "a logical value",
										integer = "an integer", double = if (is.infinite(x)) show_infinites(x) else "a number",
										complex = "a complex number", character = if (nzchar(x)) "a string" else "\"\"",
										raw = "a raw value", .rlang_stop_unexpected_typeof(x)))
		}
		if (length(x) == 0) {
			return(switch(typeof(x), logical = "an empty logical vector",
										integer = "an empty integer vector", double = "an empty numeric vector",
										complex = "an empty complex vector", character = "an empty character vector",
										raw = "an empty raw vector", list = "an empty list",
										.rlang_stop_unexpected_typeof(x)))
		}
	}
	vec_type_friendly(x)
}

.rlang_as_friendly_type<-function (type) {
	switch(type, list = "a list", `NULL` = "`NULL`", environment = "an environment",
				 externalptr = "a pointer", weakref = "a weak reference",
				 S4 = "an S4 object", name = , symbol = "a symbol", language = "a call",
				 pairlist = "a pairlist node", expression = "an expression vector",
				 char = "an internal string", promise = "an internal promise",
				 ... = "an internal dots object", any = "an internal `any` object",
				 bytecode = "an internal bytecode object", primitive = ,
				 builtin = , special = "a primitive function", closure = "a function",
				 type)
}

.rlang_stop_unexpected_typeof<-function (x, call = caller_env()){
	abort(sprintf("Unexpected type <%s>.", typeof(x)), call = call)
}

vec_type_friendly<-function (x, length = FALSE)
{
	if (!is_vector(x)) {
		abort("`x` must be a vector.")
	}
	type <- typeof(x)
	n_dim <- length(dim(x))
	add_length <- function(type) {
		if (length && !n_dim) {
			paste0(type, sprintf(" of length %s", length(x)))
		}
		else {
			type
		}
	}
	if (type == "list") {
		if (n_dim < 2) {
			return(add_length("a list"))
		}
		else if (is.data.frame(x)) {
			return("a data frame")
		}
		else if (n_dim == 2) {
			return("a list matrix")
		}
		else {
			return("a list array")
		}
	}
	type <- switch(type, logical = "a logical %s", integer = "an integer %s",
								 numeric = , double = "a double %s", complex = "a complex %s",
								 character = "a character %s", raw = "a raw %s", type = paste0("a ",
								 																															type, " %s"))
	if (n_dim < 2) {
		kind <- "vector"
	}
	else if (n_dim == 2) {
		kind <- "matrix"
	}
	else {
		kind <- "array"
	}
	out <- sprintf(type, kind)
	if (n_dim >= 2) {
		out
	}
	else {
		add_length(out)
	}
}

# Environment that holds various global variables and settings for ggplot,
# such as the current theme. It is not exported and should not be directly
# manipulated by other packages.
ggplot_global <- new.env(parent = emptyenv())

# The current theme. Defined here only as placeholder, and defined properly
# in file "theme-current.R". This setup avoids circular dependencies among
# the various source files.
ggplot_global$theme_current <- list()

# Element tree for the theme elements. Defined here only as placeholder, and
# defined properly in file "theme-elements.r".
ggplot_global$element_tree <- list()

# List of all aesthetics known to ggplot
# (In the future, .all_aesthetics should be removed in favor
# of direct assignment to ggplot_global$all_aesthetics, see below.)
.all_aesthetics <- c(
	"adj", "alpha", "angle", "bg", "cex", "col", "color",
	"colour", "fg", "fill", "group", "hjust", "label", "linetype", "lower",
	"lty", "lwd", "max", "middle", "min", "pch", "radius", "sample", "shape",
	"size", "srt", "upper", "vjust", "weight", "width", "x", "xend", "xmax",
	"xmin", "xintercept", "y", "yend", "ymax", "ymin", "yintercept", "z"
)

ggplot_global$all_aesthetics <- .all_aesthetics

# Aesthetic aliases
# (In the future, .base_to_ggplot should be removed in favor
# of direct assignment to ggplot_global$base_to_ggplot, see below.)
.base_to_ggplot <- c(
	"col"   = "colour",
	"color" = "colour",
	"pch"   = "shape",
	"cex"   = "size",
	"lty"   = "linetype",
	"lwd"   = "size",
	"srt"   = "angle",
	"adj"   = "hjust",
	"bg"    = "fill",
	"fg"    = "colour",
	"min"   = "ymin",
	"max"   = "ymax"
)

ggplot_global$base_to_ggplot <- .base_to_ggplot

