
rsilhouette <- function(clustering, n = 1000, ...) {
    vapply(seq_len(n),
        function(x, ...) {
            rand_clust <- sample(clustering, length(clustering))
            rsil <- cluster::silhouette(as.integer(rand_clust), ...)
            return(summary(rsil)$avg.width)
        }, 0, ...)
}

condcomp_cluster <-
    function(cond, n = 1000, rand.sil = NULL, return.full = FALSE, ...) {
    if (length(levels(cond)) > 2) {
        stop("Too many conditions. Maximum is 2.")
    }
    true_sil_full <- cluster::silhouette(as.integer(cond), ...)
    if (length(true_sil_full) > 1) {
        true_sil <- summary(true_sil_full)$avg.width
        if (is.null(rand.sil)) {
            rand.sil <- rsilhouette(cond, n=n, ...)
        }
        pval <- sum(rand.sil > true_sil) / length(rand.sil)
        all.sils <- c(rand.sil, true_sil)
        iqr <- "Same"
        if (length(all.sils) %in%
            which(outliers::scores(all.sils, type="iqr", lim=TRUE))) {
            iqr <- "Diff"
        }
        zscore <-
            (true_sil - mean(rand.sil)) / stats::sd(rand.sil)
    } else {
        true_sil <- NA
        zscore <- NA
        pval <- NA
        iqr <- NA
        rand.sil <- NA
    }
    res <- list()
    for (c in levels(cond)) {
        res[paste0(c, "_perc")] <- sum(cond == c) / length(cond)
    }
    for (c in levels(cond)) {
        res[paste0(c, "_ratio")] <- sum(cond == c) / sum(cond != c)
    }
    res$true_sil <- true_sil
    res$zscore <- zscore
    res$pval <- pval
    res$iqr <- iqr
    res$rand.sil <- rand.sil
    if (return.full) {
        res$true_sil_full <- true_sil_full
    }
    return(res)
}

#' Comparison of data conditions in a clustering
#'
#' Performs a condition comparison on a given \code{clustering}. The comparison
#' is performed on each cluster separately between each condcition
#' (\code{cond}). Several statistics are used and, when analysed in conjunction,
#' they might give some insight regarding the heterogeneity of some of the
#' clusters.
#'
#' For a given cluster, several metrics are computed, see the 'Return' section
#' for details about each metric. Some metrics make use of Random Silhouettes,
#' which is defined as follows: given a labeled data set, assign a random label
#' (from the set of labels) to each data point without changing the original
#' ratio of groups. Then compute the \code{\link{silhouette}} index for this
#' data considering these randomly assigned labels, the average silhouette width
#' is the Random Silhouette for the data (with randomly assigned labels). Being
#' an stochastic process, the Monte Carlo approach is applied which gives a
#' vector of several Random Silhouettes.
#'
#' @param clustering A clustering of the data.
#' @param cond A factor indicating the condition which each data point is
#'     subject to.
#' @param dmatrix A distance matrix describing the data to be analysed.
#' @param n The number of random silhouettes to be performed. Keep in mind that
#'     the computation of several random silhouettes is the bottleneck of this
#'     process.
#' @param remove.na Logical. Remove lines with NA (i.e. clusters which the
#'     silhouette could not be computed).
#' @return A data frame with various statistics regarding data heterogeneity
#'     inside each cluster.
#'
#' Each row of the data frame contains several metrics regarding the conditions
#' found in an specific cluster. For now only a maximum of two conditions are
#' supported. These metrics are described below:
#'
#' \describe{
#' \item{x_perc}{Numeric. The percentage of data points belonging to condition
#' 'x'.}
#' \item{x_ratio}{Numeric. The ratio of data points belonging to condition 'x'.
#' For example, considering another condition 'y', the 'x_ratio' would be
#' computed as \eqn{x_perc / y_perc}.}
#' \item{true_sil}{Numeric. True silhouette. The silhouette for the data in
#' this cluster considering the conditions, as defined by the parameter
#' \code{cond}, as groups.}
#' \item{zscore}{Numeric. The Z-score computed based on the
#' \code{\link[cluster]{silhouette}}. See the 'Details' section.}
#' \item{pval}{Numeric. The p-value for 'true_sil'. Computed from the number of
#' Random Silhouettes (see 'Details') that are greater than the 'true_sil' for
#' this cluster.}
#' \item{iqr}{Factor. Interquartile Range (\link[outliers:scores]{IQR}) based
#' outlier detection. Considering the vector including the random silhouettes
#' (see 'Details') and the 'true_sil', the method checks whether 'true_sil' is
#' an outlier in said vector. This will be set to 'Diff' in case 'true_sil' is
#' an outlier or 'Same' otherwise.}
#' }
#' @examples
#' clustering <- iris$Species
#' dmatrix <- as.matrix(dist(iris[-length(iris)]))
#' # Suppose the conditions are 'young' and 'old' fish
#' cond <- sample(c("young", "old"), nrow(iris), replace=TRUE)
#' comp <- condcomp(clustering, cond, dmatrix=dmatrix, n=10)
#' @export
condcomp <-
    function(
        clustering,
        cond,
        dmatrix,
        n = 1000,
        remove.na = TRUE) {
    if (length(levels(cond)) > 2) {
        stop("Too many conditions. Maximum is 2.")
    }
    if (n <= 0) {
        stop("'n' must be positive.")
    }
    if (!is.matrix(dmatrix)) {
        dmatrix <- as.matrix(dmatrix)
    }
    if (!is.factor(cond)) {
        cond <- as.factor(cond)
    }
    if (!is.factor(clustering)) {
        clustering <- as.factor(clustering)
    }
    res <- data.frame(row.names=levels(clustering))
    res[paste0(levels(cond), "_perc")] <-
        lapply(seq_along(length(levels(cond))), function(x) numeric())
    res[paste0(levels(cond), "_ratio")] <-
        lapply(seq_along(length(levels(cond))), function(x) numeric())
    res["true_sil"] <- numeric()
    res["zscore"] <- numeric()
    res["pval"] <- numeric()
    res["iqr"] <- factor(levels=c("Diff", "Same"))
    for (c in levels(clustering)) {
        message("\nProcessing cluster ", c)
        cells_use <- clustering == c
        cond.sub <- cond[cells_use]
        dmtx <- dmatrix[cells_use, cells_use]
        clust.comp <-
        condcomp_cluster(cond.sub, n=n, dmatrix=dmtx)
        res[c, ] <- clust.comp[-length(clust.comp)]
    }
    if (remove.na) {
        res <- subset(res, !is.na(res$true_sil))
    }
    return(res)
}

#' Plots the data frame of comparison of data conditions
#'
#' This function takes the output from \code{\link{condcomp}} and plots some of
#' its attributes in a scatter plot.
#'
#' The first condition ratio that appears in the data frame will be plotted in
#' the y-axis (-log10 scale), whereas the Z-score will be plotted along the
#' x-axis. Each group will be colored by their respective IQR value as shown
#' in the legend.
#'
#' @param ccomp A data frame output from \code{\link{condcomp}}.
#' @param col Color parameter to be used. The default is to color according to
#'     the IQR column of \code{ccomp}.
#' @param main Character vector (or expression) giving plot title.
#' @param legend.title Character vector giving the legend title.
#' @return A ggplot2 object.
#' @examples
#' clustering <- iris$Species
#' dmatrix <- as.matrix(dist(iris[-length(iris)]))
#' # Suppose the conditions are 'young' and 'old' fish
#' cond <- sample(c("young", "old"), nrow(iris), replace=TRUE)
#' comp <- condcomp(clustering, cond, dmatrix=dmatrix, n=10)
#' condcompPlot(comp)
#' @export
condcompPlot <-
    function(ccomp, col = ccomp$iqr, main = NULL, legend.title = "IQR") {
    ## Will use the first 'x_ratio' that appears in the data frame column.
    ratio_name <- grep("*_ratio", colnames(ccomp), value=TRUE)[1]
    df <- data.frame(
        x=ccomp$zscore,
        y=-log10(ccomp[[ratio_name]]),
        col=col,
        row.names=rownames(ccomp)
    )
    ylim <- max(abs(df$y))
    ylim <- c(-ylim, ylim)
    ggplot2::ggplot(
        data=df,
        mapping=ggplot2::aes(
            x=x,
            y=y,
            col=col
        )
    ) + ggplot2::geom_point() +
        ggrepel::geom_text_repel(ggplot2::aes(label=rownames(df))) +
        ggplot2::theme_classic() +
        ggplot2::scale_color_brewer(palette="Dark2") +
        ggplot2::labs(
            x="Heterogeneity (z-score)",
            y=paste0(ratio_name, " -log10"),
            color=legend.title) +
        ggplot2::ggtitle(main) +
        ggplot2::ylim(ylim)
}
