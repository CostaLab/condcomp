library(condcomp)
context("Testing 'condcomp' call")

test_that("Check the type of the data frame columns", {
    x <- iris[-length(iris)]
    clustering <- iris$Species
    dmatrix <- as.matrix(dist(x))
    cond <- sample(1:2, nrow(iris), replace = TRUE)
    comp <- condcomp(clustering,
                     cond,
                     dmatrix = dmatrix,
                     n = 10)
    comp_numeric <- comp[colnames(comp) != "iqr"]
    expect_equal(all(sapply(comp_numeric, is.numeric)), TRUE)
    comp_factor <- comp["iqr"]
    expect_equal(all(sapply(comp_factor, is.factor)), TRUE)
})

test_that("Remove one of the levels of the clustering factor completely", {
    iris_sub <- subset(iris, iris$Species != "setosa")
    x <- iris_sub[-length(iris_sub)]
    clustering <- iris_sub$Species
    dmatrix <- as.matrix(dist(x))
    cond <- sample(1:2, nrow(iris_sub), replace = TRUE)
    comp <- condcomp(clustering, cond, dmatrix = dmatrix, n = 10)
    expect_equal(is.data.frame(comp), TRUE)
})

test_that("One of the clusters only has data from one condition (returning NA)",
          {
              x <- iris[-length(iris)]
              clustering <- iris$Species
              dmatrix <- as.matrix(dist(x))
              cond <- sample(1:2, nrow(iris), replace = TRUE)
              cond[iris$Species == "setosa"] <- 2
              comp <- condcomp(
                  clustering,
                  cond,
                  dmatrix = dmatrix,
                  n = 10,
                  remove.na = FALSE
              )
              expect_equal(is.data.frame(comp), TRUE)
              setosa <- comp["setosa", ]
              expect_equal(all(is.na(setosa[c("true_sil", "zscore", "pval", "iqr")])), TRUE)
              expect_equal(sum(setosa[c("1_perc", "1_ratio")]), 0)
              expect_equal(setosa[["2_perc"]], 1)
              expect_equal(setosa[["2_ratio"]], Inf)
          })

test_that("One of the clusters only has data from one condition (omitting NA)",
          {
              x <- iris[-length(iris)]
              clustering <- iris$Species
              dmatrix <- as.matrix(dist(x))
              cond <- sample(1:2, nrow(iris), replace = TRUE)
              cond[iris$Species == "setosa"] <- 2
              comp <- condcomp(
                  clustering,
                  cond,
                  dmatrix = dmatrix,
                  n = 10,
                  remove.na = TRUE
              )
              expect_equal(is.data.frame(comp), TRUE)
              expect_equal(grep("setosa", rownames(comp)), integer(0))
          })

test_that("'clustering' parameter is not factor", {
    x <- iris[-length(iris)]
    clustering <- as.character(iris$Species)
    dmatrix <- as.matrix(dist(x))
    cond <- sample(1:2, nrow(iris), replace = TRUE)
    comp <- condcomp(
        clustering,
        cond,
        dmatrix = dmatrix,
        n = 10,
        remove.na = TRUE
    )
    expect_equal(is.data.frame(comp), TRUE)
    expect_equal(any(is.na(comp)), FALSE)
})

test_that("'dmatrix' parameter is not a matrix", {
    x <- iris[-length(iris)]
    clustering <- as.character(iris$Species)
    dmatrix <- dist(x)
    cond <- sample(1:2, nrow(iris), replace = TRUE)
    comp <- condcomp(
        clustering,
        cond,
        dmatrix = dmatrix,
        n = 10,
        remove.na = TRUE
    )
    expect_equal(is.data.frame(comp), TRUE)
    expect_equal(any(is.na(comp)), FALSE)
})
