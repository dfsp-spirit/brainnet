

#' @title Compute F stats map for mass-univariate GLM analysis.
#'
#' @author C Ecker
#'
#' @inheritParams slm_effect_sizes
#'
#' @importFrom magrittr %>%
#'
#' @keywords internal
slm_F <- function(X, Y, predictors, output = c("F", "p", "p.adjust")) {

    X.full <- X
    term.index <- attr(X, "assign") + 1
    n.terms <- max(term.index)

    n <- dim(X)[1] # number of cases
    p <- vector() # number of model terms

    SS.reg <- matrix(0, nrow = n.terms, ncol = dim(Y)[2])

    for (i in 1:n.terms) {
        X <- X.full[,1:length(term.index[term.index <= i])] # reduce X

        B <- solve(t(X) %*% X) %*% t(X) %*% Y # compute model coefficients for reduced model

        Y.pred <- X %*% B # predicted values
        Y.res <- Y - Y.pred # residuals

        #SS.total <- apply(Y, 2, function(Y) {sum((Y - mean(Y))^2)})
        SS.total <- colSums(t(( t(Y) - colMeans(Y) )^2))
        SS.res <- colSums(Y.res^2)
        SS.reg[i,] <- SS.total - SS.res

        p[i] <- length(which(term.index <= i))
    }

    SS.reg <- SS.reg[-1,] - SS.reg[-dim(SS.reg)[1],]
    rownames(SS.reg) <- predictors
    df.reg <- p[-1] - p[-n.terms]
    MS.reg <- SS.reg / df.reg
    rownames(MS.reg) <- predictors

    df.res <- (n-max(p))
    MS.res <- SS.res / df.res
    F <- apply(MS.reg, 1, function(MS.reg) {MS.reg / MS.res}) %>% t
    p <- stats::pf(F, df.reg, df.res, lower.tail = FALSE)
    p.adjust <- stats::p.adjust(p, method = "fdr") # FDR adjusted p value

    # Mon Jun 15 14:49:55 2020 ------------------------------

    return.obj <- lapply(output, function(x) {get(x, inherits = TRUE)})
    names(return.obj) <- output

    return(return.obj)
}


#' @title Compute t map for mass-univariate GLM analysis.
#'
#' @inheritParams slm_effect_sizes
#'
#' @author C Ecker
#'
#' @importFrom stats pt, p.adjust
#'
#' @keywords internal
slm_t <- function(X, Y, model.term, output=c("t", "p", "p.adjust")) {

    term <- grep(model.term, colnames(X)) # looks for model term index

    n <- dim(X)[1] # number of cases
    p <- dim(X)[2] # number of predictors
    df <- n-p # degrees of freedom

    inv.XtX <- solve(t(X) %*% X)
    b <- inv.XtX %*% t(X) %*% Y

    Y.pred <- X %*% b  # predicted (i.e. fitted) values
    e <- Y - Y.pred # residuals

    e.var <- colSums(e^2/df) # Compute error variance at each vertex
    se <- sqrt(e.var) # Compute residual standard error or MSE summary(model)$sigma at each vertex

    var.b <- e.var * inv.XtX[term, term] # coefficient variance or variance of parameters
    se.b <- sqrt(var.b) # coefficient standard error

    t <- b[term,] / se.b # t value for term
    p <- 2 * stats::pt(t, df) # p-value for t
    p.adjust <- stats::p.adjust(p, method = "fdr") # FDR adjusted p value

    return.obj <- lapply(output, function(x) {get(x, inherits = TRUE)})
    names(return.obj) <- output

    return(return.obj)
}


#' @title Threshold a statistical map by setting all values with p > alpha to NaN.
#'
#' @author C Ecker
#'
#' @param x vector of values
#'
#' @param p the p-values for the x values
#'
#' @param alpha the alpha level
#'
#' @param p.adjust.method passed on to \code{stats::p.adjust}
#'
#' @importFrom stats p.adjust
#'
#' @keywords internal
slm_threshold_statistical_map <- function(x, p, alpha=0.05, p.adjust.method="none") {

    p <- stats::p.adjust(p, method = p.adjust.method)
    x[p > alpha] = NaN

    return(x)
}


#' @title Compute effect sizes for mass-univariate GLM analysis.
#'
#' @author C Ecker
#'
#' @param X numerical matrix, the design or model matrix, typically created from the demographics data using \code{stats::model.matrix}.
#'
#' @param Y numerical matrix, the target value, typically neuroimaging data
#'
#' @param predictors vector of character strings, the names of the predictors in the model matrix
#'
#' @param output vector of pre-defined character strings, defined what values to return. Leave alone if in doubt.
#'
#' @return see output parameter
#'
#' @importFrom magrittr %>%
#' @importFrom pwr pwr.f2.test
#' @importFrom stats pf
#'
#' @keywords internal
slm_effect_sizes <- function(X, Y, predictors, output = c("F", "p", "etasq", "partial.etasq", "cohens.f", "rsq", "power")) {

    X.full <- X
    term.index <- attr(X, "assign") + 1
    n.terms <- max(term.index)

    n <- dim(X)[1] # number of cases
    p <- vector() # number of model terms

    SS.reg <- matrix(0, nrow = n.terms, ncol = dim(Y)[2])

    for (i in 1:n.terms) {
        X <- X.full[,1:length(term.index[term.index <= i])] # reduce X

        B <- solve(t(X) %*% X) %*% t(X) %*% Y # compute model coefficients for reduced model

        Y.pred <- X %*% B # predicted values
        Y.res <- Y - Y.pred # residuals

        #SS.total <- apply(Y, 2, function(Y) {sum((Y - mean(Y))^2)})
        SS.total <- colSums(t(( t(Y) - colMeans(Y) )^2))
        SS.res <- colSums(Y.res^2)
        SS.reg[i,] <- SS.total - SS.res

        p[i] <- length(which(term.index <= i))
    }

    SS.reg <- SS.reg[-1,] - SS.reg[-dim(SS.reg)[1],]
    rownames(SS.reg) <- predictors
    df.reg <- p[-1] - p[-n.terms]
    MS.reg <- SS.reg / df.reg
    rownames(MS.reg) <- predictors

    df.res <- (n-max(p))
    MS.res <- SS.res / df.res
    F <- apply(MS.reg, 1, function(MS.reg) {MS.reg / MS.res}) %>% t
    p <- stats::pf(F, df.reg, df.res, lower.tail = FALSE)


    etasq <- apply(SS.reg, 1, function(SS.reg) {SS.reg / SS.total}) %>% t
    partial.etasq <- apply(SS.reg, 1, function(SS.reg) {SS.reg / (SS.reg + SS.res)}) %>% t
    cohens.f <- sqrt(partial.etasq / (1 - partial.etasq) )
    rsq <- 1 - (SS.res / SS.total)
    power <- sapply(1:dim(Y)[2], function(i) {ifelse(is.na(cohens.f[,i]), NA, pwr::pwr.f2.test(u=df.reg, v=df.res, f2=cohens.f[,i]^2)$power)} )

    return.obj <- lapply(output, function(x) {get(x, inherits = TRUE)})
    names(return.obj) <- output

    return(return.obj)
}


