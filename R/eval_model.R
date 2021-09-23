

#' @title Compute metrics for classification model evaluation
#'
#' @param actual vector of actual labels
#'
#' @param predicted vector of predicted labels
#'
#' @return Nothing, called for printing side effect.
#'
#' @export
evaluate_model <- function(actual, predicted) {
    if(length(predicted) != length(actual)) {
        stop(sprintf("Length of actual and predicted labels must match.\n"));
    }
    cm = as.matrix(table(Actual = actual, Predicted = predicted)); # confusion matrix

    n = sum(cm) # number of instances
    nc = nrow(cm) # number of classes
    mdiag = diag(cm) # number of correctly classified instances per class
    rowsums = apply(cm, 1, sum) # number of instances per class
    colsums = apply(cm, 2, sum) # number of predictions per class
    #p = rowsums / n # distribution of instances over the actual classes
    #q = colsums / n # distribution of instances over the predicted classes
    cat(sprintf("Evalutating classification results: %d testing instances assigned to %d classes.\n", length(actual), nc));

    # per class
    accuracy = sum(diag) / n;
    precision = mdiag / colsums;
    recall = diag / rowsums;
    f1 = 2 * precision * recall / (precision + recall);
    print(data.frame(precision, recall, accuracy, f1));

    # macro-averaged
    macroPrecision = mean(precision);
    macroRecall = mean(recall);
    macroF1 = mean(f1);
    print(data.frame(macroPrecision, macroRecall, macroF1));
    return(invisible(NULL));
}
