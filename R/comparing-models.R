
#' Retrieve inner AVE
#'
#' @param x The sgcca or rgcca object
#' @param pos The position of the AVE_inner to extract
#' @return A numeric value
#' @export
getAVEs <- function(x, pos = 1L) {
  x$AVE$AVE_inner[pos]
}



#' Standard error of the mean
#'
#' @param x A numeric vector
#'
#' @return The standard error of the mean
#' @export
sem <- function(x){
  sd(x)/length(x)
}


#' Distribution of inner AVE
#'
#' @param model Original sgcca object
#' @param loo A list with using the same model but using a leave-one-out
#' strategy.
#' @return A histogram with a vertical line indicating the position of the
#' original model
#' @export
plotAVEs <- function(model, loo) {
  aves <- vapply(loo, getAVEs, numeric(1L))
  hist(aves, xlim = c(0, 1), main = "model")
  abline(v = model$AVE$AVE_inner[1])
}


#' Summarize a model
#'
#' @inheritParams plotAVEs
#' @return The inner AVE, the Mean and the SEM of the `loo`
#' @seealso [sem()]
#' @export
#' @importFrom scales scientific
m_sem <- function(model, loo) {
  aves <- vapply(loo, getAVEs, numeric(1L))
  paste0(signif(model$AVE$AVE_inner[1], 3),
         " (", signif(mean(aves), 3),
         " \u00b1 ", scales::scientific(sem(aves), 3), ")")
}


#' Clean the output of a sgcca object
#'
#' @param data The matrices of `sggcca$Y[[1]]` or `sgcca$a[[1]]`
#' @param model A label of the model used
#' @param type A label to know what it
#'
#' @return A `data.frame`
#' @export
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
tidyer <- function(data, model, type) {
  if ("comp1" %in% colnames(data)) {
    if ("comp2" %in% colnames(data)) {
        m <- dplyr::mutate(as.data.frame(data), Model = model)
        d <- tidyr::gather(data = m, "Component", !!type, .data$comp1:.data$comp2)
    } else {

      m <- dplyr::mutate(as.data.frame(data), Model = model)
      d <- tidyr::gather(data = m, "Component", !!type, .data$comp1)

    }
  } else {
    m <- dplyr::mutate(as.data.frame(data), Model = model)
    d <- tidyr::gather(data = m, "Component", !!type, 1:2)
  }
  if (!is.null(rownames(data))) {
    d$Rownames <- rep(rownames(data), nrow(d)/nrow(data))
  }
  d
  # mutate(Sample = seq_len(n()))
  # Samples name could be important!!
}
