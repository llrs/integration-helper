#' Plot samples
#'
#' Make several plots of the samples, using several colors and things
#'
#' @param samples A data.frame with RNAseq, Micro, labels columns
#' @param colors Manual colors for some samples
#' @param individual Logical. Do each plot by ID and Time?
#'
#' @return Several plots
#' @export
plot_samples <- function(samples, colors, individual = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 from CRAN", call. = FALSE)
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Install ggforce from CRAN", call. = FALSE)
  }
  # Some common structure of plots
  comm <- ggplot2::ggplot(samples, ggplot2::aes(.data$RNAseq, .data$Micro)) + # It is really biopsies
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::ggtitle("All samples at all times ") +
    ggplot2::xlab("RNAseq (component 1)") +
    ggplot2::ylab("16S (component 1)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (individual) {
    for (p in seq_along(levels(samples$Time))) {
      a <- comm +
        ggplot2::geom_text(ggplot2::aes(color = .data$ID, label = .data$ID)) +
        ggplot2::geom_vline(xintercept = 0) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::guides(col = ggplot2::guide_legend(title = "Patient")) +
        ggplot2::scale_color_manual(values = colors) +
        ggforce::facet_wrap_paginate(~.data$Time, ncol = 1, nrow = 1, page = p)
      print(a)
    }

    for (p in seq_along(levels(samples$ID))) {
      a <- comm +
        ggplot2::geom_text(ggplot2::aes(color = .data$ID,
                      label = ifelse(!is.na(.data$labels),
                                     paste(.data$Time, .data$labels, sep = "_"),
                                     as.character(.data$Time)
        ))) +
        ggplot2:: guides(col = ggplot2::guide_legend(title = "Patient")) +
        ggplot2::scale_color_manual(values = colors) +
        ggforce::facet_wrap_paginate(~.data$ID, ncol = 1, nrow = 1, page = p)
      print(a)
    }
  }
  a <- comm +
    ggplot2::geom_text(ggplot2::aes(
      color = .data$ID,
      label = ifelse(!is.na(.data$labels),
                     paste(.data$Time, .data$labels, sep = "_"),
                     as.character(.data$Time)
      )
    )) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Patient")) +
    ggplot2::scale_color_manual(values = colors)
  print(a)


  a <- comm +
    ggplot2::geom_text(ggplot2::aes(
      color = .data$HSCT_responder,
      label = ifelse(!is.na(.data$labels),
                     paste(.data$ID, .data$labels, sep = "_"),
                     as.character(.data$Time)
      )
    )) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Responders"))
  print(a)


  a <- comm +
    ggplot2::geom_text(ggplot2::aes(
      color = .data$Endoscopic_Activity,
      label = ifelse(!is.na(.data$labels),
                     paste(.data$ID, .data$labels, sep = "_"),
                     as.character(.data$ID)
      )
    )) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Endoscopic Activity"))
  print(a)

  a <- comm +
    ggplot2::geom_text(ggplot2::aes(
      color = .data$Exact_location,
      label = ifelse(!is.na(.data$labels),
                     paste(.data$ID, .data$labels, sep = "_"),
                     as.character(.data$ID)
      )
    )) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Location"))
  print(a)
  a <- comm +
    ggplot2::geom_text(ggplot2::aes(
      color = ifelse(.data$Exact_location == "ILEUM", "ILEUM", "COLON"),
      label = ifelse(!is.na(.data$labels),
                     paste(.data$ID, .data$labels, sep = "_"),
                     as.character(.data$ID)
      )
    )) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Location"))
  print(a)

  a <- comm +
    ggplot2::geom_text(ggplot2::aes(color = .data$Time,
                  label = ifelse(!is.na(.data$labels),
                                 paste(.data$ID, .data$labels, sep = "_"),
                                 as.character(.data$ID)
                  ))) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Time"))
  print(a)

  a <- comm +
    ggplot2::geom_text(ggplot2::aes(color = .data$SESCD_local,
                  label = ifelse(!is.na(.data$labels),
                                 paste(.data$ID, .data$labels, sep = "_"),
                                 as.character(.data$ID)
                  ))) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "SESCD (local)"))
  print(a)

  a <- comm +
    ggplot2::geom_text(ggplot2::aes(color = .data$IBD, label = as.character(.data$ID))) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Type"))
  print(a)
}

#' Plot bullseye
#'
#' Plot the variables of 16S, RNAseq and other origins. In RNAseq and 16S
#' plots variables that in both components are above the mean.
#' @param variables A data.frame with the weight for each variable
#' @return A plot with the most important variables
#' @export
#' @seealso [variables()]
plot_variables <- function(variables) {

  # Remove the variables that in both components are above the mean
  keepComp1RNAseq <- mean(abs(variables$comp1)[variables$Origin == "RNAseq"])
  keepComp1_16S <- mean(abs(variables$comp1)[variables$Origin != "RNAseq"])

  keepComp2RNAseq <- mean(abs(variables$comp2)[variables$Origin == "RNAseq"])
  keepComp2_16S <- mean(abs(variables$comp2)[variables$Origin != "RNAseq"])

  keepComp1 <- c(
    variables$comp1[variables$Origin == "RNAseq"] > keepComp1RNAseq,
    variables$comp1[variables$Origin != "RNAseq"] > keepComp1_16S
  )
  keepComp2 <- c(
    variables$comp2[variables$Origin == "RNAseq"] > keepComp2RNAseq,
    variables$comp2[variables$Origin != "RNAseq"] > keepComp2_16S
  )
  subVariables <- variables[keepComp1 | keepComp2, ]

  a <- ggplot2::ggplot() +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
    ggplot2::geom_text(ggplot2::aes(.data$comp1, .data$comp2, color = .data$Origin,
                  label = .data$var), data = subVariables) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(title = "Variables important for the first two components",
      x = "Comp1", y = "Comp2", col = "Origin") +
    ggplot2::coord_cartesian()
  print(a)
}


#' Plot  PCA of interesting variables
#'
#' @param subVariables A data.frame
#' @param meta_r The data.frame with the information about the samples
#' @param expr A matrix with the RNAseq expression
#' @param otus_table_i A matrix with the 16S
#' @return A
#' @export
plot_interesting <- function(subVariables, meta_r, expr, otus_table_i){

  rnaseq_i <- subVariables$var[subVariables$Origin == "RNAseq"]
  if (length(rnaseq_i) >= 2) {
    pr <- prcomp(t(expr[rnaseq_i, ]), scale. = TRUE)
    prS <- summary(pr)
    a <- ggplot2::ggplot(as.data.frame(pr$x),
                         ggplot2::aes(.data$PC1, .data$PC2, color = as.factor(meta_r$HSCT_responder))) +
      ggplot2::geom_point() +
      ggplot2::xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
      ggplot2::ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
      ggplot2::ggtitle("RNAseq PCA from the important variables") +
      ggplot2::guides(col = ggplot2::guide_legend(title = "Responders"))
    print(a)
  }

  micro_i <- subVariables$var[subVariables$Origin == "16S"]
  if (length(micro_i) >= 2) {
    pr <- prcomp(t(otus_table_i[micro_i, ]), scale. = TRUE)
    prS <- summary(pr)
    a <- ggplot2::ggplot(as.data.frame(pr$x),
                         ggplot2::aes(.data$PC1, .data$PC2,
                    color = as.factor(meta_r$HSCT_responder))) +
      ggplot2::geom_point() +
      ggplot2::xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
      ggplot2::ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
      ggplot2::ggtitle("16S PCA from the important variables") +
      ggplot2::guides(col = ggplot2::guide_legend(title = "Responders"))
    print(a)
  }

}
