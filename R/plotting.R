#' @importFrom ggforce facet_wrap_paginate
NULL

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

  # Some common structure of plots
  comm <- ggplot(samples, aes(RNAseq, Micro)) + # It is really biopsies
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle("All samples at all times ") +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    theme(plot.title = element_text(hjust = 0.5))
  if (individual) {
    for (p in seq_along(levels(samples$Time))) {
      a <- comm +
        geom_text(aes(color = ID, label = ID)) +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        guides(col = guide_legend(title = "Patient")) +
        scale_color_manual(values = colors) +
        ggforce::facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
      print(a)
    }

    for (p in seq_along(levels(samples$ID))) {
      a <- comm +
        geom_text(aes(color = ID, label = ifelse(!is.na(labels),
                                                 paste(Time, labels, sep = "_"),
                                                 as.character(Time)
        ))) +
        guides(col = guide_legend(title = "Patient")) +
        scale_color_manual(values = colors) +
        ggforce::facet_wrap_paginate(~ID, ncol = 1, nrow = 1, page = p)
      print(a)
    }
  }
  a <- comm +
    geom_text(aes(
      color = ID,
      label = ifelse(!is.na(labels),
                     paste(Time, labels, sep = "_"),
                     as.character(Time)
      )
    )) +
    guides(col = guide_legend(title = "Patient")) +
    scale_color_manual(values = colors)
  print(a)


  a <- comm +
    geom_text(aes(
      color = HSCT_responder,
      label = ifelse(!is.na(labels),
                     paste(ID, labels, sep = "_"),
                     as.character(Time)
      )
    )) +
    guides(col = guide_legend(title = "Responders"))
  print(a)


  a <- comm +
    geom_text(aes(
      color = Endoscopic_Activity,
      label = ifelse(!is.na(labels),
                     paste(ID, labels, sep = "_"),
                     as.character(ID)
      )
    )) +
    guides(col = guide_legend(title = "Endoscopic Activity"))
  print(a)

  a <- comm +
    geom_text(aes(
      color = Exact_location,
      label = ifelse(!is.na(labels),
                     paste(ID, labels, sep = "_"),
                     as.character(ID)
      )
    )) +
    guides(col = guide_legend(title = "Location"))
  print(a)
  a <- comm +
    geom_text(aes(
      color = ifelse(Exact_location == "ILEUM", "ILEUM", "COLON"),
      label = ifelse(!is.na(labels),
                     paste(ID, labels, sep = "_"),
                     as.character(ID)
      )
    )) +
    guides(col = guide_legend(title = "Location"))
  print(a)

  a <- comm +
    geom_text(aes(color = Time,
                  label = ifelse(!is.na(labels),
                                 paste(ID, labels, sep = "_"),
                                 as.character(ID)
                  ))) +
    guides(col = guide_legend(title = "Time"))
  print(a)

  a <- comm +
    geom_text(aes(color = SESCD_local,
                  label = ifelse(!is.na(labels),
                                 paste(ID, labels, sep = "_"),
                                 as.character(ID)
                  ))) +
    guides(col = guide_legend(title = "SESCD (local)"))
  print(a)

  a <- comm +
    geom_text(aes(color = IBD, label = as.character(ID))) +
    guides(col = guide_legend(title = "Type"))
  print(a)
}

#' Plot bullseye
#'
#' Plot the variables of 16S, RNAseq and other origins. In RNAseq and 16S
#' plots variables that in both components are above the mean.
#' @param variables A data.frame with the weight for each variable
#' @return A plot with the most important variables
#' @export
#' @seealso \code{\link{variables}}
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

  a <- ggplot(subVariables, aes(comp1, comp2), color = Origin) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
    geom_text(aes(color = Origin, label = var)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    coord_cartesian() +
    ggtitle(
      "Variables important for the first two components",
      subtitle = "Integrating stools and mucosa samples"
    )
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
    a <- ggplot(as.data.frame(pr$x),
                aes(PC1, PC2, color = as.factor(meta_r$HSCT_responder))) +
      geom_point() +
      xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
      ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
      ggtitle("RNAseq PCA from the important variables") +
      guides(col = guide_legend(title = "Responders"))
    print(a)
  }

  micro_i <- subVariables$var[subVariables$Origin == "16S"]
  if (length(micro_i) >= 2) {
    pr <- prcomp(t(otus_table_i[micro_i, ]), scale. = TRUE)
    prS <- summary(pr)
    a <- ggplot(as.data.frame(pr$x),
                aes(PC1, PC2, color = as.factor(meta_r$HSCT_responder))) +
      geom_point() +
      xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
      ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
      ggtitle("16S PCA from the important variables") +
      guides(col = guide_legend(title = "Responders"))
    print(a)
  }

}
