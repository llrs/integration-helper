#' Integrate
#'
#' integrates and does the bootstrap of the data provided. It creates also the
#' pdfs and the files with data
#' @param A The list with the three datsets to integrate
#' @param meta The metadata
#' @param label A label to apply to the files
#' @param today A date in YYYYMMDD format in character
#' @return The sgcca output
#' @export
integration <- function(A, meta, label, today) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 from CRAN", call. = FALSE)
  }
  # We cannnot comput eht tau.estimate for A[[1]]
  # (shrinkage <- sapply(A, tau.estimate))
  shrinkage <- c(0.25670333, 0, 1) # We guess a 0.1 for the RNAseq expression
  shrinkage[2] <- tau.estimate(A[[2]])
  (min_shrinkage <- sapply(A, function(x) {
    1 / sqrt(ncol(x))
  }))
  # # Don't let the shrinkage go below the thershold allowed
  shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
  # shrinkage <- rep(1, length(A))

  ncomp <- c(2, 2, 2)

  sgcca.centroid <- sgcca(
    A, C,
    c1 = shrinkage,
    ncomp = ncomp,
    scheme = "centroid",
    scale = TRUE,
    verbose = FALSE
  )
  names(sgcca.centroid$Y) <- names(A)
  names(sgcca.centroid$a) <- names(A)
  names(sgcca.centroid$astar) <- names(A)
  names(sgcca.centroid$AVE$AVE_X) <- names(A)
  sgcca.centroid$AVE$AVE_X <- simplify2array(sgcca.centroid$AVE$AVE_X)
  sgcca.centroid$AVE

  # list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
  # sgcca.factorial = sgcca.factorial)
  save(sgcca.centroid, file = paste0("sgcca_", label, ".RData"))

  samples <- data.frame(
    "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
    "microbiota" = sgcca.centroid$Y[["16S"]][, 1],
    "metadata" = sgcca.centroid$Y[["metadata"]][, 1]
  )


  ## Grouping of the variables ####
  RNAseq1 <- samples$RNAseq
  RNAseq2 <- sgcca.centroid$Y[["RNAseq"]][, 2]
  microbiota2 <- sgcca.centroid$Y[["16S"]][, 2]
  microbiota1 <- samples$microbiota

  names(RNAseq1) <- rownames(samples)
  names(microbiota1) <- rownames(samples)
  names(RNAseq2) <- rownames(samples)
  names(microbiota2) <- rownames(samples)
  groups <- split(rownames(samples), as.factor(meta$HSCT_responder))
  # First dimension seems to capture well the
  fgsea::fgsea(groups, RNAseq1, nperm = 1000)
  fgsea::fgsea(groups, microbiota1, nperm = 1000)
  # Further dimensions
  fgsea::fgsea(groups, RNAseq2, nperm = 1000)
  fgsea::fgsea(groups, microbiota2, nperm = 1000)


  grDevices::pdf(paste0("Figures/", today, "_RGCCA_plots_", label, ".pdf"))

  km <- kmeans(samples[, c("RNAseq", "microbiota")], 2, nstart = 2)
  plot(samples[, c("RNAseq", "microbiota")],
    col = km$cluster,
    main = "K-clustering (2 groups)"
  )

  ## Plotting ####
  # Colors for the plots
  names(colors) <- unique(meta$ID)

  samples <- cbind(samples, droplevels(meta))
  samples$Patient_ID <- as.factor(samples$Patient_ID)
  samples$Sample_Code <- as.character(samples$Sample_Code)

  # Labels of the samples
  label <- strsplit(as.character(samples$`Sample Name_RNA`), split = "-")
  labels <- sapply(label, function(x) {
    if (length(x) == 5) {
      x[5]
    }
    else if (length(x) != 5) {
      x[4]
    }
  })

  samples <- cbind(samples, labels)
  samples$Time <- factor(
    samples$Time,
    levels(as.factor(samples$Time))[c(1, 2, 4, 5, 3, 6, 7, 8)]
  )
  for (p in seq_along(levels(samples$Time))) {
    a <- ggplot2::ggplot(samples, ggplot2::aes(.data$RNAseq, .data$microbiota)) +
      ggplot2::geom_text(ggplot2::aes(color = .data$ID, label = .data$ID)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::ggtitle(paste0("Samples by time")) +
      ggplot2::xlab("RNAseq (component 1)") +
      ggplot2::ylab("16S (component 1)") +
      ggplot2::guides(col = ggplot2::guide_legend(title = "Patient")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_color_manual(values = colors) +
      ggforce::facet_wrap_paginate(~ .data$Time, ncol = 1, nrow = 1, page = p)
    print(a)
  }
  for (p in seq_along(levels(samples$ID))) {
    a <- ggplot2::ggplot(samples, ggplot2::aes(.data$RNAseq, .data$microbiota)) +
      ggplot2::geom_text(ggplot2::aes(color = .data$ID,
                                      label = ifelse(!is.na(.data$labels),
                                                     paste(.data$Time, .data$labels, sep = "_"),
                                                     as.character(.data$Time)
                                      ))) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::ggtitle(paste0("Samples by patient")) +
      ggplot2::xlab("RNAseq (component 1)") +
      ggplot2::ylab("16S (component 1)") +
      ggplot2::guides(col = ggplot2::guide_legend(title = "Patient")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_color_manual(values = colors) +
      ggforce::facet_wrap_paginate(~ ID, ncol = 1, nrow = 1, page = p)
    print(a)
  }
  ggplot2::ggplot(samples, ggplot2::aes(.data$RNAseq, .data$microbiota)) +
    ggplot2::geom_text(ggplot2::aes(
      color = .data$ID,
      label = ifelse(!is.na(.data$labels),
        paste(.data$Time, .data$labels, sep = "_"),
        as.character(.data$Time)
      )
    )) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::ggtitle("All samples at all times ") +
    # xlab("RNAseq (component 1)") +
    # ylab("16S (component 1)") +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Patient")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_color_manual(values = colors)


  ggplot2::ggplot(samples, ggplot2::aes(.data$RNAseq, .data$microbiota)) +
    ggplot2::geom_text(ggplot2::aes(
      color = .data$HSCT_responder,
      label = ifelse(!is.na(.data$labels),
        paste(.data$Time, .data$labels, sep = "_"),
        as.character(.data$Time)
      )
    )) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::ggtitle("All samples at all times ") +
    ggplot2::xlab("RNAseq (component 1)") +
    ggplot2::ylab("16S (component 1)") +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Responders")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


  ggplot2::ggplot(samples, ggplot2::aes(.data$RNAseq, .data$microbiota)) +
    ggplot2::geom_text(ggplot2::aes(
      color = .data$Endoscopic_Activity,
      label = ifelse(!is.na(.data$labels),
        paste(.data$ID, .data$labels, sep = "_"),
        as.character(.data$ID)
      )
    )) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::ggtitle("All samples at all times ") +
    ggplot2::xlab("RNAseq (component 1)") +
    ggplot2::ylab("16S (component 1)") +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Endoscopic Activity")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggplot(samples, ggplot2::aes(.data$RNAseq, .data$microbiota)) +
    ggplot2::geom_text(ggplot2::aes(color = .data$Time, label = ifelse(!is.na(.data$labels),
      paste(.data$ID, .data$labels, sep = "_"),
      as.character(.data$ID)
    ))) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::ggtitle("All samples at all times ") +
    ggplot2::xlab("RNAseq (component 1)") +
    ggplot2::ylab("16S (component 1)") +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Time")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  subVariables <- select_var(sgcca.centroid)

  ggplot2::ggplot(subVariables, ggplot2::aes(.data$comp1, .data$comp2, color = .data$Origin)) +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
    ggplot2::geom_path(ggplot2::aes(.data$x, .data$y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
    ggplot2::geom_text(ggplot2::aes(color = .data$Origin, label = .data$var)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::coord_cartesian() +
    ggplot2::ggtitle(
      "Variables important for the first two components",
      subtitle = "Integrating stools and mucosa samples"
    )

  # Plot for the same component the variables of each block
  comp1 <- sapply(sgcca.centroid$a, function(x) {
    x[, 1]
  })
  variables_weight(comp1)

  # Second component
  comp2 <- sapply(sgcca.centroid$a, function(x) {
    x[, 2]
  })
  variables_weight(comp2)

  # Bootstrap of sgcca
  STAB <- boot_sgcca(A, C, shrinkage, 1000)

  save(STAB, file = paste("bootstrap_", label, ".RData"))

  # Evaluate the boostrap effect and plot
  boot_evaluate(STAB)

  PCAs_important(subVariables, A[["RNASeq"]], A[["16S"]], meta)

  dev.off()
  sgcca.centroid
}

PCAs_important <- function(important_vars, expr, otus, meta) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr from CRAN", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 from CRAN", call. = FALSE)
  }
  subVariables <- important_vars
  rnaseq_i <- subVariables$var[subVariables$Origin == "RNAseq"]
  if (length(rnaseq_i) >= 2) {
    pr <- prcomp(t(expr[rnaseq_i, ]), scale. = TRUE)
    prS <- summary(pr)
    p <- ggplot2::ggplot(as.data.frame(pr$x),
                         ggplot2::aes(.data$PC1, .data$PC2,
                    color = as.factor(meta$HSCT_responder))) +
      ggplot2::geom_point() +
      ggplot2::xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
      ggplot2::ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
      ggplot2::ggtitle("RNAseq PCA from the important variables") +
      ggplot2::guides(col = ggplot2::guide_legend(title = "Responders"))
    print(p)
  }

  micro_i <- subVariables$var[subVariables$Origin == "16S"]
  if (length(micro_i) >= 2) {
    pr <- prcomp(t(otus[micro_i, ]), scale. = TRUE)
    prS <- summary(pr)
    p <- ggplot2::ggplot(as.data.frame(pr$x),
                         ggplot2::aes(.data$PC1, .data$PC2,
                    color = as.factor(meta$HSCT_responder))) +
      ggplot2::geom_point() +
      ggplot2::xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
      ggplot2::ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
      ggplot2::ggtitle("16S PCA from the important variables") +
      ggplot2::guides(col = ggplot2::guide_legend(title = "Responders"))
    print(p)
  }
}

#' Select important variables
#' @param sgcca result of sgcca or rgcca
#' @return a data.frame with variables and their location
#' @export
select_var <- function(sgcca) {
  variables <- data.frame(
    Origin = rep(names(sgcca$Y), sapply(sgcca$a, nrow)),
    comp1 = unlist(sapply(
      sgcca$a,
      function(x) {
        x[, 1]
      }
    )),
    comp2 = unlist(sapply(
      sgcca$a,
      function(x) {
        x[, 2]
      }
    ))
  )
  variables$var <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
  rownames(variables) <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
  variables$var <- gsub("^RNAseq\\.(ENSG.*)$", "\\1", rownames(variables))
  rownames(variables) <- gsub("^.*\\.(ENSG.*)$", "\\1", rownames(variables))
  rownames(variables) <- gsub("^metadata\\.(.*)$", "\\1", rownames(variables))
  variables$var <- gsub("^metadata\\.(.*)$", "\\1", rownames(variables))

  # Remove the variables that in both components are 0
  keepComp1RNAseq <- mean(abs(variables$comp1)[variables$Origin == "RNAseq"])
  keepComp1_16S <- mean(abs(variables$comp1)[variables$Origin == "16S"])
  keepComp1_metadata <- mean(abs(variables$comp1)[variables$Origin == "metadata"])

  keepComp2RNAseq <- mean(abs(variables$comp2)[variables$Origin == "RNAseq"])
  keepComp2_16S <- mean(abs(variables$comp2)[variables$Origin == "16S"])
  keepComp2_metadata <- mean(abs(variables$comp2)[variables$Origin == "metadata"])

  keepComp1 <- c(
    abs(variables$comp1[variables$Origin == "RNAseq"]) > keepComp1RNAseq,
    abs(variables$comp1[variables$Origin == "16S"]) > keepComp1_16S,
    abs(variables$comp1[variables$Origin == "metadata"]) > keepComp1_metadata
  )
  keepComp2 <- c(
    abs(variables$comp2[variables$Origin == "RNAseq"]) > keepComp2RNAseq,
    abs(variables$comp2[variables$Origin == "16S"]) > keepComp2_16S,
    abs(variables$comp2[variables$Origin == "metadata"]) > keepComp2_metadata
  )

  variables[keepComp1 & keepComp2, ]
}
