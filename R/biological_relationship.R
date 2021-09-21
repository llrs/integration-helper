
#' Performs the calculation of biological information
#'
#' Looks if the score is really associated in the permutations
#' Does the ORA enrichment and fgsea analysis on the genes
#' @param otus_tax matrix as output of `taxonomy` function
#' @param sgcca.centroid SGCCA output
#' @param STAB Output of bootstrap by `boot_sgcca` function.
#' @param label Name of the output files
#' @param epithelium Data from a file from the lab
#' @param today Date as in character format
#' @export
biological_relationships <- function(sgcca.centroid, STAB, label, otus_tax,
                                     epithelium, today) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 from CRAN", call. = FALSE)
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Install clusterProfiler from Bioconductor", call. = FALSE)
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install org.Hs.eg.db from Bioconductor", call. = FALSE)
  }
  if (!requireNamespace("reactome.db", quietly = TRUE)) {
    stop("Install reactome.db from Bioconductor", call. = FALSE)
  }
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Install fgsea from Bioconductor", call. = FALSE)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Install AnnotationDbi from Bioconductor", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table from CRAN", call. = FALSE)
  }

  # RNAseq ####
  b <- STAB[["RNAseq"]]
  d <- sgcca.centroid$a[["RNAseq"]][, 1]

  # Remove duplicated if sgcca failed due to LAPACK subroutine
  boot_NA <- is.na(b[, 1])
  b <- b[!boot_NA, , drop = FALSE]
  d_rm <- d[!boot_NA]

  if (sum(boot_NA) > 1) {
    warning(sum(boot_NA), " iterations failed.")
  }
  if (nrow(b) == 0) {
    stop("Not enough data to continue")
  }

  pdf(paste0("Figures/", today, "_enrichments_", label, ".pdf"))

  # Test if the gene is significant by comparing to how much times is different
  # from 0 (because CCA tends to compensate itself)
  count <- apply(b, 2, function(x) sum(x != 0, na.rm = TRUE))
  freq <- count / nrow(b)

  # Filter to match length
  if (length(d_rm) > length(freq)) {
    freq <- freq[names(d_rm)]
  } else if (length(freq) > length(d_rm)) {
    d_rm <- d_rm[names(freq)]
  }

  p <- ggplot2::ggplot(as.data.frame(cbind(freq, d_rm))) +
    ggplot2::geom_point(ggplot2::aes(d_rm, freq)) +
    ggplot2::xlab("Score") +
    ggplot2::ylab("Frequency (!= 0)") +
    ggplot2::ggtitle("RNAseq")
  print(p)

  # Select those genes that are significant
  # Selected more than 50% of the time
  significant <- names(d_rm)[freq > 0.5]
  significant <- trimVer(significant)

  loadings <- sgcca.centroid$a[["RNAseq"]]

  # Convert the ids to the entrezIds
  ensemblID <- rownames(loadings)
  ensemblID <- trimVer(ensemblID)
  rownames(loadings) <- trimVer(rownames(loadings))
  entrezID <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensemblID, keytype = "ENSEMBL",
    column = "ENTREZID"
  )
  comp1 <- loadings[, 1]
  names(comp1) <- entrezID
  comp1 <- comp1[diff0(comp1)]

  # Extract the information of the pathways
  paths2genes <- access_reactome()
  genes <- unlist(paths2genes, use.names = FALSE)
  pathways <- rep(names(paths2genes), lengths(paths2genes))

  ## Compute the hypergeometric/enrichment analysis ####
  message("Calculating the enrichment")
  T2G <- cbind.data.frame(pathways, genes)
  enrich <- clusterProfiler::enricher(
    gene = entrezID[significant], TERM2GENE = T2G
  )
  enrich <- as.data.frame(enrich)
  if (nrow(enrich) >= 1) {
    enrich$Description <- AnnotationDbi::mapIds(reactome.db::reactome.db, keys = rownames(enrich),
                                 keytype = "PATHID", column = "PATHNAME")
  }
  write.csv(enrich,
    file = paste0("enrichment_RNAseq_", label, ".csv"), row.names = FALSE
  )

  # Store the entrezid
  entrezSig <- entrezID[significant]
  entrezSig <- entrezSig[!is.na(entrezSig)]
  paths2genes[["significant"]] <- entrezSig
  paths2genes[["Epithelium"]] <- epithelium

  message("Calculating the gene set enrichment")
  ## Compute the GSEA for the size effect ####
  gseaSizeEffect <- fgsea::fgsea(paths2genes, comp1, nperm = length(comp1))

  # Get the name of the pathway
  namesPaths <- AnnotationDbi::select(
    reactome.db::reactome.db,
    keys = gseaSizeEffect$pathway,
    keytype = "PATHID", columns = "PATHNAME"
  )
  NES <- padj <- pval <- 0
  # Remove the homo sapiens part
  namesPaths$PATHNAME <- gsub("Homo sapiens: (.*)", "\\1", namesPaths$PATHNAME)
  # Add a column
  gseaSizeEffect[, namesPaths := namesPaths$PATHNAME]
  # Order the dataframe by size effect
  gseaSizeEffect <- gseaSizeEffect[order(-abs(NES), padj, -size)]
  if (sum(gseaSizeEffect$padj < 0.05) == 0) {
    warning("GSEA didn't result in any pathway")
  }
  # Store the output
  data.table::fwrite(gseaSizeEffect[pval < 0.05, ],
    file = paste0("gsea_pathways_", label, ".csv")
  )

  ## 16S ####
  b <- STAB[["16S"]]
  d <- sgcca.centroid$a[["16S"]][, 1]

  # Remove if sgcca failed due to LAPACK subroutine
  b <- b[!boot_NA, ]
  d_rm <- d[!boot_NA]

  # Test if the gene is significant by comparing to how much times is different
  # from 0 (because CCA tends to compensate itself)
  count <- apply(b, 2, function(x) sum(x != 0, na.rm = TRUE))
  freq <- count / nrow(b)

  # Filter to match length
  if (length(d_rm) > length(freq)) {
    freq <- freq[names(d_rm)]
  } else if (length(freq) > length(d_rm)) {
    d_rm <- d_rm[names(freq)]
  }

  p <- ggplot2::ggplot(as.data.frame(cbind(freq, d_rm))) +
    ggplot2::geom_point(ggplot2::aes(d_rm, freq)) +
    ggplot2::xlab("Score") +
    ggplot2::ylab("Frequency (!= 0)") +
    ggplot2::ggtitle("16S")
  print(p)
  otus <- names(d_rm)[freq > 0.5]


  ## Split the taxonomy
  Taxon2Class <- as.list(as.data.frame(otus_tax))

  grouping <- split(Taxon2Class$Genus, Taxon2Class$Genus)
  grouping <- sapply(grouping, names)

  term2gene <- data.frame(
    "Gene" = otus_tax[, "Genus"],
    "Term" = rownames(otus_tax)
  )
  term2name <- data.frame(
    "Name" = otus_tax[, "Genus"],
    "Term" = rownames(otus_tax)
  )
  enrich <- as.data.frame(clusterProfiler::enricher(
    gene = otus, universe = rownames(otus_tax),
    minGSSize = 1, TERM2GENE = term2gene,
    TERM2NAME = term2name
  ))

  write.csv(enrich, file = paste0("enrichment_otus_genus_", label, ".csv"))

  # GSEA
  comp1 <- sgcca.centroid$a[["16S"]][, 1]
  comp1 <- comp1[diff0(comp1)]
  gseaSizeEffect <- fgsea::fgsea(grouping, comp1, nperm = 20000)
  gseaSizeEffect <- gseaSizeEffect[order(-abs(NES), padj, -size)]
  if (sum(gseaSizeEffect$padj < 0.05) == 0) {
    warning("GSEA didn't result in any pathway")
  }
  data.table::fwrite(gseaSizeEffect[pval < 0.05],
    file = paste0("gsea_otus_genus_", label, ".csv")
  )
  dev.off()
}

#' Return data from reactome
#'
#'
#' Looks in reactome for the pathways associated with each gene
#' @return A list of paths and human genes in the list
#' @export
#' @examples
#' paths2genes <- access_reactome()
access_reactome <- function(){
  if (!requireNamespace("reactome.db", quietly = TRUE)) {
    stop("Install reactome.db from Bioconductor", call. = FALSE)
  }
  genes2Pathways <- as.list(reactome.db::reactomeEXTID2PATHID)
  pathways <- unlist(genes2Pathways, use.names = FALSE)
  genes <- rep(names(genes2Pathways), lengths(genes2Pathways))
  paths2genes <- split(genes, pathways)
  human <- grep("R-HSA-", names(paths2genes))
  paths2genes <- paths2genes[human]
  paths2genes[lengths(paths2genes) >= 2]
}
