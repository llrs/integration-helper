#' Enrichment by microorganisms
#'
#' Function to split by microorganism and compare if they are enriched in some
#' pathways
#' @param all_genes All the genes that are found relevant (weight != 0)
#' @param x The table of microoganisms, genes, correlations and adjusted p-values.
#' @return A list of tables with the pathways enriched for every microorganism
#' @export
pathsPerMicro <- function(x, all_genes){
  genesOrig <- trimVer(names(all_genes))
  x[, "Gene"] <- trimVer(x[, "Gene"])
  perMicro <- split(x, x[, "Microorganism"])
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  requireNamespace("reactome.db", quietly = TRUE)
  entrezID <- mapIds(org.Hs.eg.db, keys = genesOrig, keytype = "ENSEMBL",
                     column = "ENTREZID")
  entrezID <- entrezID[!is.na(entrezID)]
  paths2genes <- access_reactome()
  genes <- unlist(paths2genes, use.names = FALSE)
  pathways <- rep(names(paths2genes), lengths(paths2genes))
  message("Calculating the enrichment")
  T2G <- cbind.data.frame(pathways, genes)
  lapply(perMicro, function(x) {
    significant <- x$Gene
    relevant <- entrezID[significant]
    relevant <- relevant[!is.na(relevant)]
    if (length(relevant) <= 10) {
      return(NULL)
    }
    enrich <- clusterProfiler::enricher(gene = relevant ,
                                        # universe = entrezID,
                                        TERM2GENE = T2G)
    as.data.frame(enrich)
    enrich <- as.data.frame(enrich)
    requireNamespace("reactome.db", quietly = TRUE)
    if (nrow(enrich) >= 1) {
      enrich$Description <- mapIds(reactome.db, keys = rownames(enrich),
                                   keytype = "PATHID", column = "PATHNAME")
    }
    enrich
  })
}
#' Plot correlations
#'
#' @param x Gene expression matrix
#' @param gene Name of the gene to plot
#' @param y Microorganism matrix
#' @param microorganism Name of the microorganism
#' @param colr Factor to use for colors
#' @param case Factor to use for points shape
#' @param cor_val Title of the correlation (usually the correlation value)
#' @export
#' @importFrom graphics legend
plot_single_cor <- function(x, gene, y, microorganism, colr, case, cor_val) {
  genes <- trimVer(gene)
  symbol <- tryCatch({mapIds(
    org.Hs.eg.db, keys = genes, keytype = "ENSEMBL",
    column = "SYMBOL"
  )},  error = function(e){NA})

  if (is.na(symbol)) {
    return(NA)
  }

  x_s <- x[gene, ]
  x_s[x_s == 0] <- NA

  y_s <- y[microorganism, ]
  y_s[y_s == 0] <- NA
  pch_start <- 15

  # main <- cor(x_s, y_s, method = "spearman", use = "pairwise.complete.obs")
  plot(x_s, y_s, xlab = paste(symbol, collapse = " "), ylab = microorganism, main = cor_val, col = colr,
       pch = pch_start + as.numeric(case))
  legend("bottomleft", fill = as.factor(levels(colr)), legend = levels(colr))
  legend("topright", pch = pch_start + as.numeric(as.factor(levels(case))),
         legend = levels(case))
  # legend("topleft", lty = "solid", col = as.factor(levels(bg)),
  # legend = levels(bg), lwd = 2)
}


#' Bootstrapping to asses how probable is to have such enrichment in those p-values matrix
#' @param pvalues Vector of pvalues
#' @param component The named component different from 0 from a SGCCA object
#' @param iter The number of permutations
#' @return The numeric pvalue of the bootstrap.
#' @export
boots_corr <- function(pvalues, component, iter = 10000) {
  b <- vapply(seq_len(iter), function(x){
    sel <- sample(x = seq_len(ncol(pvalues)), size = length(component))
    sum(pvalues[, sel] < 0.05, na.rm = TRUE)
  }, numeric(1L))
  n <- sum(pvalues[, colnames(pvalues) %in% names(component)] < 0.05, na.rm = TRUE)
  sum(b >= n)/length(b)
}

#' Read component
#'
#' @param file The .RDS file were the sgcca object is stored
#' @param set The set for which we want the component
#' @param component The number of the component we want.
#' @return The vector of the component without those variables with weight 0
#' @export
readSGCCA <- function(file, set = 1, component = 1) {
  # Load data from all the patients
  sgcca.centroid <- readRDS(file)

  # Find outliers/important genes
  comp1 <- sgcca.centroid$a[[set]][, component]
  outliers <- comp1 != 0
  comp1[outliers]
}

#' Filter correlations
#'
#' Filter the correlations based on the pvalue
#' @param comp1 A component with names
#' @param cors A matrix with the correlations
#' @param pval A matrix with the pvalues
#' @param threshold The numeric threshold of the filter
#' @return A list with matrices of correlations and pvalues that pass the filter
#' @export
#' @seealso [sign_cor()], [relevant()]
filter_values <- function(comp1, cors, pval, threshold) {

  keepGenes <- colnames(cors) %in% names(comp1)
  cors <- cors[, keepGenes]
  pval <- pval[, keepGenes]

  cors <- cors[, !is.na(colnames(cors))]
  pval <- pval[, !is.na(colnames(pval))]

  message("Dimensions ", paste0(dim(cors), collapse = ", "))

  # Genes below the threshold
  keepCols <- apply(pval, 2, function(x){any(x < threshold)})
  keepRows <- apply(pval, 1, function(x){any(x < threshold)})

  keepCols[is.na(keepCols)] <- FALSE
  keepRows[is.na(keepRows)] <- FALSE

  if (sum(keepCols) == 0 || sum(keepRows) == 0) {
    stop("No relevant correlations with this threshold")
  }

  cors <- cors[keepRows, keepCols]
  pval <- pval[keepRows, keepCols]
  if (is.null(pval) || is.null(cors)) {
    stop("No relevant correlations with this threshold")
  }
  message("Dimensions ", paste0(dim(cors), collapse = ", "))

  list(cors = cors, pval = pval)
}



#' List the correlations
#'
#' Prepare the correlations to be shared
#' @inheritParams filter_values
#' @return A data.frame with microorganisms, genes, their correlations and the
#' pvalue
#' @export
#' @seealso [filter_values()], [sign_cor()]
relevant <- function(comp1, cors, pval, threshold = 0.05) {
  l <- filter_values(comp1, cors, pval, threshold)
  pval <- l$pval
  cors <- l$cors
  if (ncol(pval) == 0) {
    stop("No relevant correlations with this threshold")
  }
  ind <- as.data.frame(which(pval < threshold, arr.ind = TRUE),
                       stringAsFactors = FALSE)
  rownames(ind) <- seq_len(nrow(ind)) # TODO test
  cor_pval <- apply(ind, 1, function(x){
    c("cors" = cors[x[1], x[2]],
      "pvalue" = pval[x[1], x[2]])
  })
  ind$row <- rownames(cors)[ind$row]
  ind$col <- colnames(cors)[ind$col]
  ind <- cbind.data.frame(ind, t(cor_pval))
  colnames(ind) <- c("Microorganism", "Gene", "Correlation", "pvalue")

  ind <- ind[!duplicated(ind), ]
  ind <- ind[order(ind$Microorganism, ind$pvalue, decreasing = c(TRUE, FALSE)), ]
  rownames(ind) <- seq_len(nrow(ind))
  ind
}


#' List the correlations
#'
#' Prepare the correlations to be shared (without filtering them)
#' @inheritParams filter_values
#' @return A data.frame with microorganisms, genes, their correlations and the
#' pvalue
#' @export
#' @seealso [filter_values()], [relevant()]
sign_cor <- function(cors, pval, threshold = 0.05) {
  if (sum(pval < threshold, na.rm = TRUE) == 0) {
    stop("No relevant correlations with this threshold")
  }

  ind <- as.data.frame(which(pval < threshold, arr.ind = TRUE),
                       stringAsFactors = FALSE)
  rownames(ind) <- seq_len(nrow(ind)) # TODO test
  cor_pval <- apply(ind, 1, function(x){
    c("cors" = cors[x[1], x[2]],
      "pvalue" = pval[x[1], x[2]])
  })
  ind$row <- rownames(cors)[ind$row]
  ind$col <- colnames(cors)[ind$col]
  ind <- cbind(ind, t(cor_pval))
  colnames(ind) <- c("Microorganism", "Gene", "Correlation", "pvalue")

  ind <- ind[!duplicated(ind), ]
  ind <- ind[order(ind$Microorganism, ind$pvalue, decreasing = c(TRUE, FALSE)), ]
  rownames(ind) <- seq_len(nrow(ind))
  ind
}

#' Plot a correlation with ggplot
#'
#' Use filter_values to decide how to and create an html file
#' @param file The name of the file with the model
#' @inheritParams filter_values
#' @param label Label for the output html file
#' @note Expects genes in rows and species at the columns in `cors` and
#' `pvalue`
#' @return Create a html file at Figures/heatmap...html
#' @export
#' @importFrom heatmaply heatmaply
plot_cor <- function(file, cors, pval, threshold, label) {
  l <- filter_values(file, cors, pval, threshold)
  cors <- l$cors

  cors <- cors[!duplicated(rownames(cors)), ]
  colors_g <- ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                                            midpoint = 0, limits = c(-1, 1))
  heatmaply::heatmaply(cors, name = "Cor",
            ylab = "Genes",
            xlab = "Genus",
            scale_fill_gradient_fun = colors_g,
            file = paste0("Figures/heatmap", label,".html"))

}


#' Translate ensmbl to symbols
#'
#' @param x A data.frame with a column named Gene
#' @return The same data.frame with gene names without duplications and empty
#' symbols
#' @export
ensembl2symbol <- function(x) {
  all_samples_symbol <- x
  all_samples_symbol$Gene <- trimVer(all_samples_symbol$Gene)
  all_samples_symbol$Gene <- mapIds(org.Hs.eg.db, keys = all_samples_symbol$Gene,
                                    keytype = "ENSEMBL", column = "SYMBOL")
  all_samples_symbol <- all_samples_symbol[!is.na(all_samples_symbol$Gene), ]
  all_samples_symbol[!duplicated(all_samples_symbol), ]
}

#' Translate symbols and store
#'
#' @inheritParams ensembl2symbol
#' @param file The name of the file were to store them
#' @seealso [ensembl2symbol()]
#' @return NULL
#' @export
write_cor <- function(x, file){
  all_samples_symbol <- ensembl2symbol(x)
  write.csv(all_samples_symbol, file = file, row.names = FALSE, na = "")
}


#' Store result by microorganism
#'
#' @param x Name of microorganism
#' @param o2 Filter samples
#' @param label Label for the files
#' @return NULL
#' @export
store_micro <- function(x, o2, label) {
  if (nrow(o2[[x]]) > 1 && !is.null(o2[[x]])) {
    write.csv(o2[[x]], row.names = FALSE, na = "",
              file = paste0(x, label, ".csv"))
  }
}

#' Filter genes
#'
#' @param file Name of the file with folder
#' @param expr Matrix of expression
#' @return A Matrix with the filtered genes
#' @export
select_genes_int <- function(file, expr) {
  # Load data from all the patients

  comp1 <- readSGCCA(file)

  keepGenes <- rownames(expr) %in% names(comp1)
  expr[keepGenes, ]
}
