
# Code from the RHEA pipeline

# Calculate the species richness in a sample
Species.richness <- function(x) {
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  count <- sum(x > 0.5)
  return(count)
}

# Calculate the Shannon diversity index
Shannon.entropy <- function(x) {
  total <- sum(x)
  se <- -sum(x[x > 0] / total * log(x[x > 0] / total))
  return(se)
}

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x) {
  total <- sum(x)
  se <- round(exp(-sum(x[x > 0] / total * log(x[x > 0] / total))), digits = 2)
  return(se)
}

# Calculate the Simpson diversity index
Simpson.concentration <- function(x) {
  total <- sum(x)
  si <- sum((x[x > 0] / total)^2)
  return(si)
}

# Calculate the effective number of species for Simpson
Simpson.effective <- function(x) {
  total <- sum(x)
  si <- round(1 / sum((x[x > 0] / total)^2), digits = 2)
  return(si)
}
