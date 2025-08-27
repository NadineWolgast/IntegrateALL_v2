#!/usr/bin/env Rscript
# Install Bioconductor packages for CNV analysis

cat("Installing Bioconductor packages for CNV analysis...\n")

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.r-project.org/")
}

# List of required Bioconductor packages
bioc_packages <- c(
    "GenomicRanges",
    "IRanges", 
    "GenomicFeatures",
    "rtracklayer",
    "Biostrings",
    "BSgenome",
    "GenomeInfoDb",
    "vcfR"
)

# Install each package
for (pkg in bioc_packages) {
    cat(paste("Installing", pkg, "...\n"))
    tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
        cat(paste(pkg, "installed successfully\n"))
    }, error = function(e) {
        cat(paste("Warning: Failed to install", pkg, ":", e$message, "\n"))
    })
}

# Install RNAseqCNV from the local copy
cat("Installing RNAseqCNV from local source...\n")
tryCatch({
    # First try to install from CRAN/Bioconductor
    BiocManager::install("RNAseqCNV", ask = FALSE, update = FALSE)
    cat("RNAseqCNV installed from Bioconductor\n")
}, error = function(e) {
    cat("Could not install RNAseqCNV from Bioconductor, will use local copy\n")
})

cat("Bioconductor package installation completed!\n")