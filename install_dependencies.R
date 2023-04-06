if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c(
    "dplyr",
    "ggplot2",
    "ggrepel",
    "impute",
    "limma",
    "Matrix",
    "nipals",
    "patchwork",
    "RColorBrewer",
    "scater",
    "scp",
    "scpdata",
    "tidyr"
))
