library("ggplot2")
library("patchwork")
library("RColorBrewer")
source("utils.R")

## Colorblind-friendly palette
cols <- brewer.pal(8, "Set2")

#### ---- Processing of the datasets ----####

datasets <- c(
    "specht2019v3", "dou2019_mouse", "zhu2019EL",
    "liang2020_hela", "schoof2021", "leduc2022",
    "derks2022", "brunner2022"
)

processedData <- lapply(datasets, prepareData)
names(processedData) <- datasets
featureTableList <- lapply(processedData, featureTables, byAssay = FALSE)
idMatrixList <- lapply(featureTableList, featureTablesToIdMatrices)

#### ---- Missingness in features ----####

pepMisDf <- lapply(names(idMatrixList), function(n) {
    x <- idMatrixList[[n]]$peptide
    mis <- rowMeans(x == 0)
    mis <- mis[order(mis)]
    data.frame(
        missingness = mis,
        index = seq_along(mis),
        propFeatures = seq_along(mis) / length(mis),
        dataset = n
    )
})
pepMisDf <- do.call(rbind, pepMisDf)
examplePoint <- pepMisDf[58700, , drop = FALSE]
(pl <- ggplot(pepMisDf) +
    aes(
        x = propFeatures * 100,
        y = missingness * 100,
        colour = dataset
    ) +
    geom_line(linewidth = 1) +
    geom_point(
        data = examplePoint, colour = "red",
        shape = 21, size = 4
    ) +
    geom_hline(
        yintercept = examplePoint$missingness * 100,
        linetype = "dashed", colour = "grey70"
    ) +
    geom_vline(
        xintercept = examplePoint$propFeatures * 100,
        linetype = "dashed", colour = "grey70"
    ) +
    scale_color_manual(values = cols) +
    ylab("% of missing values") +
    xlab("% of peptides") +
    theme_minimal())

if (!dir.exists("figs")) {
    dir.create("figs")
}
ggsave("figs/figure1.pdf", plot = pl, width = 5, height = 4)

