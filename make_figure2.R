library("ggplot2")
library("patchwork")
source("utils.R")
# devtools::install_github("Mengbo-Li/protDP")
library("protDP")

#### ---- Minimal processing of the datasets ----####

datasets <- c(
    "specht2019v3", "dou2019_mouse", "zhu2019EL",
    "liang2020_hela", "schoof2021", "leduc2022",
    "derks2022", "brunner2022"
)

processedData <- lapply(datasets, prepareData)
names(processedData) <- datasets
featureTableList <- lapply(processedData, featureTables, byAssay = FALSE)
idMatrixList <- lapply(featureTableList, featureTablesToIdMatrices)

# load("data/processedData.RData")

#### ---- Processing Brunner et al. 2022 ----####

brunner <- processedData$brunner2022
experiments(brunner) <- endoapply(experiments(brunner), function(x) {
    rownames(x) <- rowData(x)$Precursor.Id
    x
})
brunner <- joinAssays(brunner, i = names(brunner), name = "precursor")
brunner <- logTransform(brunner, i = "precursor", name = "precursor_log")
brunner <- aggregateFeatures(brunner,
                             i = "precursor_log",
                             name = "peptides",
                             fcol = "Stripped.Sequence",
                             fun = colMedians,
                             na.rm = TRUE)
peptides_brunner2022 <- getWithColData(brunner, "precursor_log")
# save(peptides_brunner2022, file = "data/peptides_brunner2022.RData")

#### ---- Derks et al. 2022 ----####

derks <- processedData$derks2022[, , 2]
derks <- logTransform(derks,
                      i = names(derks),
                      name = "log"
)
derks <- aggregateFeatures(derks,
                           i = "log",
                           name = "peptides",
                           fcol = "Stripped.Sequence",
                           fun = colMedians,
                           na.rm = TRUE
)
peptides_derks2022 <- getWithColData(derks, "peptides")
# save(peptides_derks2022, file = "data/peptides_derks2022.RData")

#### ---- Leduc et al. 2022 ----####

leduc <- processedData$leduc2022
anames <- names(leduc)
leduc <- logTransform(leduc,
                      
                      i = anames,
                      name = paste0(anames, "_log")
)
leduc <- aggregateFeatures(leduc,
                           i = paste0(anames, "_log"),
                           name = paste0(anames, "_peptides"),
                           fcol = "Sequence",
                           fun = colMedians,
                           na.rm = TRUE
)
leduc <- joinAssays(leduc,
                    i = paste0(anames, "_peptides"),
                    name = "peptides")
peptides_leduc2022 <- getWithColData(leduc, "peptides")
# save(peptides_leduc2022, file = "data/peptides_leduc2022.RData")

#### ---- Schoof et al. 2021 ----####

schoof <- processedData$schoof2021
anames <- names(schoof)
schoof <- logTransform(schoof,
                      i = anames,
                      name = paste0(anames, "_log")
)
schoof <- aggregateFeatures(schoof,
                            i = paste0(anames, "_log"),
                            name = paste0(anames, "_peptides"),
                            fcol = "Annotated.Sequence",
                            fun = colMedians,
                            na.rm = TRUE
)
schoof <- joinAssays(schoof,
                     i = paste0(anames, "_peptides"),
                     name = "peptides")
peptides_schoof2021 <- getWithColData(schoof, "peptides")
# save(peptides_schoof2021, file = "data/peptides_schoof2021.RData")

#### ---- Create plot ----####

# load("data/peptides_brunner2022.RData")
# load("data/peptides_derks2022.RData")
# load("data/peptides_leduc2022.RData")
# load("data/peptides_schoof2021.RData")

peps <- c(
    "peptides_brunner2022", "peptides_derks2022", 
    "peptides_leduc2022", "peptides_schoof2021"
)
pl <- lapply(peps, function(i) {
    x <- assay(get(i))
    
    allc <- union(rownames(datasetB$prec), rownames(x))
    table(dpc = allc%in%rownames(datasetB$prec), scpdata = allc %in% rownames(x))
    allc <- union(colnames(datasetB$prec), colnames(x))
    colnames(x) <- sapply(colnames(x), function(x) paste0(head(tail(strsplit(x, "")[[1]], 12), 10), collapse = ""))
    table(dpc = allc %in% colnames(datasetB$prec), scpdata = allc %in% colnames(x))
    x <- x[rownames(datasetB$prec), colnames(datasetB$prec)]
    
    x <- table()
    
    df <- data.frame(
        intensity = rowMeans(x, na.rm = TRUE),
        completeness = 1 - rowMeans(is.na(x))
    )
    dpcfit <- dpc(x, b1.upper = 1000)
    df$dpc <- plogis(dpcfit$beta[1] + dpcfit$beta[2] * df$intensity)
    ggplot(df) +
        aes(x = intensity, y = completeness) +
        geom_point(size = 0.2, alpha = 0.05) +
        geom_line(aes(y = dpc), color = "red") +
        # stat_smooth(
        #     method = "glm",
        #     formula = "y ~ 1 + x",
        #     method.args = list(family = quasibinomial(link = "logit")),
        #     se = FALSE) +
        ggtitle(sub("peptides_", "", i)) +
        xlab("Log peptide intensity") +
        ylab("Proportion detected") +
        theme_minimal()
})
(pl <- wrap_plots(pl) +
        plot_annotation(tag_levels = "a"))
ggsave("figs/figure2.pdf", plot = pl, width = 6, height = 6)
