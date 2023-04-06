library("ggplot2")
library("patchwork")
library("ggrepel")
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

#### ---- Cumulative sensitivity curves ----####

featLevel <- "peptide"
i <- "schoof2021"
cd <- colData(processedData[[i]])
cumSensPerSampleTypeSchoof <- lapply(unique(cd$SampleType), function(st) {
    selSamples <- rownames(cd)[cd$SampleType == st]
    id <- idMatrixList[[i]][[featLevel]]
    out <- cumulativeSensitivity(id[, colnames(id) %in% selSamples],
        niters = 10, nsteps = 50
    )
    out$SampleType <- st
    out
})
dfCumSensitivity <- do.call(rbind, cumSensPerSampleTypeSchoof)

## Plot sensitivity
dfCumSensitivity <- filter(
    dfCumSensitivity,
    SampleType %in% c("BLAST CD38-", "LSC")
)
## Asymptotic regression model
dfLSC <- dfCumSensitivity[dfCumSensitivity$SampleType == "LSC", ]
AsymptoticRegressionModel <- nls(
    formula = Sensitivity ~ SSasymp(NumberSamples, Asym, R0, lrc), 
    data = dfLSC, 
    ## weight more observation with more samples
    weights = dfLSC$NumberSamples^4 
)
predictedSensitivity <- data.frame(
    NumberSamples = 1:max(dfCumSensitivity$NumberSamples) 
)
predictedSensitivity$Sensitivity <- predict(
    AsymptoticRegressionModel,
    newdata = predictedSensitivity
)
localSens <- median(filter(dfCumSensitivity, NumberSamples == 1)$Sensitivity)
totalSens <- median(filter(dfCumSensitivity, NumberSamples == max(NumberSamples))$Sensitivity)
(plCumS <- ggplot(dfCumSensitivity) +
    aes(
        x = NumberSamples,
        y = Sensitivity,
        colour = SampleType
    ) +
    geom_point(size = 1) +
    geom_line(
        data = predictedSensitivity,
        colour = "red",
        linewidth = 0.5
    ) +
    geom_hline(
        yintercept = c(localSens, totalSens),
        linetype = c("dashed", "solid"),
        color = "grey60"
    ) +
    scale_colour_manual(
        name = element_blank(),
        values = c(cols[[6]], "gold4")
    ) +
    ylim(0, 1.025 * totalSens) +
    labs(
        y = "Sensitivity",
        x = "Number of cells"
    ) +
    theme_minimal() +
    theme(
        legend.position = "top",
        axis.title.x = element_text(vjust = 0)
    ))

#### ---- Metrics ----####

featLevel <- "peptide"
df <- list()
for (dataset in names(processedData)) {
    cd <- colData(processedData[[dataset]])
    for (st in unique(cd$SampleType)) {
        selSamples <- rownames(cd)[cd$SampleType == st]
        id <- idMatrixList[[dataset]][[featLevel]][, selSamples]
        df <- rbind(
            df,
            data.frame(
                LocalSensitivityMedian = median(colSums(id)),
                LocalSensitivitySd = sd(colSums(id)),
                TotalSensitivity = sum(rowSums(id) != 0),
                Completeness = mean(id != 0),
                NumberCells = ncol(id),
                Dataset = dataset,
                SampleType = st
            )
        )
    }
}

localSensitivityDf <- list()
for (dataset in names(processedData)) {
    cd <- colData(processedData[[dataset]])
    id <- idMatrixList[[dataset]][[featLevel]]
    localSensitivityDf <- rbind(
        localSensitivityDf,
        data.frame(
            LocalSensitivity = colSums(id),
            SampleType = cd[colnames(id), "SampleType"],
            Dataset = dataset
        )
    )
}

JaccardDf <- list()
for (dataset in names(processedData)) {
    cd <- colData(processedData[[dataset]])
    for (st in unique(cd$SampleType)) {
        selSamples <- rownames(cd)[cd$SampleType == st]
        id <- idMatrixList[[dataset]][[featLevel]][, selSamples]
        jaccMat <- jaccardIndexFromIdMatrix(id)
        JaccardDf <- rbind(
            JaccardDf,
            data.frame(
                JaccardIndex = jaccMat@x,
                SampleType = st,
                Dataset = dataset
            )
        )
    }
}

(plTotalSens <- ggplot(df) +
    aes(
        y = TotalSensitivity,
        x = SampleType,
        fill = Dataset
    ) +
    geom_col() +
    facet_grid(~Dataset, scales = "free_x") +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Total sensitivity") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ))

(plLocalSens <- ggplot(localSensitivityDf) +
    aes(
        y = LocalSensitivity,
        x = SampleType,
        fill = Dataset
    ) +
    geom_violin() +
    facet_grid(~Dataset, scales = "free_x") +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Local sensitivity") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ))

(plJaccard <- ggplot(JaccardDf) +
    aes(
        y = JaccardIndex,
        x = SampleType,
        fill = Dataset
    ) +
    geom_violin() +
    facet_grid(~Dataset, scales = "free_x") +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Jaccard index") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ))

(plMetrics <- ggplot(df) +
    aes(
        x = Completeness * 100,
        y = LocalSensitivityMedian,
        colour = Dataset
    ) +
    geom_errorbar(aes(
        ymin = LocalSensitivityMedian - LocalSensitivitySd,
        ymax = LocalSensitivityMedian + LocalSensitivitySd
    )) +
    geom_point(aes(size = NumberCells)) +
    geom_text_repel(aes(label = SampleType), size = 3.25) +
    geom_point(
        data = . %>%
            group_by(Dataset) %>%
            mutate(
                Completeness = mean(Completeness),
                LocalSensitivityMedian = mean(LocalSensitivityMedian)
            ),
        size = 5,
        shape = 1
    ) +
    labs(
        x = "% data completeness",
        y = "Local sensitivity",
        size = "Number of cells"
    ) +
    scale_colour_manual(values = cols, guide = NULL) +
    scale_size_continuous(breaks = c(10, 50, 100, 500, 1000)) +
    theme_minimal())

#### ---- Final plot ----####

(fig <- plCumS +
    plMetrics +
    plTotalSens +
    plJaccard + theme(legend.position = "none") +
    plot_layout(design = "122
                           122
                           122
                           122
                           333
                           444") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 20)))

if (!dir.exists("figs")) {
    dir.create("figs")
}
ggsave("figs/figure2.pdf", plot = fig, width = 10, height = 9)
