library("scpdata")
library("scp")
library("Matrix")
library("tidyr")
library("dplyr")
library("tidyr")
library("dplyr")
library("nipals")
library("scater")

#### ---- Functions to prepare datasets from scpdata

filterSingleCellAssays <- function(object) {
    ## Keep only assays that contain PSM at single cell data
    ## - Remove all "peptides" assays
    ## - Remove all "proteins" assays
    ## - (only derks2022) remove: precusor data (with "plexDIA") and "bulk"
    ## - (only cong2020AC) remove: "HeLa cell" (=bulk) and "Blank"
    idx <- grep("pep|protein|plexDIA|bulk|HeLa[ ]cell|Blank",
        names(object),
        ignore.case = TRUE
    )
    removeAssay(object, idx)
}

identifyProcessingSoftware <- function(object) {
    rdn <- rowDataNames(object)[[1]]
    if ("Sequence" %in% rdn) {
        "maxquant"
    } else if ("Annotated.Sequence" %in% rdn) {
        "pd"
    } else if ("Stripped.Sequence" %in% rdn) {
        "diann"
    } else if ("sequence" %in% rdn) {
        "msgf"
    } else {
        stop("Software not identified")
    }
}

contaminantFilter <- function(ppSoftware) {
    if (ppSoftware == "maxquant") {
        "Potential.contaminant != '+'"
    } else if (ppSoftware %in% c("pd", "msgf")) {
        "!isContaminant"
    } else {
        NULL
    }
}

decoyFilter <- function(ppSoftware) {
    if (ppSoftware == "maxquant") {
        "Reverse != '+'"
    } else if (ppSoftware == "msgf") {
        "!isDecoy"
    } else {
        NULL
    }
}

qvalueFilter <- function(ppSoftware) {
    if (ppSoftware == "maxquant") {
        "PEP < 0.01"
    } else if (ppSoftware == "pd") {
        "Percolator.q.Value < 0.01"
    } else if (ppSoftware == "diann") {
        "Lib.Q.Value < 0.01"
    } else if (ppSoftware == "msgf") "MS.GF.QValue < 0.01"
}

filterFormula <- function(ppSoftware) {
    str <- paste(
        c(
            contaminantFilter(ppSoftware),
            decoyFilter(ppSoftware),
            qvalueFilter(ppSoftware)
        ),
        collapse = "&"
    )
    str <- gsub("[&]{2,}", "&", str)
    as.formula(paste("~", str))
}

cellFilter <- function(cd, dataset) {
    switch(dataset,
        specht2019v3 = grepl("Macro|Mono", cd$SampleType),
        dou2019_lysates = cd$SampleType == "Lysate",
        dou2019_boosting = cd$SampleType %in% c("RAW", "C10", "SVEC"),
        dou2019_mouse = cd$SampleType %in% c("RAW", "C10", "SVEC"),
        zhu2019EL = cd$Cells.per.well == 1,
        liang2020_hela = cd$nbCells == 1,
        schoof2021 = cd$nbCells == 1,
        leduc2022 = grepl("Mono|Mela", cd$SampleType),
        derks2022 = cd$Real_single_cell & cd$Celltype != "Neg",
        brunner2022 = TRUE,
        stop("dataset not recognized")
    )
}

filterSingleCells <- function(object, dataset) {
    sel <- cellFilter(colData(object), dataset)
    if (!all(sel)) object <- subsetByColData(object, sel)
    object[, , ncols(object) > 0]
}

minimalProcessing <- function(object) {
    idx <- names(object)
    object <- zeroIsNA(object, idx)
    object <- filterNA(object, idx, pNA = 0.9999)
    object[, , nrows(object) > 0]
}


## Important note: the derks2022 data will not be treated like other
## dataset. Instead of processing the fragment data, we process the
## "extracted" data to take advantage of the within run ID propagation
## algorithm (see Derks et Slavov 2023).
prepareData <- function(dataset) {
    object <- get(dataset)()
    ppSoftware <- identifyProcessingSoftware(object)
    object <- filterSingleCellAssays(object)
    if (dataset != "derks2022") {
        object <- filterFeatures(object, filterFormula(ppSoftware),
            keep = TRUE
        )
    }
    object <- filterSingleCells(object, dataset)
    object <- formatSampleType(object, dataset)
    minimalProcessing(object)
}


#### ---- Functions to generate an identification table from assays

featureTablesToIdMatrices <- function(featureTables) {
    out <- lapply(
        featureTables,
        featureTableToSparseIdMatrix
    )
    names(out) <- names(featureTables)
    out
}

featureTableToSparseIdMatrix <- function(featureTable) {
    coln <- unique(featureTable[[sampleColName()]])
    rown <- unique(featureTable[[featureColName()]])
    featInSample <- split(
        featureTable[[featureColName()]],
        featureTable[[sampleColName()]]
    )
    out <- Matrix(0,
        nrow = length(rown),
        ncol = length(coln),
        dimnames = list(rown, coln),
        sparse = TRUE
    )
    for (j in seq_along(coln)) {
        out[, j] <- as.numeric(rown %in% featInSample[[j]])
    }
    out
}

jaccardIndexFromIdMatrix <- function(idMatrix) {
    vectorSizes <- colSums(idMatrix)
    pwVectorSizes <- sapply(vectorSizes, function(x) vectorSizes + x)
    unionSize <- crossprod(idMatrix)
    jacc <- unionSize / (pwVectorSizes - unionSize)
    jacc[lower.tri(jacc, diag = TRUE)] <- 0
    as(jacc, "CsparseMatrix")
}

peptideSequenceCol <- function(ppSoftware) {
    switch(ppSoftware,
        maxquant = "Sequence",
        pd = "Annotated.Sequence",
        diann = "Stripped.Sequence",
        msgf = "sequence",
        stop("software not recognized")
    )
}

proteinNameCol <- function(ppSoftware) {
    switch(ppSoftware,
        maxquant = "Leading.proteins",
        pd = "Protein.Accessions",
        diann = "Protein.Names",
        msgf = "DatabaseAccess",
        stop("dataset not recognized")
    )
}

sampleColName <- function() "name"
featureColName <- function() "feature"
featureLevels <- function() c("peptide", "protein")

featureTableByCell <- function(object, pepcol, protcol) {
    out <- lapply(names(object), function(nameAssay) {
        x <- object[[nameAssay]]
        rd <- rowData(x)[, c(pepcol, protcol)]
        featInCellList <- lapply(colnames(x), function(nameCol) {
            isPresent <- !is.na(assay(x)[, nameCol])
            rdPresent <- rd[isPresent, ]
            rdPresent$name <- nameCol
            rdPresent
        })
        do.call(rbind, featInCellList)
    })
    do.call(rbind, out)
}

featureTableToFeatureList <- function(featureDf, pepcol, protcol) {
    featureDf <- featureDf[, c(sampleColName(), pepcol, protcol)]
    out <- lapply(c(pepcol, protcol), function(fcol) {
        subDf <- featureDf[, c(sampleColName(), fcol)]
        colnames(subDf)[2] <- featureColName()
        subDf
    })
    names(out) <- featureLevels()
    out
}

featureTables <- function(object, byAssay) {
    ppSoftware <- identifyProcessingSoftware(object)
    pepcol <- peptideSequenceCol(ppSoftware)
    protcol <- proteinNameCol(ppSoftware)
    if (!byAssay) {
        featureDf <- featureTableByCell(object, pepcol, protcol)
    } else {
        featureDf <- rbindRowData(object, names(object))
        colnames(featureDf)[colnames(featureDf) == "assay"] <- sampleColName()
    }
    featureTableToFeatureList(featureDf, pepcol, protcol)
}

getSteps <- function(assays, nsteps) {
    maxn <- length(assays)
    steps <- seq(1, maxn, length.out = min(maxn, nsteps))
    round(steps)
}

cumulativeSensitivity <- function(idMatrix, niters, nsteps) {
    out <- list()
    nSeq <- getSteps(colnames(idMatrix), nsteps)
    pb <- txtProgressBar(min = 0, max = length(nSeq), style = 3)
    for (n in nSeq) {
        for (i in 1:niters) {
            sel <- sample(colnames(idMatrix), n)
            s <- sum(rowSums(idMatrix[, sel, drop = FALSE]) != 0)
            out <- c(
                out,
                list(data.frame(
                    i = i,
                    NumberSamples = n,
                    Sensitivity = s
                ))
            )
        }
        setTxtProgressBar(pb, which(nSeq == n))
    }
    close(pb)
    do.call(rbind, out)
}

formatSampleType <- function(object, dataset) {
    if (dataset == "brunner2022") {
        object$SampleType <- recode(object$CellCycleStage,
            TB = "G1-S",
            NB = "G2-M",
            UB = "untreated"
        )
    } else if (dataset == "derks2022") {
        object$SampleType <- recode(object$Celltype,
            Melanoma_t = "Melanoma",
            PDAC_t = "PDAC",
            `U-937_t` = "U-937"
        )
    } else if (dataset == "schoof2021") {
        object$SampleType <- object$Population
    } else if (dataset == "zhu2019EL") {
        object$SampleType <- recode(object$FM1.43.signal,
            High = "FM143High",
            Low = "FM143Low"
        )
    } else if (dataset %in% c(
        "dou2019_mouse", "leduc2022",
        "liang2020_hela", "specht2019v3"
    )) {
        NULL
    } else {
        stop(dataset, " not recognized")
    }
    object
}

correctBatch <- function(object, i, name,
                         batchCol1 = NULL, batchCol2 = NULL,
                         biolCol) {
    sce <- getWithColData(object, i)
    form <- as.formula(paste("~", biolCol))
    model <- model.matrix(form, data = colData(sce))
    batch1 <- if (!is.null(batchCol1)) {
        colData(sce)[, batchCol1]
    } else {
        NULL
    }
    batch2 <- if (!is.null(batchCol2)) {
        colData(sce)[, batchCol2]
    } else {
        NULL
    }
    assay(sce) <- removeBatchEffect(
        x = assay(sce),
        batch = batch1,
        batch2 = batch2,
        design = model
    )
    colData(sce) <- NULL
    object <- addAssay(object, y = sce, name = name)
    addAssayLinkOneToOne(object, from = i, to = name)
}


addAnnotationToCorrelation <- function(x, groupBy) {
    if (!all(unique(x$i) %in% rownames(groupBy))) {
        warning(
            "x and groupBy do not match, returning ",
            "correlation without annotations"
        )
        return(x)
    }
    groupBy <- as.data.frame(groupBy)
    gnames <- colnames(groupBy)
    for (ii in gnames) {
        x[[ii]] <- paste("between", ii)
        for (jj in unique(groupBy[, ii])) {
            sel <- groupBy[x$i, ii] == jj & groupBy[x$j, ii] == jj
            x[[ii]][sel] <- paste("within", jj)
        }
    }
    return(x)
}

correlationTable <- function(x, MARGIN = 2, groupBy = NULL) {
    if (MARGIN == 1) x <- t(x)
    corx <- as(
        cor(x, use = "pairwise.complete.obs"),
        "TsparseMatrix"
    )
    corx <- data.frame(
        cor = corx@x,
        i = colnames(x)[corx@i + 1],
        j = colnames(x)[corx@j + 1]
    )
    corx <- corx[corx$i != corx$j, ]
    if (!is.null(groupBy)) {
        corx <- addAnnotationToCorrelation(corx, groupBy)
    }
    return(corx)
}

binProteinsByPercentMissing <- function(object, ngroups = 5) {
    nNArows <- nNA(object)$nNArows
    probs <- seq(0, 1, length.out = ngroups + 1)
    breaks <- quantile(nNArows$pNA, probs = probs)
    out <- cut(nNArows$pNA, breaks = breaks)
    names(out) <- nNArows$name
    out
}
