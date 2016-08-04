
#!/usr/bin/env Rscript
library(optparse)
library(parallel)


optionList =  list(
  make_option(c('-g', '--genes'), type='character', default=NULL,
              help='List of tenes to check', metavar='character'),
  make_option(c('-c', '--cancer_types'),        type='character', default=NULL,
              help='Cancer types to check', metavar = 'character'),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-o", "--out"), type="character", default="./",
              help="output folder name [default= %default]", metavar="character")
);

optParser = OptionParser(option_list=optionList);
opt = parse_args(optParser);

# opt$cancer_types = '../data/cancer_types_to_check.sample.txt'
# opt$genes = '../data/genes_to_check.sample.txt'

if (is.null(opt$cancer_types)){
  print_help(optParser)
  stop("Cancer types catalog should be supplied", call.=FALSE)
}
if (is.null(opt$genes)){
  print_help(optParser)
  stop("Genes catalog should be supplied", call.=FALSE)
}


INPUT_DATA_FOLDER   = 'data'
EXPRESSION_FILENAME = 'Cell_line_RMA_proc_basalExp.txt'
EXPRESSIONDATA_FILENAME = 'expressionData.Rdata'
AUC                 = 'mm5_auc.csv'
IC50                = 'mm5_ic50.csv'
COMPOUNDS_FILENAME  = 'mmc2-compounds.csv'
TREATMENT_FILES     = c(AUC, IC50)


SAMPLES_BY_SUBSET = 'mm5_samples_by_cancer_type.csv'

CANCER.TYPE = 'Cancer.Type'
TISSUE.DESCRIPTION = 'Tissue.description'

GENES_TO_CHECK_SUBSETS = opt$genes
CANCER_TYPES_TO_CHECK_SUBSETS = opt$cancer_types

talk <- function(message){
  if(opt$verbose){
    print(message)
  }
}


getSamplesBySubset <- function(samplesData, subsetRules, type = CANCER.TYPE){

  subsets = unique(subsetRules$SubsetName)
  validRuleColumns =  c(CANCER.TYPE, TISSUE.DESCRIPTION)
  if(!type %in% validRuleColumns){
    stop('Invalid rule column')
  }
  samples = sapply(subsets, function(group, subsetRules, samplesData){
    groupToRetrieve = as.vector(
      subsetRules[subsetRules$SubsetName == group,][type][[1]]
    )

    if(groupToRetrieve == 'ALL') {
      identifiers = samplesData$Identifier
    } else {
      identifiers = samplesData[samplesData[,type]  %in% groupToRetrieve, ]$Identifier
    }
    return(as.vector(identifiers))


  }, subsetRules = subsetRules,
  samplesData = samplesData)
  names(samples) = subsets
  return(samples)
}

getGenesBySubset  <- function(genesToCheckSubsets){

  subsets = unique(genesToCheckSubsets$SubsetName)
  genes = lapply(subsets, function(subset, genesToCheckSubsets){
    return(
      as.vector(genesToCheckSubsets[genesToCheckSubsets$SubsetName == subset,]$Gene)
    )
  }, genesToCheckSubsets = genesToCheckSubsets)
  names(genes) = subsets
  return(genes)
}


loadTreatment  <- function(treatmentFileName){

  treatment = read.csv(treatmentFileName, sep=';')

  drugs = data.frame(
    Drug_id   = treatment$Drug_id,
    Drug_name = treatment$Drug_name
  )

  treatment = treatment[-1:-2]
  treatment = as.data.frame(t(treatment))

  colnames(treatment) = paste(drugs$Drug_id, drugs$Drug_name, sep='#')
  rownames(treatment) = as.vector(
    sapply(rownames(treatment), function(name){
      return (strsplit(name, 'X')[[1]][2])}
    )
  )
  return(treatment)
}

normalizeCellLines <- function(cellLines){
  as.vector(
    sapply(cellLines, function(cellLine){
      return(strsplit(cellLine, '\\.')[[1]][2])
    })
  )
}

normalizeGeneNames  <- function(geneNames){
  # returns string w/o trailing whitespace
  trim.trailing <- function (x) sub("\\s+$", "", x)
  as.vector(
    sapply(geneNames, function(name){
      return(
        gsub(';', ',',
             trim.trailing(
               strsplit(as.character(name), '\\[')[[1]][1])
        )
      )
    })
  )
}


getGeneSymbols <- function(expressionData){
  ifelse(as.character(expressionData$GENE_SYMBOLS) != '',
         as.character(expressionData$GENE_SYMBOLS),
         paste('NA',substr(normalizeGeneNames(expressionData$GENE_title), 0, 15 ), sep='-'))
}

cleanExpressionData <- function(expressionData){
  expressionData[(expressionData$GENE_title != '' |
                    expressionData$GENE_SYMBOLS != ''),]
}



asExpressionMatrix <- function(expressionData){

  geneSymbols = getGeneSymbols(expressionData)
  expressionDataMat = expressionData[-1:-2]
  colnames(expressionDataMat) = normalizeCellLines(colnames(expressionDataMat))
  expressionDataMat = as.data.frame(t(expressionDataMat))
  colnames(expressionDataMat) = geneSymbols

  return(expressionDataMat)
}

#
# loopNormalParallel <- function(treatmentsSort, expressionDataMatSort){
#   noCores = 7
#   cl = makeCluster(noCores)
#   correlationsPMatrix = parSapply(cl, 1:ncol(treatmentsSort), function(treatmentI, treatmentsSort, expressionDataMatSort){
#     treatment = as.numeric(treatmentsSort[[treatmentI]])
#     sapply(1:ncol(expressionDataMatSort), function(expressionI, expressionDataMatSort, treatment){
#       signif(cor.test(treatment,
#                       as.numeric(expressionDataMatSort[,expressionI]))$p.value)
#     }, expressionDataMatSort = expressionDataMatSort,
#     treatment = treatment)
#   }, treatmentsSort = treatmentsSort,
#   expressionDataMatSort = expressionDataMatSort
#   )
#   stopCluster(cl)
#
#   colnames(correlationsPMatrix) = colnames(treatmentsSort)
#   rownames(correlationsPMatrix) = colnames(expressionDataMatSort)
#   return(correlationsPMatrix)
# }


makeCorTest<- function(geneValues, treatmentValues, geneName, treatmentName){
  getTreatmentId <- function(treatmentCode){
    return(strsplit(treatmentCode, '#')[[1]][1])
  }

  tryCatch({
    corTest = cor.test(geneValues, treatmentValues)
  }, error   = function(e){
    err <<- e
  }, warning = function(e){
    err <<- e
  } )

  return(
    data.frame(
      gene             = geneName,
      # # gene.title = gene[['GENE_title']],
      treatmentId      = getTreatmentId(treatmentName),
      treatmentName    = treatmentName,
      statistic        = ifelse(exists('corTest') == TRUE, signif(corTest$estimate),NA),
      p.value          = ifelse(exists('corTest') == TRUE, signif(corTest$p.value),NA),
      n.gene           = sum(!is.na(geneValues)),
      n.treatment      = sum(!is.na(treatmentValues)),
      error            = ifelse(exists('corTest') == FALSE, err[['message']], NA),
      stringsAsFactors = FALSE
    )
  )
}

loopNormalSequentialDetailed <- function(treatmentsSort, expressionDataMatSort, talk){


  correlationsPMatrix = sapply((1:ncol(treatmentsSort)), function(treatmentI,
                                                                  treatmentsSort,
                                                                  expressionDataMatSort,
                                                                  makeCorTest){
    talk(paste('Check for [', treatmentI, ']', sep=''))
    treatment = treatmentsSort[treatmentI]

    treatmentData = sapply(1:ncol(expressionDataMatSort), function(expressionI,
                                                                   expressionDataMatSort,
                                                                   treatment){
      makeCorTest(as.numeric(treatment[[1]]),
                  as.numeric(expressionDataMatSort[,expressionI]),
                  colnames(expressionDataMatSort[expressionI]),
                  colnames(treatment[1]))
    }, expressionDataMatSort = expressionDataMatSort,
    treatment = treatment)


  }, treatmentsSort = treatmentsSort,
  expressionDataMatSort = expressionDataMatSort,
  makeCorTest = makeCorTest
  )

  # TODO: Read from output
  columns = c('GeneTitle', 'TreatmentId', 'Treatment', 'Correlation',
              'P.value', 'n.gene', 'n.treatment', 'error')
  correlationsPMatrix = t(as.data.frame(matrix(unlist(correlationsPMatrix), nrow=length(columns))))
  colnames(correlationsPMatrix) = columns
  return(correlationsPMatrix)
}

sortCorrelations <- function(correlations){
  correlations[,'Correlation'] = as.numeric(as.character(correlations[,'Correlation']))
  correlations = correlations[is.na(correlations[,'Correlation']) ==FALSE,]
  correlations[order(correlations[,'Correlation']),]
}

getTreatmentName <- function(treatmentCode){
  return(strsplit(treatmentCode, '#')[[1]][2])
}

checkCorrelations <- function(corrToCheck,
                              expressionDataMat,
                              treatment,
                              group,
                              outputDir){

  geneToCheck = as.character(corrToCheck['GeneTitle'])
  treatmentToCheck = as.character(corrToCheck['Treatment'])
  # treatmentToCheck = getTreatmentName(treatmentToCheck)

  geneData = expressionDataMat[geneToCheck]

  treatmentData = treatment[treatmentToCheck]

  merged = merge(geneData, treatmentData, by=0)

  (reg = lm(merged[[treatmentToCheck]] ~ merged[[geneToCheck]] ))
  name = paste(group, '_gene_',geneToCheck,'_treatment_', getTreatmentName(treatmentToCheck), '.png', sep='')
  name = gsub('/', '', name)
  plotsDir = paste(outputDir, '/plots/', sep='')
  if(dir.exists(plotsDir) == FALSE){
    dir.create(plotsDir)
  }
  # cor.test(merged[[geneToCheck]], merged[[treatmentToCheck]])
  png(paste(plotsDir, name, sep='/'))
  plot(merged[[geneToCheck]], merged[[treatmentToCheck]],
       main=paste( "Correlation of", geneToCheck, "~", treatmentToCheck, ":",  corrToCheck['Correlation']),
       xlab= paste('gene',geneToCheck),
       ylab = paste('treatment', treatmentToCheck)
  )
  abline(reg)
  dev.off()
}

addCompounds <- function(correlations, compounds){
  return(merge(correlations, compounds, by.x = 'TreatmentId', by.y='Identifier'))
}

talk('Load gene subset rules')
genesToCheckSubsets = read.csv(GENES_TO_CHECK_SUBSETS,
                               header = TRUE,
                               sep=';')
genesBySubset = getGenesBySubset(genesToCheckSubsets)


talk('Load cancer subset rules')
cancerTypesSubsets =  read.csv(CANCER_TYPES_TO_CHECK_SUBSETS,
                               header = TRUE,
                               sep=';')


samplesData = read.csv(paste(INPUT_DATA_FOLDER, SAMPLES_BY_SUBSET, sep='/' ),
                       sep=';')
samplesBySubset = getSamplesBySubset(samplesData, cancerTypesSubsets, CANCER.TYPE)

compounds = read.csv(paste(INPUT_DATA_FOLDER, COMPOUNDS_FILENAME, sep='/'), sep=';')



talk('Load expression data')
# expressionData = read.csv(paste(INPUT_DATA_FOLDER, EXPRESSION_FILENAME, sep='/'),
#                          sep='\t')
# expressionData = cleanExpressionData(expressionData)
# expressionDataMat = asExpressionMatrix(expressionData)
# expressionDataSignif = signif(expressionData)
# save(expressionDataMatSig, file='data/expressionData.Rdata', compress = 'xz', compression_level = 9)


load(paste(INPUT_DATA_FOLDER, EXPRESSIONDATA_FILENAME, sep='/'))


expressionDataMatSort = expressionDataMatSig[order(rownames(expressionDataMatSig)),]

#Filtering selected genes in order to reduce memory consumption
allGenesToAnalyze = colnames(expressionDataMatSort) %in% unique(genesToCheckSubsets$Gene)

expressionDataMatSort = expressionDataMatSort[ allGenesToAnalyze ]

groupsAndTreatments = expand.grid(group = unique(names(genesBySubset)), treamentId = TREATMENT_FILES)
groupsAndTreatments = groupsAndTreatments[is.na(groupsAndTreatments$group) == FALSE,]
time1 = proc.time()
talk('Start loop')
apply(groupsAndTreatments, 1, function(groupAndTreatment,
                                       samplesBySubset,
                                       genesBySubset,
                                       expressionDataMatSort,
                                       compounds){
  talk(groupAndTreatment)
  groupToRetrieve = groupAndTreatment[['group']]
  treatmentId = groupAndTreatment[['treamentId']]

  treatments = loadTreatment(paste(INPUT_DATA_FOLDER, treatmentId, sep='/'))
  treatmentsSort = treatments[rownames(treatments) %in% samplesBySubset[[groupToRetrieve]],]
  treatmentsSort = treatmentsSort[order(rownames(treatmentsSort)),]

  treatmentsSort = treatmentsSort[rownames(treatmentsSort) %in% rownames(expressionDataMatSort),,
                                  drop=FALSE]

  talk(paste('Lets filter genes', genesBySubset[[groupToRetrieve]], sep =':', collapse=','))

  expressionDataMatSort = expressionDataMatSort[
    ,  colnames(expressionDataMatSort) %in% genesBySubset[[groupToRetrieve]], drop= FALSE
  ]

  expressionDataMatSort = expressionDataMatSort[
    c(rownames(expressionDataMatSort) %in% rownames(treatmentsSort)),, drop=FALSE
  ]

  correlationMatrix = loopNormalSequentialDetailed(treatmentsSort, expressionDataMatSort, talk)

  correlationMatrix <- addCompounds(correlationMatrix, compounds)

  outputDir = paste(opt$out,'/',groupToRetrieve, sep = '')
  if(dir.exists(outputDir)==FALSE){
    dir.create(outputDir)
  }

  talk(paste('save', paste(outputDir, '/corr_', treatmentId, sep='')))

  write.table(correlationMatrix,
              paste(outputDir, '/corr_', treatmentId, sep=''),
              quote = TRUE,
              row.names = FALSE,
              sep=','
  )

  correlationsSorted  = sortCorrelations(correlationMatrix)
  highestCorrelations = head(correlationsSorted)
  lowestCorrelations  = tail(correlationsSorted)

  imagesStores = apply(highestCorrelations, 1, checkCorrelations,
                       expressionDataMat = expressionDataMatSort,
                       treatment         = treatments,
                       group             = as.character(treatmentId),
                       outputDir         = outputDir)

  imagesStores = apply(lowestCorrelations, 1, checkCorrelations,
                       expressionDataMat = expressionDataMatSort,
                       treatment         = treatments,
                       group             = as.character(treatmentId),
                       outputDir         = outputDir)

}, samplesBySubset       = samplesBySubset,
genesBySubset         = genesBySubset,
expressionDataMatSort = expressionDataMatSort,
compounds             = compounds)

print("Done")
proc.time() - time1
