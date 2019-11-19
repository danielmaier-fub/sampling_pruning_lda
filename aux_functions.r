############################################################################################
## Title: How document sampling and vocabulary pruning affect the results of topic models ##
## Authors: Maier, D., Niekler, A., Wiedemann, G., Stoltenberg, D. #########################
## Journal: Computational Communication Research ###########################################
## Last Update: 2019-11-19 #################################################################
## Code: Wiedemann, G., Niekler, A., Maier, D. #############################################
############################################################################################

# install.packages("lda", "topicmodels", "htmltools", "Matrix", "xtable", "cluster", "igraph","slam","data,table","quanteda")

library("methods")
require(lda)
require(topicmodels)
require(htmltools)
require(xtable)
require(igraph)
require(cluster)
require(Matrix)
require(slam)
require(data.table)
require(tictoc)
require(magrittr)
require(ggplot2)
require(plyr)
options(stringsAsFactors = F)


#####################################################
## functions to run model evaluation ################
#####################################################

runLDAModelEvalution <- function(dtm, run_number, K, iterations, alphaPriorsToTest, etaPriorsToTest, blacklist = NULL, fraction = NULL, fixedSample=FALSE, storageDir = storageDir) {
  # Read input data
  fullDTM <- dtm
  if (!is.null(fraction)) {
    sampleSize <- fraction * nrow(dtm)
    if (fixedSample == TRUE){
      set.seed(123)
    }
    docSubset <- sample(1:nrow(dtm), sampleSize)
    if (fixedSample == TRUE){
      set.seed(sample(1:10000000, 1))
    }
    dtm <- dtm[docSubset,]
  }
  
  corpusData <- c()
  corpusVocab <- c()
  
  if (!is.null(blacklist)) {
    corpusDTM <- dtm[, !((colnames(dtm) %in% blacklist))]
    corpusDTM <- corpusDTM[,colSums(corpusDTM)>0]
    corpusDTM <- corpusDTM[rowSums(corpusDTM)>0,]
    corpusLDA <- dtm2ldaformat(makeSimpleTripletMatrix(corpusDTM))
    corpusData <- corpusLDA$documents
    corpusVocab <- corpusLDA$vocab
  }
  else {
    corpusDTM <- dtm
    corpusDTM <- corpusDTM[,colSums(corpusDTM)>0]
    corpusDTM <- corpusDTM[rowSums(corpusDTM)>0,]
    corpusLDA <- dtm2ldaformat(makeSimpleTripletMatrix(corpusDTM))
    corpusData <- corpusLDA$documents
    corpusVocab <- corpusLDA$vocab
  }
  
  # Prepare output data
  modelEvaluationData <- data.frame(modelID = integer(), K = integer(), alpha = double(), eta = double(), modelLikelihood = double(), modelCoherence = integer(), topicsHTML = character(), modelFileName = character())
  modelEvaluationHTML <- data.frame()
  # Run evaluation
  modelID <- run_number
  for (alpha in alphaPriorsToTest) {
    for (eta in etaPriorsToTest) {
      #modelID <- modelID + 1
      t <- as.numeric(Sys.time())
      seed <- 1e8 * (t - floor(t))
      set.seed(seed); print(seed)
      topicModel <- lda.collapsed.gibbs.sampler(corpusData, K, vocab = corpusVocab, num.iterations = iterations, alpha = alpha, eta = eta, initial = NULL, trace = 1L, compute.log.likelihood = T)
      topicTerms <- top.topic.words(topicModel$topics, 25, by.score=TRUE)
      topicProportions <- t(topicModel$document_sums) / colSums(topicModel$document_sums)
      modelLikelihood <- tail(as.vector(topicModel$log.likelihoods[2, ]), 1)
      topicCoherenceForAllTopics <- topicCoherence(ldaformat2Matrix(corpusData, corpusVocab), topicModel)
      modelCoherence <- mean(topicCoherenceForAllTopics)
      
      # Prob
      tProportions <- colSums(topicProportions) / nrow(topicProportions)
      oProportions <- order(tProportions, decreasing = T)
      # Rank 1
      firstDocTopics <- apply(topicProportions, 1, FUN=function(x) order(x, decreasing=TRUE)[1])
      primaryDocTopics <- factor(firstDocTopics, 1:K)
      nRanks1 <- table(primaryDocTopics)
      tRanks1 <- as.integer(nRanks1)
      oRanks1 <- order(tRanks1, decreasing = T)
      # Coherence
      tCoherence <- topicCoherenceForAllTopics
      oCoherence <- order(topicCoherenceForAllTopics, decreasing = T)
      topics <- data.frame(
        TopicID = sprintf("%02d", 1:K),
        #r_r1 = oRanks1,
        Rank1 = tRanks1,
        #r_pr = oProportions,
        Prob = tProportions,
        #r_c = oCoherence,
        Coherence = tCoherence,
        Terms = paste0("<pre>", apply(topicTerms, 2, paste0, collapse=" "), "</pre>")
      )
      topics <- topics[oRanks1, ]
      topicsHTML <- print(xtable(topics, digits=5), print.results = F, type="html", sanitize.text.function=function(x){x}, include.rownames = F, html.table.attributes = paste0('id="T', modelID, '" class="sortable"'))
      
      if (is.null(fraction) == FALSE) {
        ftype <- "sample" 
        fractionName <- fraction
      }else{
        ftype <- "reference"
        fractionName <- ""
      }
      
      storageDirectory <- paste0(storageDir, "/", ftype, "/", fractionName, "/") 
      fname <- paste0(corp , "_", pmode, "_", ftype, "_", fractionName, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),".RData" )
      
      modelFileName <- paste0(storageDirectory, fname)
      save(topicModel, file = modelFileName)
      
      modelEvaluationHTML <- rbind(modelEvaluationHTML, data.frame(modelID, K, alpha, eta, modelLikelihood, modelCoherence, topicsHTML, modelFileName))
      print(paste0("LDA Model ", modelID, " has been saved to ", modelFileName))
      
      storageDirectoryEval <- paste0(storageDir, "/eval/") 
      evalFileName <- paste0(storageDirectoryEval,"eval_", corp , "_", pmode, "_", ftype, "_", fractionName, ".html" )
      
      sink(evalFileName)
      print(xtable(modelEvaluationHTML, digits = 4), file=evalFileName, append=T, sanitize.text.function=function(x){x}, include.rownames = F)
      sink()
    }
  }
  
}


topicCoherence <- function(DTM, ldaModel, N = 25) {
  
  # Ensure matrix or Matrix-format (convert if SparseM)
  if (is.simple_triplet_matrix(DTM)) {
    DTM <- sparseMatrix(i=DTM$i, j=DTM$j, x=DTM$v, dims=c(DTM$nrow, DTM$ncol), dimnames = dimnames(DTM))
  }
  
  DTMBIN <- DTM
  
  if (any(DTMBIN > 1)) {
    DTMBIN <- DTMBIN >= 1 + 0
  }
  
  #DTMBIN[DTMBIN > 0] <- 1
  
  documentFrequency <- colSums(DTMBIN)
  names(documentFrequency) <- colnames(DTMBIN)
  
  K <- nrow(ldaModel$topics)
  
  topNtermsPerTopic <- top.topic.words(ldaModel$topics, N, by.score=TRUE)
  allTopicModelTerms <- unique(as.vector(topNtermsPerTopic))
  
  DTMpreprocessed <- DTMBIN[, allTopicModelTerms]
  DTMpreprocessedCooc <- t(DTMpreprocessed) %*% DTMpreprocessed
  DTMpreprocessedCooc <- t((DTMpreprocessedCooc + 1) / colSums(DTMpreprocessed))
  DTMpreprocessedCooc <- log(DTMpreprocessedCooc)
  DTMpreprocessedCooc <- as.matrix(DTMpreprocessedCooc)
  
  coherence <- rep(0, K)
  pb <- txtProgressBar(max = K)
  for (topicIdx in 1:K) {
    setTxtProgressBar(pb, topicIdx)
    topWordsOfTopic <- topNtermsPerTopic[, topicIdx]
    coherence[topicIdx] <- 0
    for (m in 2:length(topWordsOfTopic)) {
      for (l in 1:(m-1)) {
        mTerm <- topWordsOfTopic[m]
        lTerm <- topWordsOfTopic[l]
        coherence[topicIdx] <- coherence[topicIdx] + DTMpreprocessedCooc[mTerm, lTerm]
      }
    }
  }
  close(pb)
  return(coherence)
}


ldaformat2Matrix <- function (documents, vocab) {
  #require(slam)
  stm <- simple_triplet_matrix(
    i = rep(seq_along(documents), sapply(documents, ncol)), 
    j = as.integer(unlist(lapply(documents,"[", 1, )) + 1L), 
    v = as.integer(unlist(lapply(documents,"[", 2, ))), 
    nrow = length(documents), 
    ncol = length(vocab), 
    dimnames = list(names(documents), vocab))
  dtm <- sparseMatrix(i=stm$i, j=stm$j, x=stm$v, dims=c(stm$nrow, stm$ncol), dimnames = dimnames(stm))
}


peekIntoModelCorpus <- function(sampleCorpusFulltext, topicModel, topicToInvestigate, topicThresholdInDocument = NULL, n = 1) {
  topicProportions <- t(topicModel$document_sums) / colSums(topicModel$document_sums)
  if (!is.null(topicThresholdInDocument)) {
    idx <- which(topicProportions[, topicToInvestigate] > topicThresholdInDocument)
    docSampledIds <- sample(idx, min(n, length(idx)), n)
  } else {
    docSampledIds <- order(topicProportions[, topicToInvestigate], decreasing = T)[1:n]
  }
  sampleTexts <- sampleCorpusFulltext[docSampledIds]
  html_print(HTML(paste0(sampleTexts, collapse = "<hr><br><br/><br/>")))
}


getModelAlignments <- function(modelEvaluationData, topWordsToMatch = 100, similarityThreshold = 0.2, verbose = F) {
  numModels <- nrow(modelEvaluationData)
  if (numModels < 2) stop("Nothing to compare, got just than one model!")
  cat(c("Parameters: nWords", topWordsToMatch, "| threshold", similarityThreshold, "\n"))
  pairs <- combn(as.character(modelEvaluationData$modelFileName), 2, simplify = F)
  allReliabilities <- rep(0, length(pairs))
  i <- 0
  for (pair in pairs) {
    i <- i + 1
    cat(c("------", "\n"))
    cat(c(pair[1], "\n"))
    cat(c(pair[2], "\n"))
    tm1 <- get(load(file = pair[1]))
    tm2 <- get(load(file = pair[2]))
    alignment <- alignTopicModels(tm1, tm2, topWordsToMatch, similarityThreshold)
    printAlignedTopics(alignment, verbose = verbose)
    allReliabilities[i] <- alignment$reliability
  }
  return(allReliabilities)
}


#####################################################
## functions to compare topic models ################
#####################################################

topicProbability <- function(topics) {
  token <- sum(topics)
  prob <- c()
  for(i in 1:nrow(topics)) {
    ts <- sum(topics[i, ]) / token
    prob <- c(prob, ts)
  }
  return(prob)
}


toComparable <- function(topics, topWordsToMatch) {
  for(i in 1:nrow(topics)) {
    ts <- topics[i, ] / sum(topics[i, ])
    left <- names(ts[order(-ts)][1:topWordsToMatch])
    topics[i, ] <- 0
    topics[i, left] <- ts[left]
  }
  return(topics)
}


TM_Aligner <- function(topics1, topics2, thres, probs1, probs2, topWordsToMatch) {
  K <- nrow(topics1)
  reliability <- 0
  cosineDists <- as.matrix(1 - topics1 %*% t(topics2) / (sqrt(rowSums(topics1 ^ 2) %*% t(rowSums(topics2 ^ 2)))))
  minIndexes <- apply(cosineDists, 2, which.min)
  mins <- apply(cosineDists, 2, min)
  
  alignedTopics <-list()
  alignedTopics$ids <- matrix(0, nrow = K, ncol = 2)
  alignedTopics$probabilities <- matrix(0, nrow = K, ncol = 2)
  alignedTopics$sharedTerms <- vector("list", K)
  alignedTopics$distance <- rep(0, K)
  
  for (i in 1:K) {
    
    index <- arrayInd(which.min(cosineDists), .dim = c(K, K))
    value <- min(cosineDists)
    if (value > thres)
      break
    
    cosineDists[index[1], ] <- 1
    cosineDists[, index[2]] <- 1
    
    alignedTopics$ids[i, ] <- c(index[1], index[2])
    alignedTopics$probabilities[i, ] <- c( probs1[index[1]],  probs1[index[2]])
    alignedTopics$sharedTerms[[i]] <- intersect(
      names(topics1[index[1], ][order(-topics1[index[1], ])][1:topWordsToMatch]), 
      names(topics2[index[2], ][order(-topics2[index[2], ])][1:topWordsToMatch]))
    
    reliability = reliability + 1
    
    alignedTopics$distance[i] <- value
    
  }
  alignedTopics$reliability <- reliability / K
  return(alignedTopics)
}


alignTopicModels <- function(tm1, tm2, topWordsToMatch = 50, similarityThreshold = 0.2) {
  c_topics1 <- toComparable(tm1$topics, topWordsToMatch)
  c_topics_p1 <- topicProbability(tm1$topics)
  c_topics2 <- toComparable(tm2$topics, topWordsToMatch)
  c_topics_p2 <- topicProbability(tm2$topics)
  alignedTopics <- TM_Aligner(c_topics1, c_topics2, similarityThreshold, c_topics_p1, c_topics_p2, topWordsToMatch)
  return(alignedTopics)
}


printAlignedTopics <- function(alignedTopics, verbose = F) {
  if (verbose) {
    for (i in 1:length(alignedTopics$sharedTerms)) {
      if (length(alignedTopics$sharedTerms[[i]]) > 0) {
        cat(c("___________________________________________________________","\n"))
        cat(c("Shared terms:", alignedTopics$sharedTerms[[i]], "\n"))
        cat(c("Distance:", sprintf("%.4f", alignedTopics$distance[i]),"\n"))
        cat(c("Alignment:", alignedTopics$ids[i, 1], "TO", alignedTopics$ids[i, 2], "\n"))
        cat(c("Probabilities:", sprintf("%.4f", alignedTopics$probabilities[i, 1]), "TO:", sprintf("%.4f", alignedTopics$probabilities[i, 2]), "\n"))
      }
    }
    cat(c("===========================================================", "\n"))
  }
  cat(c("RELIABILITY:", sprintf("%.4f", alignedTopics$reliability), "\n"))
}

RUN_AS_MAIN <- FALSE
if (RUN_AS_MAIN) {
  # Load Models
  tm1 <- get(load(file = "../online_climate_uk/model-1_K40_a0.002_e0.025_i10000rnd_d2015-12-18_16-28-51.RData"))
  tm2 <- get(load(file = "../online_climate_uk/model-1_K40_a0.002_e0.025_i10000fix_d2015-12-16_01-55-11.RData"))
  alignment <- alignTopicModels(tm1, tm2, topWordsToMatch = 50, similarityThreshold = 0.3)
  printAlignedTopics(alignment)
}

#####################################################
## MISCELLANEOUS ####################################
#####################################################

create_required_folders <- function(wd = getwd(), download.data = FALSE, create.folders = FALSE){
  if (download.data == TRUE){
    dir.create("data")
    download.file(url = "https://maier.userpage.fu-berlin.de/CCR_2019/dtm_coab.Rdata", destfile = paste0(getwd(),"/data/dtm_coab.Rdata"))
    download.file(url = "https://maier.userpage.fu-berlin.de/CCR_2019/dtm_guardian.Rdata", destfile = paste0(getwd(),"/data/dtm_guardian.Rdata"))
    download.file(url = "https://maier.userpage.fu-berlin.de/CCR_2019/dtm_twitter.Rdata", destfile = paste0(getwd(),"/data/dtm_twitter.Rdata"))
  }
  if(create.folders == TRUE){
    dir.create("analysis")
    dir.create("models")
    corpora <- c("coab", "guardian", "twitter")
    pmode <- c("pruned", "unpruned")
    ftype <- c("eval", "reference", "sample")
    ssize <- c("0.01", "0.05", "0.1", "0.2", "0.5")
    
    for (i in 1:length(corpora)){
      path <- paste0("models/", corpora[i])
      dir.create(path)
      for (j in 1:length(pmode)){
        path <- paste0("models/", corpora[i],"/", pmode[j])
        dir.create(path)
        for (k in 1:length(ftype)){
          path <- paste0("models/", corpora[i],"/", pmode[j], "/", ftype[k])
          dir.create(path)
          if (ftype[k] == "sample"){
            for (n in 1:length(ssize)){
              path <- paste0("models/", corpora[i],"/", pmode[j], "/", ftype[k], "/", ssize[n])
              dir.create(path)
            }
          }
        }
      }
    }
  }
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



TM_Aligner2 <- function(topics1, topics2, thres, probs1, probs2, topWordsToMatch) {
  K <- nrow(topics1)
  reliability <- 0
  topics2 <- topics2[,intersect(colnames(topics1), colnames(topics2))] # ??????????????????
  topics1 <- topics1[,intersect(colnames(topics1), colnames(topics2))] # ??????????????????
  cosineDists <- as.matrix(1 - topics1 %*% t(topics2) / (sqrt(rowSums(topics1 ^ 2) %*% t(rowSums(topics2 ^ 2))))) # Problem des Modellvergleichs bei unterschiedlicher Vokabularlänge
  minIndexes <- apply(cosineDists, 2, which.min)
  mins <- apply(cosineDists, 2, min)
  
  alignedTopics <-list()
  alignedTopics$ids <- matrix(0, nrow = K, ncol = 2)
  alignedTopics$probabilities <- matrix(0, nrow = K, ncol = 2)
  alignedTopics$sharedTerms <- vector("list", K)
  alignedTopics$distance <- rep(0, K)
  
  for (i in 1:K) {
    
    index <- arrayInd(which.min(cosineDists), .dim = c(K, K))
    value <- min(cosineDists)
    if (value > thres)
      break
    
    cosineDists[index[1], ] <- 1
    cosineDists[, index[2]] <- 1
    
    alignedTopics$ids[i, ] <- c(index[1], index[2])
    alignedTopics$probabilities[i, ] <- c( probs1[index[1]],  probs1[index[2]])
    alignedTopics$sharedTerms[[i]] <- intersect(
      names(topics1[index[1], ][order(-topics1[index[1], ])][1:topWordsToMatch]), 
      names(topics2[index[2], ][order(-topics2[index[2], ])][1:topWordsToMatch]))
    
    reliability = reliability + 1
    
    alignedTopics$distance[i] <- value
    
  }
  alignedTopics$reliability <- reliability / K
  return(alignedTopics)
}

alignTopicModels2 <- function(tm1, tm2, topWordsToMatch, similarityThreshold) {
  c_topics1 <- toComparable(tm1$topics, topWordsToMatch)
  c_topics_p1 <- topicProbability(tm1$topics)
  c_topics2 <- toComparable(tm2$topics, topWordsToMatch)
  c_topics_p2 <- topicProbability(tm2$topics)
  alignedTopics <- TM_Aligner2(c_topics1, c_topics2, similarityThreshold, c_topics_p1, c_topics_p2, topWordsToMatch)
  return(alignedTopics)
}

makeSimpleTripletMatrix <- function(sparseDTMMatrix){
  tmp <- Matrix::summary(sparseDTMMatrix)
  return(slam::simple_triplet_matrix(i=tmp[,1], j=tmp[,2], v=tmp[,3],
                                     dimnames=dimnames(sparseDTMMatrix),nrow =
                                       nrow(sparseDTMMatrix), ncol =
                                       ncol(sparseDTMMatrix)))
}
