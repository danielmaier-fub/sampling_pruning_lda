############################################################################################
## Title: How document sampling and vocabulary pruning affect the results of topic models ##
## Authors: Maier, D., Niekler, A., Wiedemann, G., Stoltenberg, D. #########################
## Journal: Computational Communication Research ###########################################
## Last Update: 2019-11-19 #################################################################
## Code: Wiedemann, G., Niekler, A., Maier, D. #############################################
############################################################################################

source("aux_functions.r")
ppath <- "models"

corpus <- c("/coab", "/guardian", "/twitter")
mode <- c("/unpruned/", "/pruned/")
sSize <- c("0.01", "0.05", "0.1", "0.2", "0.5")
topWordsToMatch <- 20
similarityThreshold <- 0.5

start <- Sys.time()
for (g in 1:length(corpus)){
  for (h in 1:length(mode)){
    cNr <- g 
    mNr <- h
    path <- paste0(ppath, corpus[cNr], mode[mNr])
    
    # Load sample models
    tms <- list()
    print(paste0("loading sample models ... corpus: ", corpus[cNr], " mode: ", mode[mNr]))
    start <- Sys.time()
    for(j in 1:length(sSize)){
      pathToTM <- paste0(path, "sample/", sSize[j], "/")
      files <- list.files(pathToTM)
      tms[[j]] <- list()
      for (o in 1:length(files)){
        tms[[j]][[o]] <- get(load(file = paste0(pathToTM, files[o])))
      }
    }
    print(start-Sys.time())
    
    # Load reference models
    pathToTMref <- paste0(path, "reference/")
    files <- list.files(pathToTMref)
    tmr <- list()
    print(paste0("loading reference models ... corpus: ", corpus[cNr], " mode: ", mode[mNr]))
    for (o in 1:length(files)){
      tmr[[o]] <- get(load(file = paste0(pathToTMref, files[o])))
    }
    
    # start of analysis
    sSize <- c("0.01", "0.05", "0.1", "0.2", "0.5")
    modelCount <- matrix(NA, nrow = length(sSize))
    
    for(j in 1:length(sSize)){
      sizes <- j
      pathToTM <- paste0(path, "sample/", sSize[sizes])
      nModels <- length(list.files(pathToTM))
      modelCount[j,] <- nModels
    }
    
    rownames(modelCount) <- sSize
    folder <- "reference/"
    modelCountFull <- matrix(NA, nrow = length(folder))
    pathToTM <- paste0(path, folder)
    Full <- length(list.files(pathToTM))
    modelCount <- rbind(modelCount, Full)
    colnames(modelCount) <- paste0(corpus[cNr], mode[mNr])
    print(modelCount)
    
    # Internal Reliability of Sample Models
    sSize <- c("0.01/", "0.05/", "0.1/", "0.2/", "0.5/")
    resList <- list()
    for(z in 1:length(sSize)){
      pathToTM <- paste0(path, "sample/", sSize[z])
      nModels <- length(list.files(pathToTM))
      vec <- 1:nModels 
      comb <- expand.grid(vec,vec)
      ncomb <- comb[!duplicated(t(apply(comb, 1, sort))),]
      ncomb <- as.matrix(ncomb[-which(ncomb[,1] == ncomb[,2]),])
      y <- c()
      for(j in 1:length(ncomb[,1])){
        alignment <- alignTopicModels2(tms[[z]][[ncomb[j,1]]], tms[[z]][[ncomb[j,2]]], 
                                       topWordsToMatch = topWordsToMatch, 
                                       similarityThreshold = similarityThreshold)
        y <- c(y, alignment$reliability)
        printAlignedTopics(alignment)
      }
      resList[[z]] <- y
    }
    
    m <- c()
    for (z in 1:length(resList)){
      zeile <- resList[[z]]
      m <- rbind(m, zeile)
    }
    rownames(m) <- c("0.01", "0.05", "0.1", "0.2", "0.5")
    
    # convert to long format
    Size <- rep(c(0.01, 0.05, 0.1, 0.2, 0.5), each = length(m[1,]))
    Reli <- as.numeric(t(m))
    Corpus <- rep(corpus[cNr], length(Reli))
    Prune_Mode <- rep(mode[mNr], length(Reli))
    CorpMode <- paste0(Corpus, Prune_Mode)
    df_samp <- data.frame(CorpMode, Size, Reli)
    df_samp_sum <- summarySE(df_samp, measurevar = "Reli", groupvar = c("Size", "CorpMode"))
    
    # Internal Reliability of Full Models	
    folder <- c("reference/")
    resListRef <- list()
    pathToTM <- paste0(path, folder)
    nModels <- length(list.files(pathToTM))
    
    vec <- 1:nModels 
    comb <- expand.grid(vec,vec)
    ncomb <- comb[!duplicated(t(apply(comb, 1, sort))),]
    ncomb <- as.matrix(ncomb[-which(ncomb[,1] == ncomb[,2]),])
    y <- c()
    
    for(j in 1:length(ncomb[,1])){
      alignment <- alignTopicModels(tmr[[ncomb[j,1]]], tmr[[ncomb[j,2]]], 
                                    topWordsToMatch = topWordsToMatch, 
                                    similarityThreshold = similarityThreshold)
      y <- c(y, alignment$reliability)
      printAlignedTopics(alignment)
    }
    resListRef <- y
    
    m <- resListRef
    Size <- rep(1, length(m))
    Reli <- as.numeric(t(m))
    Corpus <- rep(corpus[cNr], length(Reli))
    Prune_Mode <- rep(mode[mNr], length(Reli))
    CorpMode <- paste0(Corpus, Prune_Mode)
    df_ref <- data.frame(CorpMode, Size, Reli)
    df_ref_sum <- summarySE(df_ref, measurevar = "Reli", groupvar = c("Size", "CorpMode"))
    
    df_int <- rbind(df_samp_sum, df_ref_sum)
    
    save(df_int, file = paste0("analysis/", c("chall_", "guard_", "twitt_")[cNr],
                               c("_unpruned", "_pruned_")[mNr], "internal_reliability", ".RData"))
    
    # Comparison of Sample Models and Full Corpus Models	
    folder <- "sample/"
    sSize <- c("0.01/", "0.05/", "0.1/", "0.2/", "0.5/")
    ref <- c("reference/")
    
    resListRep <- list() 
    for(z in 1:length(sSize)){
      pathToTM <- paste0(path, folder, sSize[z])
      pathToTMref <- paste0(path, ref)
      nModels <- length(list.files(pathToTM))
      vec <- 1:nModels 
      ncomb <- expand.grid(vec,vec)
      y <- c()
      for(j in 1:length(ncomb[,1])){ 
        alignment <- alignTopicModels2(tms[[z]][[ncomb[j,1]]], tmr[[ncomb[j,2]]], 
                                       topWordsToMatch = topWordsToMatch, 
                                       similarityThreshold = similarityThreshold)
        y <- c(y, alignment$reliability)
        printAlignedTopics(alignment)
      }
      resListRep[[z]] <- y
    }
    
    comp_rep <- resListRep
    
    df_rep <- data.frame()
    for(z in 1:length(sSize)){
      m <- comp_rep[[z]]
      Size <- rep(sSize[z], length(m))
      Reli <- as.numeric(t(m))
      Corpus <- rep(corpus[cNr], length(Reli))
      Prune_Mode <- rep(mode[mNr], length(Reli))
      CorpMode <- paste0(Corpus, Prune_Mode)
      tmp_df_rep <- data.frame(CorpMode, Size, Reli)
      df_rep <- rbind(df_rep, tmp_df_rep)
    }
    
    df_rep_sum <- summarySE(df_rep, measurevar = "Reli", groupvar = c("Size", "CorpMode"))
    df_rep_sum$Size <- rep(c(0.01, 0.05, 0.1, 0.2, 0.5))
    
    save(df_rep_sum, file = paste0("analysis/", c("chall_", "guard_", "twitt_")[cNr],
                                   c("unpruned_", "pruned_")[mNr], "reproducibility",".RData"))
  }
}
Sys.time()-start

# Construct final data list "fin_df"
# fin_df is a list, where list item 1
# contains reproducibility values
# list item 2 contains reliability

fin_df <- list()

#which(str_detect(list.files("analysis"), "reproducibility"))
#list.files("analysis")
rep_path <- "analysis"
#files <- list.files(rep_path)
files <- list.files("analysis")[which(str_detect(list.files("analysis"), "reproducibility"))]
#files <- files[-which(files == "plots")]
for (i in 1:length(files)){
  rep_df <- get(load(paste0(paste0("analysis/",files[i]))))
  if (i == 1){
    fin_df[[1]] <- rep_df
  }else{
    fin_df[[1]] <- rbind(fin_df[[1]], rep_df)
  }
}


rel_path <- "analysis"
#files <- list.files(rep_path)
files <- list.files("analysis")[which(str_detect(list.files("analysis"), "reliability"))]
#rel_path <- "results/reliability/"
#files <- list.files(rel_path)
#files <- files[-which(files == "plots")]
for (i in 1:length(files)){
  rel_df <- get(load(paste0("analysis/",files[i])))
  if (i == 1){
    fin_df[[2]] <- rel_df
  }else{
    fin_df[[2]] <- rbind(fin_df[[2]], rel_df)
  }
}

fin_df[[1]]$CorpMode <- stringr::str_replace(fin_df[[1]]$CorpMode, "/unpruned/", " (unpruned)")
fin_df[[1]]$CorpMode <- stringr::str_replace(fin_df[[1]]$CorpMode, "/pruned/", " (pruned)")
fin_df[[2]]$CorpMode <- stringr::str_replace(fin_df[[2]]$CorpMode, "/unpruned/", " (unpruned)")
fin_df[[2]]$CorpMode <- stringr::str_replace(fin_df[[2]]$CorpMode, "/pruned/", " (pruned)")
fin_df[[1]]$CorpMode <- stringr::str_replace(fin_df[[1]]$CorpMode, "/twitter", "Twitter")
fin_df[[1]]$CorpMode <- stringr::str_replace(fin_df[[1]]$CorpMode, "/guardian", "News")
fin_df[[1]]$CorpMode <- stringr::str_replace(fin_df[[1]]$CorpMode, "/challenger", "Websites")
fin_df[[2]]$CorpMode <- stringr::str_replace(fin_df[[2]]$CorpMode, "/twitter", "Twitter")
fin_df[[2]]$CorpMode <- stringr::str_replace(fin_df[[2]]$CorpMode, "/guardian", "News")
fin_df[[2]]$CorpMode <- stringr::str_replace(fin_df[[2]]$CorpMode, "/challenger", "Websites")

names(fin_df[[1]])[2] <- "Corpus"
names(fin_df[[2]])[2] <- "Corpus"

# Figures

sequence <- list()
sequence[[1]] <- list()
sequence[[1]][[1]] <- c(1:10)
sequence[[1]][[2]] <- c(11:20)
sequence[[1]][[3]] <- c(21:30)
sequence[[2]] <- list()
sequence[[2]][[1]] <- c(1:12)
sequence[[2]][[2]] <- c(13:24)
sequence[[2]][[3]] <- c(25:36)

res_path <- "analysis/"

for (i in 1:3){
  pdf(file = paste0(res_path, c("challenger", "news", "twitter")[i], "_reliability.pdf"),
      width = 9, height = 6)
  
  ggplot(fin_df[[2]][sequence[[2]][[i]],], aes(x=Size, y=Reli, color = Corpus)) + 
    geom_errorbar(aes(ymin=Reli-se, ymax=Reli+se), width=.025) +
    geom_line() +
    geom_point() +
    coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(0.57, 0.92), expand = F) +
    scale_y_discrete(limits = c(0.6, 0.7, 0.8, 0.9), 
                     labels = c("0.6", "0.7", "0.8", "0.9")) +
    scale_x_discrete(limits = c(0.01, 0.05, 0.1, 0.2, 0.5, 1), 
                     labels = c("1", "5", "10", "20", "50", "100")) +
    
    labs(title = "Reliability of Topics",
         subtitle = paste0(topWordsToMatch," Topwords, t = ", similarityThreshold),
         x = "Sample Size in %", y = "Share of matched Topics") 
  
  dev.off()
  
}

for (i in 1:3){
  pdf(file = paste0(res_path, c("challenger", "news", "twitter")[i], "_reproducibility.pdf"),
      width = 9, height = 6)
  
  ggplot(fin_df[[1]][sequence[[1]][[i]],], aes(x=Size, y=Reli, color = Corpus)) + 
    geom_errorbar(aes(ymin=Reli-se, ymax=Reli+se), width=.025) +
    geom_line() +
    geom_point() +
    coord_cartesian(xlim = c(-0.02, 0.52), ylim = c(0.37, 0.92), expand = F) +
    scale_y_discrete(limits = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 
                     labels = c("0.4", "0.5", "0.6", "0.7", "0.8", "0.9")) +
    scale_x_discrete(limits = c(0.01, 0.05, 0.1, 0.2, 0.5), 
                     labels = c("1", "5", "10", "20", "50")) +
    
    labs(title = "Reproducibility of Models by Sample Size",
         subtitle = paste0(topWordsToMatch," Topwords, t = ", similarityThreshold),
         x = "Sample Size", y = "Share of matched Topics")
  
  dev.off()
}

# Pruning Reproducibility

resListPrun <- list()
for(s in 1:length(corpus)){
  resListPrun[[s]] <- list()
  pathToTMprun <- paste0(ppath, corpus[s], mode[1])
  pathToTMUnprun <- paste0(ppath, corpus[s], mode[2])
  for(z in 1:length(sSize)){
    pathToTM1 <- paste0(pathToTMprun, "sample/", sSize[z])
    pathToTM2 <- paste0(pathToTMUnprun, "sample/", sSize[z])
    nModels <- 2
    vec <- 1:nModels 
    ncomb <- expand.grid(vec,vec)
    y <- c()
    for(j in 1:length(ncomb[,1])){ 
      tm1 <- get(load(file = paste0(pathToTM1, list.files(pathToTM1)[ncomb[j,1]])))
      tm2 <- get(load(file = paste0(pathToTM2, list.files(pathToTM2)[ncomb[j,2]])))
      alignment <- alignTopicModels2(tm1, tm2, topWordsToMatch = topWordsToMatch, similarityThreshold = similarityThreshold)
      y <- c(y, alignment$reliability)
      printAlignedTopics(alignment)
    }
    resListPrun[[s]][[z]] <- y
  }
}

save(resListPrun, file = "analysis/resListPrun.RData")

sSize <- c("0.01", "0.05", "0.1", "0.2", "0.5")
df_prun <- data.frame()
for (s in 1:length(corpus)){
  for(z in 1:length(sSize)){
    m <- resListPrun[[s]][[z]]
    Size <- rep(sSize[z], length(m))
    Reli <- as.numeric(t(m))
    if(s == 1){
      Corpus <- rep("Websites", length(Reli))
    }
    if(s == 2){
      Corpus <- rep("News", length(Reli))
    }
    if(s == 3){
      Corpus <- rep("Twitter", length(Reli))
    }
    tmp_df_prun <- data.frame(Corpus, Size, Reli)
    df_prun <- rbind(df_prun, tmp_df_prun)
  }
}

df_prun <- summarySE(df_prun, measurevar = "Reli", groupvar = c("Size", "Corpus"))
df_prun$Size <- rep(c(0.01, 0.05, 0.1, 0.2, 0.5), each = 3)


res_path <- "analysis/"

pdf(file = paste0(res_path, "pruning.pdf"),
    width = 9, height = 6)

ggplot(df_prun, aes(x=Size, y=Reli, color = Corpus)) + 
  geom_errorbar(aes(ymin=Reli-se, ymax=Reli+se), width=.025) +
  geom_line() +
  geom_point() +
  
  labs(title = "Reproducibility Unpruned from Pruned Models",
       subtitle = paste0(topWordsToMatch," Topwords, t = ", similarityThreshold),
       x = "Sample Size", y = "Share of matched topics")

dev.off()
