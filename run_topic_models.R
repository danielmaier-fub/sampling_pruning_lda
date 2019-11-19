############################################################################################
## Title: How document sampling and vocabulary pruning affect the results of topic models ##
## Authors: Maier, D., Niekler, A., Wiedemann, G., Stoltenberg, D. #########################
## Journal: Computational Communication Research ###########################################
## Last Update: 2019-11-19 #################################################################
## Code: Wiedemann, G., Niekler, A., Maier, D. #############################################
############################################################################################


#####################################################
## run topic models #################################
#####################################################

source("aux_functions.r")

# Create infrastructure and download data
create_required_folders(wd = getwd(), download.data = FALSE, create.folders = TRUE)

# define basic terminology
corpora <- c("coab", "guardian", "twitter")
prmode <- c("pruned", "unpruned")

# Choose corpus: 1 = Web, 2 = News, 3 = Twitter
for (g in 1:length(corpora)){
  thisCorpus <- g
  
  for (m in 1:length(prmode)){
    # Choose pruning mode: 1 = unpruned, 2 = pruned
    thisPmode <- m
    
    # Load data 
    corp <- corpora[thisCorpus]
    pmode <- prmode[thisPmode]
    dtmFileName <- paste0("dtm_", corp, ".RData")
    dtm <- get(load(file = paste0("data/", dtmFileName)))  
    
    # path to storage Directory
    storageDir <- paste0("models/", corp , "/", pmode)
    
    # Set number of topics
    K <- 50
    
    # Iterations should be at least 1000 for final evaluation
    iterations <- 1000
    alphaPriorsToTest <- c(0.5)
    etaPriorsToTest <- c(1 / K)
    runs <- 5
    fraction_sizes <- c(0.01, 0.05, 0.1, 0.2, 0.5)
    results <- list()
    
    # Blacklist of some tokens to be ignored
    blacklistBadTokens <- readLines("./blacklist_badtokens.txt", encoding = "UTF-8")
    tokensToIgnore <- blacklistBadTokens %>%
      tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE)  %>%  tokens_tolower() 
    
    if (pmode == "pruned"){
      dtm %<>% dfm_trim(min_docfreq = nrow(dtm)*0.005, max_docfreq = nrow(dtm)*99)
      dtm <- dtm[rowSums(dtm)>5,]
    }
    
    # Sample Models 5 runs, Random Initialization, for 1%,5%,10%, 20%, 50% sample size
    for (frac in fraction_sizes){
      tic(paste(frac))
      for (i in 1:runs){
        run_number <- i
        print(paste0("run: ", run_number, ", frac: ", frac, " IniMeth: rnd"))
        modelEvaluationData <- runLDAModelEvalution(dtm, run_number, K, iterations, alphaPriorsToTest, etaPriorsToTest, 
                                                    blacklist = tokensToIgnore,  
                                                    storageDir = storageDir,
                                                    fraction = frac, fixedSample = TRUE)
      }
      results[[length(results)+1]] <- toc()
    }
    
    #Reference Model 5 runs, Random Initialization
    for (i in 1:runs){
      run_number <- i
      print(paste0("run: ", run_number, " Reference Model", " IniMeth: rnd"))
      tic(paste("1"))
      modelEvaluationData <- runLDAModelEvalution(dtm, run_number, K, iterations, alphaPriorsToTest, etaPriorsToTest, 
                                                  blacklist = tokensToIgnore, 
                                                  storageDir = storageDir,
                                                  fraction = NULL, fixedSample = TRUE)
      results[[length(results)+1]] <- toc()
    }
    
  }
  
}
