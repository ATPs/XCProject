## define a function to generate a tree vector
GetRandomTreeIn01 <- function(n){
  #given a interger n, return a vector of string with '0' and '1'
  ## check if n is a valid input.
  if (!(n > 0 && n %% 1 ==0)) {
    print('input n is invalid')
    return('')
  } 
  if (n == 1) return ('0')
  t = c('0','1')
  if (n == 2) return (t)
  for (i in 3:n) {
    l = length(t)
    leaf.id <- sample(l,1)
    leaf <- t[leaf.id]
    leaf.l <- paste(leaf,'0',sep = '')
    leaf.r <- paste(leaf,'1',sep = '')
    t[leaf.id] <- leaf.l
    t[l+1] <- leaf.r
  }
  return(sample(t))
}



## function read the input Tree file used for tree comparison
ReadInputTree <- function(filename){
  return(read.table(filename,sep = '\t',header = TRUE, as.is = TRUE, colClasses = "character", quote = ''))
}

ReadTreeComparison <- function(filename){
  f.lines <- readLines(filename, file.info(filename)$size)
  f.length <- length(f.lines)
  t.count <- f.length %/% 13
  id <- rep(NA,t.count)
  Score <- rep(NA,t.count)
  RootS <- rep(NA,t.count)
  RootT <- rep(NA,t.count)
  PruneS <- rep(NA,t.count)
  PruneT <- rep(NA,t.count)
  MatchS <- rep(NA,t.count)
  MatchT <- rep(NA,t.count)
  alignments <- data.frame(id, Score, RootS, RootT, PruneS, PruneT, MatchS, MatchT)
  for (i in 1:t.count){
    l <- i * 9 - 9
    id <- i
    for (j in 1:8){
      alignments[i,j] = f.lines[l+j]
    }
  }
  alignments$id <- sub('{','',alignments$id,fixed = TRUE)
  alignments$Score <- sub('Score:', '', alignments$Score, fixed = TRUE)
  alignments$RootS <- sub('RootS:', '', alignments$RootS, fixed = TRUE)
  alignments$RootT <- sub('RootT:', '', alignments$RootT, fixed = TRUE)
  alignments$PruneS <- sub('PruneS:', '', alignments$PruneS, fixed = TRUE)
  alignments$PruneT <- sub('PruneT:', '', alignments$PruneT, fixed = TRUE)
  alignments$MatchS <- strsplit(sub('MatchS:', '', alignments$MatchS, fixed = TRUE), ' ')
  alignments$MatchT <- strsplit(sub('MatchT:', '', alignments$MatchT, fixed = TRUE), ' ')
  row.names(alignments) <- alignments$id
  alignments <- alignments[,-1]
  return(alignments)
}

## function read the tree comparison file, return a list
ReadTreeComparison <- function(filename){
  f.lines <- readLines(filename, file.info(filename)$size)
  f.length <- length(f.lines)
  t.count <- f.length %/% 13
  id <- rep(NA,t.count)
  Score <- rep(NA,t.count)
  RootS <- rep(NA,t.count)
  RootT <- rep(NA,t.count)
  PruneS <- rep(NA,t.count)
  PruneT <- rep(NA,t.count)
  MatchS <- rep(NA,t.count)
  MatchT <- rep(NA,t.count)
  alignments <- data.frame(id, Score, RootS, RootT, PruneS, PruneT, MatchS, MatchT)
  for (i in 1:t.count){
    l <- i * 9 - 9
    id <- i
    for (j in 1:8){
      alignments[i,j] = f.lines[l+j]
    }
  }
  alignments$id <- sub('{','',alignments$id,fixed = TRUE)
  alignments$Score <- sub('Score:', '', alignments$Score, fixed = TRUE)
  alignments$RootS <- sub('RootS:', '', alignments$RootS, fixed = TRUE)
  alignments$RootT <- sub('RootT:', '', alignments$RootT, fixed = TRUE)
  alignments$PruneS <- sub('PruneS:', '', alignments$PruneS, fixed = TRUE)
  alignments$PruneT <- sub('PruneT:', '', alignments$PruneT, fixed = TRUE)
  alignments$MatchS <- strsplit(sub('MatchS:', '', alignments$MatchS, fixed = TRUE), ' ')
  alignments$MatchT <- strsplit(sub('MatchT:', '', alignments$MatchT, fixed = TRUE), ' ')
  row.names(alignments) <- alignments$id
  alignments <- alignments[,-1]
  return(alignments)
}

## function create n radom tree file for a given input tree
GenerateRandomTreeFile <- function(inputTree, referenceTree, fileNumber = 1, outfilePrefix = './randomTree'){
  #inputTree, the input tree to be randomized. it is in the format of a vector of character in 01 coding the leaves
  #referenceTree, a data.frame, with leaves 01 as row.name, and other information to be recognized by the HSA program
  #fileNumber, the number of random Trees to be generated
  #outfilePrefix, the prefix for the output files
  outfileNames <- paste(outfilePrefix,1:fileNumber,sep = '')
  inputTree.vector <- inputTree[[1]]#convert list to vector
  leavesKeep <- referenceTree[,1] %in% inputTree.vector
  for (i in 1:fileNumber){
    inputTree.df <- referenceTree[leavesKeep,]
    randLineage <- GetRandomTreeIn01(sum(leavesKeep))
    inputTree.df[,1] <- randLineage
    write.table(inputTree.df,sep = '\t', quote = FALSE, row.names = FALSE, file = outfileNames[i])
    a <- write.table(inputTree.df,sep = '\t', quote = FALSE, row.names = FALSE)
  }
  return(NA)
}

