# Functions to parse the trees.
## modified from Meng Yuan's code


GetParentInBinary <- function(childBinary){
  # childBinary is like "01011", return its parent, which should be "0101"
  parentBinary <- substr(childBinary,1, nchar(childBinary)-1)
  if (parentBinary == "") parentBinary <- "Root"
  return(parentBinary)
}

GetParentsInBinary <- function(childrenBinary){
  # childrenBinary is a vector of childBinary, like c("0101", "0000")
  # return a vector of each of their parents. c("010","000")
  return(unlist(lapply(childrenBinary, GetParentInBinary)))
}

Binary2Decimal <- function(nodeBinary){
  # nodeBinary is the node code with 01, such as "001".
  # convert nodeBinary to decimal number. By first add "1" to the beginning.
  # "001" will be "1001", and return 9.
  # 9 is the ID of the node in a complete binary tree
  if (nodeBinary == "Root" || nodeBinary == "root") return("1")
  return(strtoi(paste0("1", nodeBinary), base = 2))
}

Decimal2Binary <- function(decimalID){
  # the revert function of Binary2decimal
  if (decimalID == 1) return("Root")
  binaryID <- as.integer(rev(intToBits(decimalID)))
  return (sub("^0*1","", paste0(binaryID, collapse = "")))
}

CountNodeDirectBranches <- function(binaryIDs, nodeBinary){
  # binaryIDs is a vector with all IDs of nodes, like "Root", "0","1","00"...
  # nodeBinary is a node code in 01, like "0"
  # return the number of direct children of nodeBinary in binaryIDs.
  # for example, for node "0", if only "01" in binaryIDs, return 1
  n <- 0
  nodeBinary.0 <- paste0(nodeBinary,"0")
  nodeBinary.1 <- paste0(nodeBinary,"1")
  if (nodeBinary.0 %in% binaryIDs) n <- n + 1
  if (nodeBinary.1 %in% binaryIDs) n <- n + 1
  return(n)
}

GetMissingChild <- function(binaryIDs, nodeBinary){
  # if nodeBinary has 1 child in binaryIDs, return the other.
  if (CountNodeDirectBranches(binaryIDs, nodeBinary) == 1){
    nodeBinary.0 <- paste0(nodeBinary,"0")
    nodeBinary.1 <- paste0(nodeBinary,"1")
    return(setdiff(c(nodeBinary.0, nodeBinary.0), binaryIDs))
  }
}

GetMissingChildren <- function(binaryIDs){
  # return a vector with missing children found by GetMissingChild
  MissingChildren <- unlist(lapply(binaryIDs, GetMissingChild, binaryIDs = binaryIDs))
  return(unique(MissingChildren))
}

GetOutsideParent <- function(binaryIDs, nodeBinary){
  # return the parent of nodeBinary if the parent is not in binaryIDs and nodes with the same depth of parent exist.
  # "Root" is not stored in binaryIDs, nodeBinary cannot be "Root"
  minLevel <- min(nchar(binaryIDs))
  if (nchar(nodeBinary)>minLevel & nodeBinary != "Root"){
    parentBinary <- substr(nodeBinary,1, nchar(nodeBinary)-1)
    if (!(parentBinary %in% binaryIDs)) return(parentBinary)
  }
}

GetOutsideParents <- function(binaryIDs){
  # return a vector of parents if the parents is not in binaryIDs and nodes with the same depth of parent exist.
  # use the function GetOutsideParent
  parentsBinary <- lapply(binaryIDs, GetOutsideParent, binaryIDs = binaryIDs)
  return(unique(unlist(parentsBinary)))
}

GetAllAncestorsOne <- function(binaryID){
  # return a vector of all Ancestors of one leaf coded in 01
  # "" not included
  # binaryID = "010", return "0"  "01"
  return(unlist(lapply(1:(nchar(binaryID)-1), function(x,s) substr(s,1,x), s = binaryID)))
}

GetAllUniqueAncestors <- function(binaryIDs){
  # binaryIDs, leaves coded in 01
  # return a sorted vector of all possible parents. "" stands for "Root"
  return(sort(unique(unlist(lapply(binaryIDs,GetAllAncestorsOne)))))
}

String2Vector <- function(s){
  # s is a string like "00 01 11", return a vector ["00","01","11"]
  return(unlist(strsplit(s, split = " ")))
}

ReadInputTree <- function(filename){
  # read input complete tree to data.frame
  return(read.table(filename,sep = '\t',header = TRUE, as.is = TRUE, colClasses = "character", quote = ''))
}

BuildFullPhyloWithLineage <- function(binaryIDs){
  # binaryIDs is a vector of leaves coded in 01
  
  ## n: number of tips. terminal nodes, leaves, node of degree 1
  n <- length(binaryIDs)
  ## m: number of internal nodes.
  m <- n-1
  tip.label <- binaryIDs
  ## optional
  node.label <- c("Root",GetAllUniqueAncestors(binaryIDs))
  all.label <- c(binaryIDs, node.label)
  all.labelID <- 1:length(all.label)
  names(all.labelID) <- all.label
  edgeBinary.tip <- c(binaryIDs, node.label[-1])
  edgeBinary.node <- GetParentsInBinary(edgeBinary.tip)
  edgeID.tip <- all.labelID[edgeBinary.tip]
  edgeID.node <- all.labelID[edgeBinary.node]
  edgeID <- matrix(c(edgeID.node, edgeID.tip), ncol = 2, byrow = FALSE)

  phyloTree = list(edge = edgeID, tip.label = tip.label, Nnode = m, node.label = node.label, Ntip = n, name = "Root")
  class(phyloTree) = "phylo"
  return(phyloTree)
}

## function read the tree comparison file, return a list
ReadLocaTreeComparison <- function(filename,t.count = NULL){
  if (is.null(t.count)){
    f.lines <- readLines(filename, file.info(filename)$size)
    f.length <- length(f.lines)
    t.count <- f.length %/% 13
  } else {
    f.lines <- readLines(filename, 9*t.count)
  }
  alignments.matrix <- matrix(data = f.lines, ncol = 9, byrow = TRUE)
  alignments <- as.data.frame(alignments.matrix)
  colnames(alignments) <- c("id", "Score", "RootS", "RootT", "PruneS", "PruneT", "MatchS", "MatchT", "MatchLength")
  alignments$id <- sub('{','',alignments$id,fixed = TRUE)
  alignments$Score <- sub('Score:', '', alignments$Score, fixed = TRUE)
  alignments$RootS <- sub('RootS:', '', alignments$RootS, fixed = TRUE)
  alignments$RootT <- sub('RootT:', '', alignments$RootT, fixed = TRUE)
  alignments$PruneS <- sub('PruneS:', '', alignments$PruneS, fixed = TRUE)
  alignments$PruneT <- sub('PruneT:', '', alignments$PruneT, fixed = TRUE)
  alignments$MatchS <- strsplit(sub('MatchS:', '', alignments$MatchS, fixed = TRUE), ' ')
  alignments$MatchT <- strsplit(sub('MatchT:', '', alignments$MatchT, fixed = TRUE), ' ')
  alignments$MatchLength <- unlist(lapply(alignments$MatchS, length))
  row.names(alignments) <- alignments$id
  return(alignments)
}

ReadTreeAlignmentsWithTest <- function(filename){
  # read the TreeAlignmentFile, return a data.frame. convert the numeric values.
  f.lines <- readLines(filename)
  f.length <- length(f.lines)
  t.count <- f.length %/% 13
  alignments <- read.csv(col.names = c("id","Score","RootS","RootT","PruneS","PruneT","MatchS", "MatchT","Random","Random.min","Random.max","Random.avg","pvalue"),text = "")
  for (i in 1:t.count){
    l <- i * 9 - 9
    id <- i
    for (j in 1:8){
      alignments[i,j] = f.lines[l+j]
    }
  }
  for (i in 1:t.count){
    l <- 9 * t.count + 4 * i - 4
    alignments[i, "Random"] <- f.lines[l+2]
    alignments[i, "Random.min"] <- strsplit(strsplit(f.lines[l+3], ' ')[[1]][1],':')[[1]][2]
    alignments[i, "Random.max"] <- strsplit(strsplit(f.lines[l+3], ' ')[[1]][2],':')[[1]][2]
    alignments[i, "Random.avg"] <- strsplit(strsplit(f.lines[l+3], ' ')[[1]][3],':')[[1]][2]
    alignments[i, "pvalue"] <- strsplit(strsplit(f.lines[l+3], ' ')[[1]][4],':')[[1]][2]
  }
  
  alignments$id <- sub('{','',alignments$id,fixed = TRUE)
  alignments$Score <- sub('Score:', '', alignments$Score, fixed = TRUE)
  alignments$RootS <- sub('RootS:', '', alignments$RootS, fixed = TRUE)
  alignments$RootT <- sub('RootT:', '', alignments$RootT, fixed = TRUE)
  alignments$PruneS <- strsplit(sub('PruneS:', '', alignments$PruneS, fixed = TRUE), ' ')
  alignments$PruneT <- strsplit(sub('PruneT:', '', alignments$PruneT, fixed = TRUE), ' ')
  alignments$MatchS <- strsplit(sub('MatchS:', '', alignments$MatchS, fixed = TRUE), ' ')
  alignments$MatchT <- strsplit(sub('MatchT:', '', alignments$MatchT, fixed = TRUE), ' ')
  alignments$Random <- strsplit(sub('PValue:', '', alignments$Random, fixed = TRUE), ' ')
  
  alignments$Score <- as.numeric(alignments$Score)
  for (i in 1:t.count){
    alignments$Random[[i]] <- as.numeric(alignments$Random[[i]])
  }
  alignments$Random.min <- as.numeric(alignments$Random.min)
  alignments$Random.max <- as.numeric(alignments$Random.max)
  alignments$Random.avg <- as.numeric(alignments$Random.avg)
  alignments$pvalue <- as.numeric(alignments$pvalue)
  
  return(alignments)
}









