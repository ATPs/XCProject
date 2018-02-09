
library(data.tree)

library(parallel)

F.sigmoid <- function(x, a=100){
  return (2 / (1 + exp(-a*x)) -1)
}

F.geneOF <- function(gene.level){
    if (any(gene.level > 1) || any(gene.level < -1)){
        print("warning! some gene.level is outside the value limit of -1<= gene.level <= 1")
    }
    return (as.numeric(gene.level > 0))
}

F.divide <- function(gene.level){
    if (gene.level[1] > 0){
        L <-  gene.level
        R <- gene.level
        L[1] <- -1
        R[1] <- -1
        if (gene.level[2] >0) L[2] <- -1
        return (data.frame(L,R))
    }
    return (NULL)
}

F.grow <- function(gene.level, gene.network, a = 100){
    return (F.sigmoid(as.vector(gene.network %*% gene.level), a))
}

F.randNetwork <- function(N,K,digits = NULL){
    gene.interactions = rnorm(N*K)
    gene.network = c(gene.interactions, rep(0, N*(N-K)))
    if (is.null(digits))
        return (matrix(sample(gene.network), ncol = N, nrow = N))
    else
        return (round((matrix(sample(gene.network), ncol = N, nrow = N)), digits))
}

F.randInit <- function(N){
    randInit <- c(-1,sample(c(1,-1), size = N-1, replace = TRUE))
    return (randInit)
}

F.develop <- function(gene.init, gene.network, tmax = 50, Dmax = 10, a = 100){
    gene.OF <- F.geneOF(gene.init)
    cell.tree  <- Node$new(name = "Root", DAge = 0, tAge = 0, Lineage="", Gene.level = gene.init, Gene.OF = gene.OF)
    tAge = 0
    DAge = 0
    while (TRUE){
        if (tAge > tmax || DAge > Dmax) break
        for (cell.leaf in cell.tree$leaves){
            cell.divide = F.divide(cell.leaf$Gene.level)
            if (is.null(cell.divide)){
                tAge  <- cell.leaf$tAge + 1
                if (tAge > tmax || DAge > Dmax) break
                name  <- "LR"
                DAge <- cell.leaf$DAge
                Lineage <- cell.leaf$Lineage
                Gene.level <- F.grow(cell.leaf$Gene.level, gene.network = gene.network, a = a)
                Gene.OF <- F.geneOF(Gene.level)
#                 cell.leaf$AddChild(name = name, DAge = DAge, tAge = tAge, 
#                                    Lineage = Lineage, Gene.level = Gene.level, Gene.OF = Gene.OF)
                
#                 cell.leaf$name <- name
                cell.leaf$DAge <- DAge
                cell.leaf$tAge <- tAge
                cell.leaf$Lineage <- Lineage
                cell.leaf$Gene.level <- Gene.level
                cell.leaf$Gene.OF <- Gene.OF
            }
            else {
                tAge <- cell.leaf$tAge + 1
                DAge <- cell.leaf$DAge + 1
                if (tAge > tmax || DAge > Dmax) break
                Lname <- "L"
                LLineage <- paste0(cell.leaf$Lineage,"0")
                Rname <- "R"
                RLineage <- paste0(cell.leaf$Lineage,"1")
                LRinit <- F.divide(cell.leaf$Gene.level)
                LGene.level <- F.grow(LRinit$L, gene.network = gene.network, a = a)
                RGene.level <- F.grow(LRinit$R, gene.network = gene.network, a = a)
                LGene.OF  <- F.geneOF(LGene.level)
                RGene.OF <- F.geneOF(RGene.level)
                cell.leaf$AddChild(name = Lname, DAge = DAge, tAge = tAge, 
                                   Lineage = LLineage, Gene.level = LGene.level, Gene.OF = LGene.OF)
                cell.leaf$AddChild(name = Rname, DAge = DAge, tAge = tAge, 
                                   Lineage = RLineage, Gene.level = RGene.level, Gene.OF = RGene.OF)
            }
        }
    }
    return (cell.tree)
}

F.saveRandomTree <- function(cell.tree, network = NULL, score = "common", 
                             class = "Gene.OF", outfilePrefix = "",silent = FALSE) {
    cell.tree.dfAll <- ToDataFrameTree(cell.tree,"pathString","level","name",'Lineage', 'DAge','tAge',
                                     Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''),
                                     Gene.level = function(node) paste0(node$Gene.level, collapse = ';'),
                                     isLeaf = function(node) node$isLeaf)
    outfile.cell.tree.all <- paste0(outfilePrefix, "cell.tree.all")
    write.table(cell.tree.dfAll, outfile.cell.tree.all, sep="\t", row.names = FALSE)
    if (!silent) print(paste0("write the all details for the tree to file ", outfile.cell.tree.all))
    
    cell.tree.dfLeaf <- ToDataFrameTable(cell.tree,"pathString","level","name",'Lineage', 'DAge','tAge',
                                         Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''))
    cell.tree.Leaf <- cell.tree.dfLeaf[,c("Lineage","name","Gene.OF")]
    outfile.cell.tree.leaves <- paste0(outfilePrefix, "cell.tree.leaves")
    write.table(cell.tree.Leaf, outfile.cell.tree.leaves,sep = "\t", row.names = FALSE, quote = FALSE)
    #remove the last newline symbol in outfile.cell.tree.leaves
    txt <- readChar(outfile.cell.tree.leaves,file.info(outfile.cell.tree.leaves)$size,useByte = TRUE)
    txt <- gsub("\r","",txt)
    if (substr(txt,nchar(txt),nchar(txt)) == "\n")
        cat(substr(txt,1,nchar(txt)-1), file = outfile.cell.tree.leaves)
    if (!silent) print(paste0("write the leaves file for HSA to ", outfile.cell.tree.leaves))
    
    outfile.cell.tree.score <- paste0(outfilePrefix, "cell.tree.score")
    gene.OF <- cell.tree.dfLeaf$Gene.OF
    if (length(unique(gene.OF)) > 2)
        gene.OF.pairs <- combn(unique(gene.OF),2)
    else
        gene.OF.pairs <- NULL
    dfscore <- data.frame(S=c(unique(gene.OF),gene.OF.pairs[1,]), T=c(unique(gene.OF),gene.OF.pairs[2,]))
    if (score == "common")
        dfscore$Score <- apply(dfscore,1,function(x) sum(unlist(strsplit(x["S"],"")) == unlist(strsplit(x["T"],""))))
    else if (score == "same")
        dfscore$Score <- (dfscore$S == dfscore$T)*2
    else if (score == "norm"){
        dfscore$Score <- apply(dfscore,1,function(x) sum(unlist(strsplit(x["S"],"")) == unlist(strsplit(x["T"],""))))
        dfscore$Score <- round((dfscore$Score - mean(dfscore$Score))/sd(dfscore$Score)) + 1
    } else if (!silent) print("invalid value for score, score can be same, common or norm")
        
    write.table(dfscore, outfile.cell.tree.score, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    if (!silent) print(paste0("write the score file to ", outfile.cell.tree.score))
    
    if (!is.null(network)) {
        outfile.network  <- paste0(outfilePrefix, "cell.tree.network")
        write.table(network, outfile.network, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        print(paste0("write the interaction network to file ", outfile.network))
    }
    else if (!silent) print("You need to provide Network in order to save it!")
}

F.generateRandTree <- function(a, N, K, Dmax, tmax, digits = 3){
    gene.init  <- F.randInit(N)
    gene.network <- F.randNetwork(N, K, digits = digits)
    cell.tree <- F.develop(gene.init, gene.network, tmax = tmax, Dmax = Dmax, a = a)
    return(list(gene.init = gene.init, gene.network = gene.network, cell.tree = cell.tree))
}

F.treeCellType <- function(tree1){
    gene.OF <- lapply(tree1$leaves, function(node) paste0(node$Gene.OF,collapse = ""))
    return(length(unique(unlist(gene.OF))))
}

F.generateRandTreeWithFilter  <- function(a, N, K, Dmax, tmax, digits = 3, 
                                          leafTypeMin = 6, leafCountMin = 500, tryMax = 100){
    for (i in 1:tryMax){
        celltree <- F.generateRandTree(a = a, N = N, K=K, Dmax = Dmax, tmax = tmax, digits = digits)
        leafTypeCount <- F.treeCellType(celltree[["cell.tree"]])
        leafCount <- celltree[["cell.tree"]]$leafCount
        celltree$leafTypeCount <- leafTypeCount
        celltree$leafCount <- leafCount
        if (leafCount >= leafCountMin && leafTypeCount >= leafTypeMin){
            return(celltree)
        }
    }
    print(paste0("tried ", tryMax, " times. Still No good. return NULL"))
    return(NULL)
}

N = 8
K = 4
Dmax = 10
tmax = 50
a = 100

gene.init <-  F.randInit(N)
gene.network <- F.randNetwork(N,K)
gene.OF <- F.geneOF(gene.init)
cell.tree  <- Node$new(name = "Root", DAge = 0, tAge = 0, Lineage="1", Gene.level = gene.init, Gene.OF = gene.OF)

for (cell.leaf in cell.tree$leaves) print(cell.leaf)

gene.network

gene.init

gene.init <-  F.randInit(N)
gene.network <- F.randNetwork(N,K)
tree1 <- F.develop(gene.init = gene.init, gene.network = gene.network,tmax = 10, Dmax = 5)
print(tree1,"Lineage","DAge","tAge",Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''))
tree1$leafCount

ToDataFrameTree(tree1,"pathString","level","name",'Lineage', 'DAge','tAge',
                 Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''),
                 Gene.level = function(node) paste0(node$Gene.level, collapse = ';'),
                 isLeaf = function(node) node$isLeaf)

tree1.df <- ToDataFrameTree(tree1,"pathString","level","name",'Lineage', 'DAge','tAge',
                 Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''),
                 Gene.level = function(node) paste0(node$Gene.level, collapse = ';'),
                 isLeaf = function(node) node$isLeaf)

getwd()
# write.table(tree1.df,"test.cell.tree.all",sep = "\t",row.names = FALSE)

ToDataFrameTable(tree1,"pathString","level","name",'Lineage', 'DAge','tAge',
                 Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''))

F.treeCellType(tree1)

tree.dfLeaf <- ToDataFrameTable(tree1,"pathString","level","name",'Lineage', 'DAge','tAge',
                 Gene.OF = function(node) paste0(node$Gene.OF, collapse = ''))
tree.dfLeaf[,c("Lineage","name","Gene.OF")]

# F.saveRandomTree(tree1, network = gene.network,score = "same")

tree1

df.settingD12 = data.frame(a = c(1,10,100,100,100,100,100), 
                        N = c(16,16,16,16,16,8,32),
                        K = c(4,4,2,4,8,4,4),
                        Dmax = rep(12,7),
                        tmax = rep(50,7))

df.settingD12

leafTypeMin = 6
leafCountMin = 500

s <- df.settingD12[1,]
s

celltrees = list()
for (i in row.names(df.settingD12)) {
    s <- df.settingD12[i,]
    celltrees[[1]] = F.generateRandTree(a = s$a, N = s$N, K = s$K, Dmax = s$Dmax, tmax = s$tmax, digits = 4)
    print(s)
    break
}
celltrees[[1]][["cell.tree"]]$leafCount


celltrees[[1]][["cell.tree"]]$Gene.OF

F.treeCellType(celltrees[[1]][["cell.tree"]])

tree1 = celltrees[[1]][["cell.tree"]]

celltrees <- lapply(1:2,
                    function(x)F.generateRandTreeWithFilter(a = s$a, N = s$N, K = s$K, Dmax = s$Dmax, 
                                        tmax = s$tmax, digits = 4,
                                        leafTypeMin = leafTypeMin,
                                        leafCountMin = leafCountMin, tryMax = 100))

for (tr in celltrees) print(c(tr$leafCount, tr$leafTypeCount))

tr <- celltrees[[1]]
# F.saveRandomTree(cell.tree = tr$cell.tree, 
#                  score = "norm", outfilePrefix = "X:\\YangLab\\2018simulation\\", silent = TRUE)



# fileTest <- "C:\\Users\\ATPs\\Documents\\2018simulation\\a1N16K4\\test.txt"
# fileTestWrite <- file(fileTest,"w")
# writeLines(text = c(paste0("ID=1;a=1;N=16;K=4;L=",tr$leafCount,";T=",tr$leafTypeCount),
#                    paste(tr$gene.init,collapse = " "),
#                    paste(as.vector(tr$gene.network), collapse = " ")), con = fileTestWrite)
# close(fileTestWrite)

tr$leafTypeCount

s <- df.settingD12[3,]

folder <- paste0("X:\\YangLab\\2018simulation\\a",s$a,"N",s$N,"K",s$K,"\\")

fileInfo <- paste0("X:\\YangLab\\2018simulation\\a",s$a,"N",s$N,"K",s$K,".info")
fileInfoWrite <- file(fileInfo,"w")
for (n in 1:1000){
    tr <- F.generateRandTreeWithFilter(a = s$a, N = s$N, K = s$K, Dmax = s$Dmax, 
                                        tmax = s$tmax, digits = 4,
                                        leafTypeMin = leafTypeMin,
                                        leafCountMin = leafCountMin, tryMax = 100)
    filenamePrefix <- paste0(folder,"a",s$a,"N",s$N,"K",s$K,"ID",n)
    F.saveRandomTree(cell.tree = tr$cell.tree, network = NULL, score = "norm", class = "Gene.OF", 
                     outfilePrefix = filenamePrefix, silent = TRUE)
    writeLines(text = c(paste0("ID=",n,";a=",s$a,";N=",s$N,";K=",s$K,";L=",tr$leafCount,";T=",tr$leafTypeCount),
                   paste(tr$gene.init,collapse = " "),
                   paste(as.vector(tr$gene.network), collapse = " ")), con = fileInfoWrite)
}
close(fileInfoWrite)

