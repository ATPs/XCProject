# 20171022
filename <- 'C:\\Users\\ATPs\\OneDrive\\Lab\\YangLab\\YangLabSharedProject_Xiaolong\\Celegans\\Lineages\\fun.alm'
referenceTree <- ReadInputTree(filename)

t1 = Sys.time()
randTree <- GetRandomTreeIn01(671)
print(Sys.time() - t1)



filename <- 'C:\\Users\\ATPs\\OneDrive\\Lab\\YangLab\\YangLabSharedProject_Xiaolong\\Celegans\\Lineages\\20171011CelegansEmbryoLineage.txtl'
alignments <- ReadTreeComparison(filename)



inputTree <- alignments[2,'MatchT']
GenerateRandomTreeFile()

#20171108
##calculate the enrichment of paired genes
filename = 'C:\\Users\\ATPs\\OneDrive\\Lab\\YangLab\\YangLabSharedProject_Xiaolong\\Celegans\\Lineages\\fun.alml'
alignments <- ReadTreeComparison(filename)
