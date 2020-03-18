
jacc <- function(x,y) {
  prod = sum(x*y)
  sumx = sum(x)
  sumy = sum(y)
  uni  = sumx + sumy - prod
  return(prod/uni)
}
jacc_comp = compiler::cmpfun(jacc)

copyLowerToUpper <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}
copyLowerToUpper_comp = compiler::cmpfun(copyLowerToUpper)

getJaccMat <- function(mat) {

  jaccMat = matrix(0,dim(mat)[2],dim(mat)[2])
  rownames(jaccMat) = colnames(mat)
  colnames(jaccMat) = colnames(mat)

  combs = t(combn(ncol(mat),2))
  jaccs = apply(combs, 1, function(x) jacc_comp(mat[,x[1]], mat[,x[2]]))
  jaccMat[combs]=jaccs
  diag(jaccMat)=1
  jaccMat = copyLowerToUpper_comp(jaccMat)
  return(jaccMat)
}

makeJaccNet <- function(jaccMat, cutJacc = 0.5) {
  jaccElementsAll = rownames(jaccMat)
  jaccTable = reshape2::melt(jaccMat)
  jaccTable = jaccTable[jaccTable[,3] >= cutJacc & jaccTable[,1] != jaccTable[,2], ]
  jaccTable[,1]=as.character(jaccTable[,1])
  jaccTable[,2]=as.character(jaccTable[,2])

  jaccElementsLeft = unique(c(jaccTable[,1], jaccTable[,2]))
  jaccElementsExcluded = jaccElementsAll[!jaccElementsAll %in% jaccElementsLeft]

  jaccNet = igraph::graph.data.frame(jaccTable[,1:2], directed=F)
  jaccNet = igraph::simplify(jaccNet, remove.multiple=T, remove.loops=T)

  jaccNet = igraph::add_vertices(jaccNet, length(jaccElementsExcluded), attr=list(name=jaccElementsExcluded))

  return(jaccNet)
}


makeModuleList <- function(corNet, debug=F) {
  fc = igraph::cluster_walktrap(corNet)
  cluster = fc$membership
  geneCluster = data.frame(gene=igraph::V(corNet)$name,
                           cluster=cluster,
                           stringsAsFactors=F)
  geneClusterList = split.data.frame(geneCluster, geneCluster$cluster)
  geneClusterList = lapply(geneClusterList, "[[", "gene")
  geneClusterSizes = do.call(c, lapply(geneClusterList, length))

  if (debug) {
    print("... done")
    print(paste("... modularity:", as.character(modularity(fc))))
    print(paste("... no. clusters:", as.character(length(geneClusterList))))
    print(paste("... no. genes of max cluster:", as.character(sort(geneClusterSizes,T)[1])))

  }
  return(geneClusterList)
}


getDescFromKoTerms <- function(koTerms) {
  desc = koDescMap$desc[match(koTerms, koDescMap$KO)]
  return(desc)
}


getHyperP <-function( targetGenes, refGenes, allGenes, debug ) {
  if (debug) {
    nonExistingTargets = as.numeric(table(targetGenes %in% allGenes)["FALSE"])
    print(paste("target genes not existed in reference record:", nonExistingTargets))
  }
  targetGenes = targetGenes[targetGenes %in% allGenes]

  numAll = length(allGenes)
  numOverlap = length( refGenes[refGenes %in% targetGenes] )
  numTrue = length(refGenes)
  numPositive = length(targetGenes)

  out = 1-phyper(numOverlap-1, numTrue, numAll - numTrue, numPositive)

  if (debug) {
    print( "[ hyper-geom. statistics ]")
    print( paste("  all items:", as.character(numAll)) )
    print( paste("  all ref. items:", as.character(numTrue)) )
    print( paste("  all target items:", as.character(numPositive)) )
    print( paste("  all overlapped:", as.character(numOverlap)) )
    print( paste("  p-value:", as.character(out)) )
  }
  return(out)
}


getPhylumName <- function(mgsName, taxo) {
  phylum = taxo[match(mgsName, rownames(taxo)),"phylum"]
  return(phylum)
}
getFamilyName <- function(mgsName, taxo) {
  family = taxo[match(mgsName, rownames(taxo)),"family"]
  return(family)
}

getClassName <- function(mgsName, taxo) {
  class = taxo[match(mgsName, rownames(taxo)),"class"]
  return(class)
}

getOrderName <- function(mgsName, taxo) {
  order = taxo[match(mgsName, rownames(taxo)),"order"]
  return(order)
}


getGenusName <- function(mgsName, taxo) {
  genus = taxo[match(mgsName, rownames(taxo)),"genus"]
  return(genus)
}

getSpeciesName <- function(mgsName, taxo) {
  species = taxo[match(mgsName, rownames(taxo)),"species"]
  return(species)
}



#' @export
getEnrichMentGsByGs <- function(refGeneSet, targetGeneSet, debug=F) {
  allRefGenes = unique.default(do.call(c, refGeneSet))
  allTargetGenes = unique.default(do.call(c, targetGeneSet))
  if (debug) print(length(allRefGenes))
  if (debug) print(length(allTargetGenes))
  if (debug) print(table(allTargetGenes %in% allRefGenes))

  hyperRes = lapply(refGeneSet, function(currRefGenes) {
    hyperPs = lapply(targetGeneSet, function(currTargetGenes){
      currTargetGenes = currTargetGenes[currTargetGenes %in% allRefGenes]
      hyperP = getHyperP( currTargetGenes, currRefGenes, allRefGenes, F )
    })
    return(hyperPs)
  })

  out = lapply(hyperRes, function(x) return(unlist(x)))
  hyperMat = do.call(cbind, out)
  if (debug) print(dim(hyperMat))
  return(hyperMat)
}

#' @export
getEnrichMentGsByList <- function(refGeneSet, targetGenes, debug=F) {
  allRefGenes = unique.default(do.call(c, refGeneSet))
  if (debug) print(table(targetGenes %in% allRefGenes))
  targetGenes = targetGenes[targetGenes %in% allRefGenes]

  hyperPs = lapply(refGeneSet, function(currRefGenes) {
    hyperP = getHyperP( targetGenes, currRefGenes, allRefGenes, F )
    return(hyperP)
  })
  hyperPs = unlist(hyperPs, use.names = F)
  names(hyperPs) = names(refGeneSet)
  return(hyperPs)
}


#' @export
makeFuncCluster <- function(mspKoMat) {
  mspKoJaccMat =  getJaccMat(mspKoMat)
  mspKoJaccNet = makeJaccNet(mspKoJaccMat)
  mspKoModules = makeModuleList(mspKoJaccNet)
  return(mspKoModules)
}

#' @export
getDescriptionOfKeggOrthologFromModules <- function(mspKoModules) {
  return(lapply(mspKoModules, getDescFromKoTerms))
}

#' @export
getCoveredMspByModules <- function(mspKoModules, mspKoMat, cutRatio=0.75, debug=F) {

  mspList = lapply(mspKoModules, function(currFuncs) {
    lenCurrFuncs = length(currFuncs)
    if (length(currFuncs)==1) {
      currMspsByFuncs = mspKoMat[,currFuncs]
    } else {
      currMspsByFuncs = rowSums(mspKoMat[,currFuncs])
    }
    currMspsByFuncs = sort(currMspsByFuncs, decreasing = T)
    currMspsByFuncs = currMspsByFuncs[currMspsByFuncs>=(lenCurrFuncs*cutRatio)]
    return(names(currMspsByFuncs))
  })
  mspTab =  reshape2::melt(mspList)[,c(2,1)]
  colnames(mspTab) = c("Cluster", "MSP")
  mspTab$species = getSpeciesName(mspTab$MSP, gutTaxo)
  return(mspTab)
}






### printing out something during installation ###

#mspKoModules = makeFuncCluster(mspKoMat)
#mspKoDescModules = getDescriptionOfKeggOrthologFromModules(mspKoModules)
#mspByModuleTab = getCoveredMspByModules(mspKoModules, mspKoMat)

###


#set.seed(1)
#koTermsNew = koTerms[sample(1:length(koTerms), 200)]
#mspKoMat = mspKoMat[,colnames(mspKoMat) %in% koTermsNew]
#use_data(mspKoMat, koDescMap, gutTaxoTab, funcModules, internal=T, overwrite = T)
#system.file("")
#print("123")

# save(koDescMap, gutTaxoTab, mspKoMat, file="./data/example.RData")
# save(koDescMap, file="./data/koDescMap.RData")
# save(gutTaxoTab, file="./data/gutTaxoTab.RData")
# save(mspKoMat, file="./data/mspKoMat.Rda")

