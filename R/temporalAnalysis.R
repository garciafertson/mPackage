

#' @export
checkDetectionLimit <- function(mgsMat, percent=0.01) {
  mgsMatScaled = mgsMat*10^9
  mgsMatScaledLof = log2(mgsMatScaled)
  mgsMatScaledLof[mgsMatScaledLof == -Inf] = 0

  mgsMatScaledNonZero = mgsMatScaled[mgsMatScaled!=0]
  mgsMatScaledNonZero = log2(mgsMatScaledNonZero)
  fit.gamma = fitdist(mergeMatScaledNonZero, "gamma")
  plot(fit.gamma)
  fit.gamma$estimate

  q01 = qgamma(percent, fit.gamma$estimate["shape"], fit.gamma$estimate["rate"])
  #4.042488

  qCutOff =(2^(q01))/10^9

  return(qCutOff)

}

#' @export
getInflowStats <- function(metaTab, mgsMat, detectionLimit=1.647821e-08) {

  mgsMat[mgsMat < detectionLimit] = 0 #suppressing under limit samples
  metaTab=metaTab[order(metaTab$time),]
  uniqueSubjects = unique(metaTab$metadataId)

  inflowStatTab = do.call(rbind, sapply(rownames(mgsMat), function(currMgs){

    markovSeqOfMsp = (sapply(uniqueSubjects, function(subj) {
      currSamples = metaTab$sampleId[metaTab$metadataId==subj]
      currMgsMat = mgsMat[,match(currSamples, colnames(mgsMat))]
      currMgsSubjVec = (currMgsMat[currMgs,] >0)*1
      currMgsSubjVec = as.character(currMgsSubjVec)

    }, simplify = F))

    myMarkovFit <- markovchain::markovchainFit(data=markovSeqOfMsp,confidencelevel = .9,method = "mle")
    myMarkovMat = myMarkovFit$estimate@transitionMatrix
    myMarkovUpper = myMarkovFit$upperEndpointMatrix
    myMarkovLower = myMarkovFit$lowerEndpointMatrix
    myMarkovSE = myMarkovFit$standardError

    if(dim(myMarkovMat)[1]==1) {
      rname = rownames(myMarkovMat)
      if (rname == "0") {
        colonized = 0
        colonizedUpper = NA
        colonizedLower = NA
        colonizedSE = NA
        result = c(colonized, colonizedUpper, colonizedLower, colonizedSE)
        return(result)
      }

      if (rname == "1") {
        colonized = 1
        colonizedUpper = NA
        colonizedLower = NA
        colonizedSE = NA
        result = c(colonized, colonizedUpper, colonizedLower, colonizedSE)

        return(result)
      }
    }
    else {
      colonized = myMarkovMat["0","1"]
      colonizedUpper = myMarkovUpper["0","1"]
      colonizedLower = myMarkovLower["0","1"]
      colonizedSE = myMarkovSE["0","1"]
      result = c(colonized, colonizedUpper, colonizedLower, colonizedSE)

      return(result)
    }

  }, simplify = F))

  rownames(inflowStatTab) = rownames(mgsMat)
  colnames(inflowStatTab) = c("inflow","upper","lower","SE")
  return(inflowStatTab)

}

#' @export
getOutflowStats <- function(metaTab, mgsMat, detectionLimit=1.647821e-08) {

  mgsMat[mgsMat < detectionLimit] = 0 #suppressing under limit samples
  metaTab=metaTab[order(metaTab$time),]
  uniqueSubjects = unique(metaTab$metadataId)

  outflowStatTab = do.call(rbind, sapply(rownames(mgsMat), function(currMgs){

    markovSeqOfMsp = (sapply(uniqueSubjects, function(subj) {
      currSamples = metaTab$sampleId[metaTab$metadataId==subj]
      currMgsMat = mgsMat[,match(currSamples, colnames(mgsMat))]
      currMgsSubjVec = (currMgsMat[currMgs,] >0)*1
      currMgsSubjVec = as.character(currMgsSubjVec)

    }, simplify = F))

    myMarkovFit<-markovchainFit(data=markovSeqOfMsp,confidencelevel = .9,method = "mle")
    myMarkovMat = myMarkovFit$estimate@transitionMatrix
    myMarkovUpper = myMarkovFit$upperEndpointMatrix
    myMarkovLower = myMarkovFit$lowerEndpointMatrix
    myMarkovSE = myMarkovFit$standardError

    if(dim(myMarkovMat)[1]==1) {
      rname = rownames(myMarkovMat)
      #if (rname == "0") return(Inf)
      #if (rname == "1") return(0)
      if (rname == "0") {
        resistant = 1
        resistantUpper = NA
        resistantLower = NA
        resistantSE = NA
        result = c(resistant, resistantUpper, resistantLower, resistantSE)
        return(result)
      }
      if (rname == "1") {
        resistant = 0
        resistantUpper = NA
        resistantLower = NA
        resistantSE = NA
        result = c(resistant, resistantUpper, resistantLower, resistantSE)
        return(result)
      }
    }
    else {
      resistant = myMarkovMat["1","0"]
      resistantUpper = myMarkovUpper["1","0"]
      resistantLower = myMarkovLower["1","0"]
      resistantSE = myMarkovSE["1","0"]
      result = c(resistant, resistantUpper, resistantLower, resistantSE)

      #preserved = myMarkovMat["0","0"]
      #relativeRisk = resistant/ preserved
      #return(relativeRisk)
      return(result)

    }

  }, simplify = F))
  rownames(outflowStatTab) = rownames(mgsMat)
  colnames(outflowStatTab) = c("outflow","upper","lower","SE")
  return(outflowStatTab)

}




