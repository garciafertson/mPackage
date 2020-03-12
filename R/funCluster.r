
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


funcJaccMat =  getJaccMat(funcMat)
funcJaccNet = makeJaccNet(funcJaccMat)
funcModules = makeModuleList(funcJaccNet)
funcModuleAnnots = annotateModulesByCC(funcJaccNet, funcModules, 3)
writeModuleNetworkAnnot(funcModuleAnnots, funcModNetNodeTableFile, funcModNetEdgeTableFile)








