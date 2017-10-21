
library(zoo)

###This file contains functions that allow for quick plotting of ROC curves!###
### .pfd, .png are unnecessary because we save in .svg => easily transferable to other formats!###

addTPcol <- function(result, pattern="UPS", TPcolname="ups"){
  result <- lapply(result, function(x){
    attr <- attr(x,"MSqRob_fdrtool")
    x <- cbind(x,grepl(pattern,rownames(x)))
    names(x)[ncol(x)] <- TPcolname
    attr(x,"MSqRob_fdrtool") <- attr
    return(x)
  })
  return(result)
}


makeSummary <- function(resultobject, upratio, qval="qval", level=0.05, TPcol="ups", estimate="estimate", addEstimate=TRUE){

  n <- length(resultobject)

  summary_object <- data.frame(bias0=rep(NA, n),bias1=rep(NA, n),sd0=rep(NA, n),sd1=rep(NA, n),mad0=rep(NA, n),mad1=rep(NA, n),RMSE0=rep(NA, n),RMSE1=rep(NA, n),TP5proc=rep(NA, n),FP5proc=rep(NA, n),TN5proc=rep(NA, n),FN5proc=rep(NA, n))

  for(i in 1:n)
  {

    significante <- resultobject[[i]][which(resultobject[[i]][,qval]<=level),]
    nonsignificante <- resultobject[[i]][which(resultobject[[i]][,qval]>level),]

    waarden3 <- na.omit(resultobject[[i]])

    if(isTRUE(addEstimate)){
    #RMSE
    summary_object$RMSE0[i] <- sqrt(mean(((subset(waarden3, waarden3[,TPcol]==0)[,estimate])-0)^2, na.rm=TRUE))
    #bias:
    summary_object$bias0[i] <- mean((subset(waarden3, waarden3[,TPcol]==0)[,estimate])-0, na.rm=TRUE)
    #standard deviation:
    summary_object$sd0[i] <- sd((subset(waarden3, waarden3[,TPcol]==0)[,estimate]), na.rm=TRUE)
    #mad:
    summary_object$mad0[i] <- mad((subset(waarden3, waarden3[,TPcol]==0)[,estimate]), na.rm=TRUE)

    #RMSE
    summary_object$RMSE1[i] <- sqrt(mean(((subset(waarden3, waarden3[,TPcol]==1)[,estimate])-upratio[i])^2, na.rm=TRUE))
    #bias:
    summary_object$bias1[i] <- mean((subset(waarden3, waarden3[,TPcol]==1)[,estimate])-upratio[i], na.rm=TRUE)
    #standard deviation:
    summary_object$sd1[i] <- sd((subset(waarden3, waarden3[,TPcol]==1)[,estimate]), na.rm=TRUE)
    #standard deviation:
    summary_object$mad1[i] <- mad((subset(waarden3, waarden3[,TPcol]==1)[,estimate]), na.rm=TRUE)
    }

    #TP:
    summary_object$TP5proc[i] <- sum(significante[,TPcol])
    #FP:
    summary_object$FP5proc[i] <- length(significante[,TPcol])-sum(significante[,TPcol])
    #FN:
    summary_object$FN5proc[i] <- sum(nonsignificante[,TPcol])
    #TN:
    summary_object$TN5proc[i] <- length(nonsignificante[,TPcol])-sum(nonsignificante[,TPcol])

  }

  rownames(summary_object) <- names(resultobject)

  return(summary_object)
}


plotROC=function(resultobjectlist, summary_objectlist, proteins, truth, colors, pAUC_cutoff=0.1, directory=getwd(), filenames="ROC_curves", signifcol="signif", main_names=NULL, relative=TRUE, plotSVG=FALSE){

  AUC <- matrix(nrow=ncol(L), ncol=length(resultobjectlist))
  pAUC <- matrix(nrow=ncol(L), ncol=length(resultobjectlist))

  dimnames(AUC)=list(colnames(L),names(resultobjectlist))
  dimnames(pAUC)=list(colnames(L),names(resultobjectlist))

  #Order everything by ordering of truth
  proteins2 <- proteins[names(truth)]
  data <- MSqRob::getData(proteins2)
  data2 <- adjustNames(data, colnames(data[[1]]))

  TPtot <- vector("list", ncol(L))
  FPtot <- vector("list", ncol(L))

  if(is.null(main_names)){main_names <- colnames(L)}

  for(i in 1:ncol(L)){

    #DA proteins with at least one peptide in each condition that is compared
    namesComp <- names(which(L[,colnames(L)[i]]!=0))

    allnamespresent <- vapply(data2, function(y){sum(unlist(lapply(y, function(x){namesComp %in% x})))}, 0)==length(namesComp)

    if(!isTRUE(relative)){allnamespresent <- rep(TRUE,length(allnamespresent))}

    #HUMAN:
    TPtot[[i]] <- sum(allnamespresent & truth)

    #YEAST:
    FPtot[[i]] <- sum(allnamespresent & !truth)

    if(isTRUE(plotSVG)){grDevices::svg(file.path(directory,paste0(filenames,i,".svg")), width=10, height=10)}

    for(j in 1:length(resultobjectlist)){

      frame <- resultobjectlist[[j]][[i]]
      #Sort frame by names(truth)
      #frame <- frame[names(truth),]
      truth2 <- truth[rownames(frame[!is.na(frame[,signifcol]),])]

      TP <- c(0,cumsum(truth2))/TPtot[[i]]
      FP <- c(0,cumsum(!truth2))/FPtot[[i]]

      if(j==1){
        plot(TP~FP, main=main_names[i], type="l", xlim=c(0,1), ylim=c(0,1), xlab="True positive rate", ylab="False positive rate",
             frame.plot=FALSE, cex.main=2, cex.lab=2, cex.axis=2, cex=2, lwd=2, pch=1, col=colors[j])
        } else{
          lines(TP~FP, col=colors[j], cex.main=2, cex.lab=2, cex.axis=2, cex=2, lwd=2, pch=1)
          }

      points(
        (summary_objectlist[[j]][i,"FP5proc"])/FPtot[[i]],
        (summary_objectlist[[j]][i,"TP5proc"])/TPtot[[i]],
        cex=2, lwd=2, col=colors[j])

      AUC[i,j] <- sum(diff(FP)*rollmean(TP,2))
      pAUC[i,j] <- sum(diff(FP[FP<pAUC_cutoff])*rollmean(TP[FP<pAUC_cutoff],2))

    }

    if(isTRUE(plotSVG)){dev.off()}

    }
  return(list(AUC=AUC, pAUC=pAUC))
}

plotMA_MSqRob=function(resultobject, proteins, quant_name="quant_value", qname="qval", level=0.05, addNames=FALSE, colors=c("green","red"), directory=getwd(), filenames="MA_plots", main_names=NULL, plotSVG=FALSE){

  if(is.null(main_names)){main_names <- names(resultobject)}

  if(isTRUE(plotSVG)){grDevices::svg(file.path(directory,paste0(filenames,i,".svg")), width=10, height=10)}

  for(i in 1:length(resultobject)){
  resultsSorted <- resultobject[[i]][as.character(getAccessions(proteins)),]

  A = vapply(MSqRob::getData(proteins), function(x) mean(x[[quant_name]]),0)
  names(A) = as.character(getAccessions(proteins))

  #Extract the log2 fold change estimate from the contrastHeLa object
  M = resultsSorted[,"estimate"]
  names(M) = as.character(getAccessions(proteins))

  plot(A,M, xlab="Average log2-expression", ylab="log2-fold-change", main=main_names[i])

  #Plots a horizontal line at 0
  abline(h=0, col=colors[1], lwd=2)

  #Plot the significant proteins in red
  topProt = as.character(getAccessions(proteins)[which(resultsSorted[,qname]<level)])
  points(A[topProt], M[topProt],col=colors[2])

  #Add the protein labels
  if(isTRUE(addNames)){
  text(A[topProt],M[topProt],labels=topProt,col=colors[2],pos=4)
  }

  if(isTRUE(plotSVG)){dev.off()}

  }

}

# plotVolcano_MSqRob <- function(){
#   #Make the volcano plot
#   plot(contrastHeLa[,"estimate"],-log10(contrastHeLa[,"pval"]), cex=0.3, main="MSqRob",xlab="log2 FC",ylab="-log10(p-value)")
#   #Draw a red circle around the significant proteins
#   points(contrastHeLa[topProt,"estimate"],-log10(contrastHeLa[topProt,"pval"]), col=2)
#   #Add labels
#   text(contrastHeLa[topProt,"estimate"],-log10(contrastHeLa[topProt,"pval"]),labels=topProt,col="red",pos=2)
# }
