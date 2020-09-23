#' Find the best shape constrained GAM for y~x
#'
#' @param x the explanatory variable
#' @param y the response variable
#' @param plot if T then a plot of the best model (lowest AIC) fit is plotted with pertinent statistics in plot titles
#' @description Uses the package scam to fit shape constrained gam models of y~ s(x). Gaussian link function. All the scam shape
#'              are fitted and the best one is selected. If the smooth term is not significant you may not want to use the model.
#'              It depends on your hypothesis and how you are using the model and residuals.
#' @keywords shape constraint, best model, best fit, AIC
#' @export
#' @examples
#'       best.model.f(USArrests$UrbanPop, USArrests$Rape, plot=T)
best.model.f= function(x, y, plot=F){
  shapes= c("mpi","mpd","cx","cv","micx","micv","mdcx","mdcv")
  modnames= c("Monotonic increasing","Monotonic decreasing","Convex", "Concave", "Monotone increasing convex",
              "Monotone increasing concave","Monotone decreasing concave","Monotone decreasing concave")
  gam.out=list()
  P.AIC= matrix(nrow=length(shapes),ncol=2)
  for (i in 1:length(shapes)){
    gam.out[[i]] <- scam(y ~ s(x,bs=shapes[i]))
    P.AIC[i,1]= summary(gam.out[[i]])$s.table[,4]
    P.AIC[i,2]= AIC(gam.out[[i]])
  }
  names(gam.out)= shapes
  row.names(P.AIC)= shapes

  # choose best model and plot
  best.model.position= match(min(P.AIC[,2]),(P.AIC[,2]))
  bm= shapes[best.model.position]
  bm.AIC= round(P.AIC[best.model.position,2],3)
  bm.P= round(P.AIC[best.model.position,1],3)
  best.model= c(bm,bm.AIC,bm.P)
  names(best.model)=c("shape","AIC","P")
  if (plot){
#    plot(gam.out[[best.model.position]], xlab="Explanatory variable", ylab="Response variable")
#    points(x, y,pch=20)
    plot(x,y,xlab="Explanatory variable", ylab="Response variable",type="n",pch=20)
    newx= data.frame(x=seq(min(x),max(x),length=1000))
    newy= predict(gam.out[[best.model.position]],newdata=newx, se.fit=T)
    confint(newx$x,newy$fit-1.96*newy$se.fit,newy$fit+1.96*newy$se.fit,col="lightgrey")
    lines(newx$x,newy$fit,col="blue",lwd=2)
    points(x,y,pch=20,col="blue")
    title(paste0(modnames[best.model.position]," (",
                 bm,", ",
                 "AIC=", bm.AIC,", ",
                 "P=", bm.P,
                 ")"))
    if (min(P.AIC[,1],na.rm=T)>0.05) title(sub= "there is no model with a significant smoother term",cex.sub=0.8,col.sub="red")
  }
  out=list(dev.expl= summary(gam.out[[best.model.position]])$dev.expl, best.model=best.model, best.scam=gam.out[[best.model.position]])
  out
}

#best.model.f(USArrests$UrbanPop, USArrests$Rape,plot=T)
#best.model.f(tmp$E,tmp$PB,plot=T)

#' Draw polygon confidence intervals
#'
#' @param
#' @keywords helper function
#' @export
#' @examples
#'
confint= function(x,ylow,yhigh,...){
	polygon(x = c(x, rev(x)), y = c(ylow,rev(yhigh)), border = NA,...)
}


# run best model over a dataset with multiple E. Call the best model function for each E, save the deviance explained, rank the
# explanatory power of the variables in terms of their explanatory power of their best fit models. Plot the top five variables in
# terms of explanatory power with their best model fits. Give an option to plot the next 5 (do this as a separate function).

#' Find the best E variable and model to describe another biological variable
#'
#' @param x.matrix a matrix of explanatory (E) variables
#' @param y the response variable
#' @description Fits multiple scam models of response~s(explanatory), choose the best scam for each explanatory variable and then
#'          produces a list of scam objects as well as a vector of deviance explained by each variable. One would likely want an E
#'          variables that explains more of the variance but you also need to consider the amount of data and the shape of the model
#'          fitted.
#'
#'          Presently, the function will fail on the inability to fit a scam model. It needs to be error proofed using tryCatch.
#' @keywords shape constraint, best model, best fit, AIC, deviance explained
#' @export
best.E.f= function(x.matrix, y){
  E.scam= list()
  deviance.explained= vector()
  for (ii in 1:ncol(x.matrix)){
    x= x.matrix[,ii]
    E.scam[[ii]]= best.model.f(x, y, plot=F)
    deviance.explained[ii]= E.scam[[ii]]$dev.expl
  }
  names(E.scam)= names(x.matrix)
  out2= list(dev.expl=round(deviance.explained*100), E.scam=E.scam)
  out2
}



# PB= PB.f(turbot,1990:2000,q=1)$PB
#
# tmp= best.E.f(x.matrix= turbot[,c(5:10,12:30)], y=PB)
# tmp$dev.expl
# xx=turbot[,5]
# taboo=c(1:5,11,31,33,47,48,52,61:753)
#E= turbot[,-taboo]
#tmp= best.E.f(x.matrix= E, y=PB)
#goodmod= match(9,order(tmp$dev.expl,decreasing=T))+5
#names(turbot)[goodmod]
#best.model.f(turbot[,goodmod],PB,plot=T)
# names(turbot)[30]
#
#
#
# x=turbot[,11]
# shite=scam(PB ~ s(x,bs="mdcv"))
# summary(shite)
#
# plot(turbot[,31],PB)
#
#
#
#
#
# best.model.f= function(x, y, plot=F){
#   shapes= c("mpi","mpd","cx","cv","micx","micv","mdcx","mdcv")
#   modnames= c("Monotonic increasing","Monotonic decreasing","Convex", "Concave", "Monotone increasing convex",
#               "Monotone increasing concave","Monotone decreasing concave","Monotone decreasing concave")
#   gam.out=list()
#   P.AIC= matrix(nrow=length(shapes),ncol=2)
#   for (i in 1:length(shapes)){
#     gam.out[[i]]= tryCatch(expr=scam(y ~ s(x,bs=shapes[i])), error="fail")
#     P.AIC[i,1]= tryCatch(expr=summary(gam.out[[i]])$s.table[,4], error="fail")
#     P.AIC[i,2]= tryCatch(expr=AIC(gam.out[[i]]), error="fail")
#   }
#   names(gam.out)= shapes
#   row.names(P.AIC)= shapes
#
#   # choose best model and plot
#   best.model.position= match(min(P.AIC[,2]),(P.AIC[,2]))
#   bm= shapes[best.model.position]
#   bm.AIC= round(P.AIC[best.model.position,2],3)
#   bm.P= round(P.AIC[best.model.position,1],3)
#   best.model= c(bm,bm.AIC,bm.P)
#   names(best.model)=c("shape","AIC","P")
#   if (plot){
# #    plot(gam.out[[best.model.position]], xlab="Explanatory variable", ylab="Response variable")
# #    points(x, y,pch=20)
#     plot(x,y,xlab="Explanatory variable", ylab="Response variable",type="n",pch=20)
#     newx= data.frame(x=seq(min(x),max(x),length=1000))
#     newy= predict(gam.out[[best.model.position]],newdata=newx, se.fit=T)
#     confint(newx$x,newy$fit-1.96*newy$se.fit,newy$fit+1.96*newy$se.fit,col="lightgrey")
#     lines(newx$x,newy$fit)
#     points(x,y,pch=20)
#     title(paste0(modnames[best.model.position]," (",
#                  bm,", ",
#                  "AIC=", bm.AIC,", ",
#                  "P=", bm.P,
#                  ")"))
#     if (min(P.AIC[,1],na.rm=T)>0.05) title(sub= "there is no model with a significant smoother term",cex.sub=0.8,col.sub="red")
#   }
#   out=list(dev.expl= summary(gam.out[[best.model.position]])$dev.expl, best.model=best.model, best.scam=gam.out[[best.model.position]])
#   out
# }
#
