# jaNAA-2018-04-08.R
# R code for the analysis of NAA samples from Jamaica as described by
# Neiman, Bloch, Galle, and Ferguson in their 2018 SAA poster:
# "Developing Reproducible Methods for Defining and Evaluating 
# Ceramic Compositional Groups Derived from NAA: A Jamaica Example"
# Started:     3.26.2018 FDN
# Last Update: 4.8.2018 FDN 
 

# load the required packages
library(Amelia)
library(RColorBrewer)
library(dplyr)
library(rrcov)
library(maptools)
library(MASS)
library(rgl)
library(cluster)
library (ellipse)


# 1. read in the data from github ####
ja<-read.csv('https://raw.github.com/fneiman/jamaica-naa/master/data/jamNAAData-2018-04-08.csv',
             sep=",",header=TRUE, stringsAsFactors=F)
# check to see what it looks like
str(ja)
# what vars are in what rows?
cbind(1:ncol(ja), colnames(ja))



# 2. Identify missing values #### 
# In this dataset, MURR has coded missing values as 0's -- these are
# observation that are below detection limits.  
# Figure out the proportion observations for element that are below 
# detection limits 
Prop0s<- colSums(ja[,41:73]==0)/nrow(ja)
Prop0s[Prop0s>0]
colSums(ja[,41:73]==0)[Prop0s>0]

# Based on this, it looks like we need to get rid of Ni 
ja$Ni<-NULL
# Now what vars are in what rows?
cbind(1:ncol(ja), colnames(ja))

# Replace 0s with NAs
missingIndex <- ja[,41:72] == 0 
ja[,41:72][missingIndex] <- NA

# do the log transformation
logJa <- with(ja,cbind(ANID,log(ja[,41:72])))

# imputation using the Amelia Package
set.seed(54897678)
a.out<-amelia(logJa, m=1,idvar="ANID",
            max.resample = 1000)
# check to see how the imputation looks
#overimpute(a.out, var="Sr")
#overimpute(a.out, var="U")
# Extract the imputed (and non-imputed) values
logJa <- a.out$imputations[[1]]


# convert to z-scores
logJa[,-1] <- scale(logJa[,-1]) 
# number of site (Estate) names
length(table(ja$Estate))
# colors for Estates
fillCol <- c('white', brewer.pal (12, name='Paired'))


# colors and symbols for plots
ja <- ja %>% mutate(estateCol = case_when(Estate =='Drax Hall' ~ fillCol[1],
                                    Estate =='Good Hope' ~ fillCol[2],
                                    Estate =='Juan de Bolas' ~ fillCol[3],
                                    Estate =='Mona' ~ fillCol[4],
                                    Estate =='Montpelier' ~ fillCol[5],
                                    Estate =='Munchie' ~ fillCol[6],
                                    Estate =='Old Kings House' ~ fillCol[7],
                                    Estate =='Old Naval Dockyard' ~ fillCol[8],
                                    Estate =='Papine' ~ fillCol[9],
                                    Estate =='Saint Peters Church' ~ fillCol[10],
                                    Estate =='Seville' ~ fillCol[11],
                                    Estate =='Stewart Castle' ~ fillCol[12],
                                    Estate =='Thetford' ~ fillCol[13]),
                    coastPch = case_when(Estate =='Drax Hall' ~ 21,
                                     Estate =='Good Hope' ~ 21,
                                     Estate =='Juan de Bolas' ~ 22,
                                     Estate =='Mona' ~ 22,
                                     Estate =='Montpelier' ~ 21,
                                     Estate =='Munchie' ~ 22,
                                     Estate =='Old Kings House' ~ 22,
                                     Estate =='Old Naval Dockyard' ~ 22,
                                     Estate =='Papine' ~ 22,
                                     Estate =='Saint Peters Church' ~ 22,
                                     Estate =='Seville' ~ 21,
                                     Estate =='Stewart Castle' ~ 21,
                                     Estate =='Thetford' ~ 22)
                    )
    
# 2.1 pairwsie scatter plots #####
help(legend)
par(mar=c(5,5,4,2), mfrow= c(1,2))

plot(logJa$Th, logJa$Hf,  
     cex=2,
     pch=ja$coastPch, 
     bg= adjustcolor(ja$estateCol, alpha=.5),
     cex.lab=1.5,
     xlab= 'ln(Thorium)',
     ylab= 'ln(Halfium)')
plot(logJa$Na, logJa$As,  
     cex=2,
     pch=ja$coastPch, 
     bg= adjustcolor(ja$estateCol, alpha=.5),
     cex.lab=1.5,
     xlab= 'ln(Sodium)',
     ylab= 'ln(Arsenic)')
Estates <- unique(cbind(ja$coastPch,ja$estateCol,ja$Estate))
legend('topleft', pch = as.numeric(Estates[,1]), pt.bg=Estates[,2], 
       legend=Estates[,3],
       ncol=2,
       cex=.75,
       pt.cex = 1.25,
       bty='o')


# 3. Outlier detection: Mahalanobis Distance  ####
# quick and dirty plot of Robust vs. classical Mahalanobis distances
mdd <- function(X, alpha, pVal){
  # Function to compute Robust vs. Classical MD
  # Uses MCD method for the robust distances, based on 1-alpha of the data
  # pVal sets the quantile for the oulier cutoff.  
  xBarS<- cov.mcd(X, quantile.used= alpha*nrow(X))
  rob.MD <- mahalanobis(X ,xBarS$center,xBarS$cov)^.5
  MD <-  mahalanobis(X ,colMeans(X),cov(X))^.5
  cutOff <- (qchisq(pVal,df=ncol(X), lower.tail = TRUE))^.5
  return(list(dists=cbind(MD,rob.MD), cutOff=cutOff) )
  }

mddResults <- mdd(X= logJa[,-1], alpha= .9, pVal = (1-(1/448)))    

# fix the plot margins
par(mar=c(5,5,4,2))
par(mfrow= c(1,1))
plot(mddResults$dists[,'MD'], mddResults$dists[,'rob.MD'],  
     cex=2,
     pch=ja$coastPch, 
     bg= adjustcolor(ja$estateCol, alpha=.5),
     cex.lab=1.5,
     xlab= 'Classical Mahalanobis Distance',
     ylab= 'Robust Mahalanobis Distance')
abline(h=mddResults$cutOff, lty=2)
abline(v=mddResults$cutOff, lty=2)
abline(0,1, col=adjustcolor('grey', alpha=.5), lwd=3)
abline(h=15, lty=2)

labelIndex<- (mddResults$dists[,'MD'] > mddResults$cutOff) &
             (mddResults$dists[,'rob.MD'] > mddResults$cutOff)

pointLabel(mddResults$dists[,'MD'][labelIndex], 
           mddResults$dists[,'rob.MD'][labelIndex], 
           logJa$ANID[labelIndex], cex=.75)

Estates <- unique(cbind(ja$coastPch,ja$estateCol,ja$Estate))
legend('bottomright', pch = as.numeric(Estates[,1]), pt.bg=Estates[,2], 
       legend=Estates[,3],
       cex=.75,
       bty='o')

# Pull out the robust MD outliers
mdOutliers<- logJa$ANID[(mddResults$dists[,'rob.MD'] >= 30)]
mdOutliers <- cbind( rep('MD', length(mdOutliers)), as.character(mdOutliers))


# 4. Outlier detection: Robust PCA al la Hubert ####
#  Hubert's RobPCA for 20 dimensions to decide on dimensionality
ctrl1<-CovControlMcd(nsamp=10000)
pca <- PcaHubert(~ As + La + Lu + Nd + Sm + U +  Yb + Ce + Co + Cr + Cs + Eu + Fe + Hf + Rb + Sb + 
                Sc + Sr + Ta + Tb + Th + Zn + Zr + Al + Ba + Ca + Dy + K +  Mn + Na + Ti + V , 
               data = logJa,
               trace = T,
               maxdir = 10000,
               scale = T, 
               alpha = .90, 
               mcd = T,
               k = 20,
               kmax = 50,
               control = ctrl1
)
# first we check out the eigenvalues in a scree plot
# define a function for the broken stick model  
broken.stick <- function(p)
  # Compute the expected values of the broken-stick distribution for 'p' pieces.
  # Example: broken.stick.out.20 = broken.stick(20)
  #             Pierre Legendre, April 2007
{
  result = matrix(0,p,2)
  colnames(result) = c("j","E(j)")
  for(j in 1:p) {
    E = 0
    for(x in j:p) E = E+(1/x)
    result[j,1] = j
    result[j,2] = E/p
  }
  return(result)
}

ProportionVariance<- attr(pca,"eigenvalues")/sum(attr(pca,"eigenvalues"))
cumsum(ProportionVariance)

par(mfrow= c(1,1))
bs<-broken.stick(length(attr(pca,"eigenvalues")))
barplot(ProportionVariance, names.arg=1:length(attr(pca,"eigenvalues")),
        lwd=1, xlab="PC Order",
        pch=21, bg="light grey", cex=1.5,
        ylab= "Proportion of Variance", cex.axis=1.5, cex.lab=1.5)
lines(bs[,1],bs[,2], col="black", lwd=2, lty=2)


#  Hubert's RobPCA for 4 dimensions
pca<-PcaHubert(~ As + La + Lu + Nd + Sm + U +  Yb + Ce + Co + Cr + Cs + Eu + Fe + Hf + Rb + Sb + 
                 Sc + Sr + Ta + Tb + Th + Zn + Zr + Al + Ba + Ca + Dy + K +  Mn + Na + Ti + V , 
               data=logJa,
               trace=T,
               maxdir=10000,
               scale=T, 
               alpha=.90, 
               mcd=T,
               k=4,
               kmax=50,
               control =ctrl1
)

# first we look at the outlier diagnostics
orthDistances<-attr(pca,"od")
scoreDistances<-attr(pca,"sd")
odCut<-attr(pca,"cutoff.od")
sdCut<-attr(pca,"cutoff.sd")
plot(scoreDistances,orthDistances,
     cex=2,
     pch=ja$coastPch, 
     bg= adjustcolor(ja$estateCol, alpha=.5),
     cex.lab=1.5,
     xlim= c(0,11),
     xlab="Score Distances",
     ylab="Orthogonal Distances")
   
abline(v=sdCut,h=odCut, lty=2)

#labelIndex<- (scoreDistances>sdCut)| (orthDistances>odCut)
labelIndex<- (orthDistances > 8)
pointLabel(scoreDistances[labelIndex],orthDistances[labelIndex],
     logJa$ANID[labelIndex],
     cex=.75,
     pos=4,
     adj=1)


Estates <- unique(cbind(ja$coastPch,ja$estateCol,ja$Estate))
legend('bottomright', pch = as.numeric(Estates[,1]), pt.bg=Estates[,2], 
       legend=Estates[,3],
       cex=.75,
       pt.cex = 1.25,
       bty='o')

# pull out ROBPCA outliers : orthogonal distances > 8
orthOutliers <- logJa$ANID[(orthDistances > 8)] 
orthOutliers <- cbind( rep('Orth', length(orthOutliers)), as.character(orthOutliers))


# 5. Remove Outliers from the data #####
# combine the MD and  ROBPCA ouliers -- the former are a subset of the latter
bothOutliers <- rbind(mdOutliers, orthOutliers)
table(bothOutliers[,1], bothOutliers[,2] )
outliers <- unique(bothOutliers[,2])


# create a data frame without outliers (.no)
logJa.no <- filter(logJa, ! ANID %in% outliers)
ja.no <- filter(ja, ! ANID %in% outliers)
# create data frame with only outliers (.oo)
logJa.oo <- filter(logJa,  ANID %in% outliers)
ja.oo <- filter(ja,  ANID %in% outliers)


# 6. Classical PCA w/o outliers  #####
pc1<-prcomp(~ As + La + Lu + Nd + Sm + U +  Yb + Ce + Co + Cr + Cs + Eu + Fe + Hf + Rb + Sb + 
              Sc + Sr + Ta + Tb + Th + Zn + Zr + Al + Ba + Ca + Dy + K +  Mn + Na + Ti + V , 
            data=logJa.no, retx = TRUE, center = TRUE, scale. = T)
# if scale.=TRUE is used, the PCs are extracted from the correlation matrix,
# vs. the covariance matrix 
 #get a summary
summary(pc1)


# first we check out the eigenvalues in a scree plot
ProportionVariance<- pc1$sd^2/(sum(pc1$sd^2))

par(mfrow=c(1,1))
bs<-broken.stick(length(pc1$sd))
barplot(ProportionVariance, names.arg=1:length(pc1$sd),
     lwd=1, xlab="PC Order",
     pch=21, bg="light grey", cex=1.5,
     ylab= "Proportion of Variance", cex.axis=1.5, cex.lab=1.5)
lines(bs[,1],bs[,2], col="black", lwd=2, lty=2)


# plot the scores of the obs -- note we scale these to unit variance
# so scatterplot distances approximate Mahalaobis D
scores<-(predict(pc1))
pcaZScores <- apply(scores, 2, scale)


par(mfrow=c(1,2))
plot(pcaZScores[,1],pcaZScores[,2], 
     pch=ja.no$coastPch, 
     bg= adjustcolor(ja.no$estateCol, alpha=.5),
     col="black", 
     xlim=c(-2,4),
     cex=2, cex.lab= 1.5, 
     xlab="PC 1 (54%)", 
     ylab="PC 2 (12%) ",
     asp=1)
Estates <- unique(cbind(ja$coastPch,ja$estateCol,ja$Estate))
legend('bottomright', pch = as.numeric(Estates[,1]), pt.bg=Estates[,2], 
       legend=Estates[,3],
       ncol=1,
       cex=.75,
       pt.cex = 1,
       bty='o')
#plot the variables on PC1 and PC2
plot(pc1$rotation[,1],pc1$rotation[,2], asp=1,type="n",
     #xlim= c(-.3,.3), ylim=c(-.3,.3),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab='PC 2')
arrows(0, 0, pc1$rotation[,1],pc1$rotation[,2], len=.2, lwd=2, col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.05* pc1$rotation[,1] , 1.05* pc1$rotation[,2] , names(pc1$rotation[,1]),
     col="black", cex=1)




plot(pcaZScores[,1],pcaZScores[,3], 
     pch=ja.no$coastPch, 
     bg= adjustcolor(ja.no$estateCol, alpha=.5),
     cex=2, cex.lab= 1.5, 
     xlab="PC 1 (54%) ", 
     ylab="PC 3 (9%)",
     asp=1)
#text(scores[,1]+1,scores[,3]+1, logJa.no$anid, cex=.5)


#Plots for the variables

plot(pc1$rotation[,1],pc1$rotation[,3], asp=1,type="n",
     xlim= c(-.3,.3), ylim=c(-.4,.5),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab='PC 3')
# Plot arrows: see ?arrows for the syntax 
arrows(0, 0, pc1$rotation[,1],pc1$rotation[,3], len=.2, lwd=2, col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.05* pc1$rotation[,1] , 1.05* pc1$rotation[,3] , 
           names(pc1$rotation[,1]),
           col="black", cex=1)

#3d plot of the classical PCA
library(rgl)
plot3d(scores[,1],
       scores[,2],
       scores[,3],
       xlab="PC 1", 
       ylab="PC 2", 
       zlab ="PC 3",
       cex.lab=1.5,
       col=ja.no$estateCol,
       type="s",
       size=1)

# text3d(scores[,1],
#        scores[,2],
#        scores[,3],
#        logECNr$anid,
#        adj = c(-.5,.5))

# 7. PAM with ouliers removed ####
# the first step is to choose the number of clusters
# clsusgap() will take a while to run!!
gapStatpam <- clusGap(pcaZScores[,1:4], FUN = pam, d.power=2,
                       K.max = 30, B = 50 )
kHat<- maxSE(gapStatpam$Tab[, "gap"], gapStatpam$Tab[, "SE.sim"], method="Tibs2001SEmax")
kHat
 
par(mfrow=c(1,1))
plot(gapStatpam, pch=21, cex=2, bg='grey',
      xlab = 'Number of Clusters (k)',
      cex.lab= 1.5,
      main= '')
abline( v=8, col='red', lty=2)
 
 

# Compute PAM for the chosen K ###########
pam.K <- pam(pcaZScores[,1:4], k=8)

# plot the clusters on the PCs 
# colors for Compositional Groups
fillCol1 <- brewer.pal (12, name='Paired')
par(mfrow=c(1,2))
plot(pcaZScores[,1], pcaZScores[,2], 
     pch=21, 
     bg= adjustcolor(fillCol1[pam.K$clustering], alpha=.75),
     col="black", 
     #xlim= c(-10,10), ylim=c(-6,10),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab="PC 2",
     asp=1)
legend('bottomright', 
       pch = rep(21, max(pam.K$clustering)) ,
       pt.bg = fillCol1[1:max(pam.K$clustering)]  ,
       legend = paste('Trial Group',1:max(pam.K$clustering)),
       pt.cex=2,
       cex=.75,
       bty='o')

plot(pcaZScores[,1], pcaZScores[,3], 
     pch=21, 
     bg= adjustcolor(fillCol1[pam.K$clustering], alpha=.75),
     col="black", 
     #xlim= c(-10,10), ylim=c(-6,10),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab="PC 3",
     asp=1)


plot3d(scores[,1],
       scores[,2],
       scores[,3],
       xlab="PC 1", 
       ylab="PC 2", 
       zlab ="PC 3",
       cex.lab=1.5,
       adjustcolor(fillCol1[pam.K$clustering], alpha=.5),
       type="s",
       size=1.5)


# 8. LDA w/o outliers  #####
library(MASS)
# LDA for ordination
LDAfit <- lda(pam.K$clustering ~ As + La + Lu + Nd + Sm + U +  Yb + Ce + Co + Cr + Cs + Eu + Fe + Hf + Rb + Sb + 
                Sc + Sr + Ta + Tb + Th + Zn + Zr + Al + Ba + Ca + Dy + K +  Mn + Na + Ti + V, 
               data=logJa.no, 
               na.action="na.omit",  method="moment")

# summarize proportion of among-group variance accounted for
propVariance <- LDAfit$svd^2/ sum(LDAfit$svd^2)
par(mfrow=c(1,1))
barplot(propVariance, names.arg=1:length(LDAfit$svd),
        lwd=1, xlab="Function Order",
        pch=21, bg="light grey", cex=1.5,
        ylab= "Proportion of Among-Group Variance", cex.axis=1.5, cex.lab=1.5)

# compute scores of samples for the trial groups
LDAScores<-predict(LDAfit)
# compute scores of the outliers
LDAPreds<-predict(LDAfit, logJa.oo[,-1] )

# plot the scores on LDF1 and LDF2
# note how we use identify to allow us to ID the points with the sample number
par(mfrow=c(1,2))
plot(LDAScores$x[,1],LDAScores$x[,2],  	
     pch= 21,
     bg= adjustcolor(fillCol1[pam.K$cluster], alpha=.75) ,
     col="black", 
     cex=2, cex.lab= 1.5, xlab="Function 1", ylab="Function 2")

points(LDAPreds$x[,1],LDAPreds$x[,2],    
       pch=23, 
       bg=adjustcolor('yellow',alpha=.5) ,
       col="black",
       cex=2) 


#pointLabel(LDAScores$x[,1],LDAScores$x[,2], as.character(pam.K$clustering), cex=.5 )
pointLabel(LDAPreds$x[,1],LDAPreds$x[,2],logJa.oo$ANID, allowSmallOverlap = F, cex=.5 )


legend('topright', ncol=2,
       pch = c(rep(21, max(pam.K$clustering)), 23),
       pt.bg = c(fillCol1[1:max(pam.K$clustering)], 'yellow'),
       legend = c( paste('Trial Group',1:max(pam.K$clustering)), 'Outliers'),
       pt.cex=1.5,
       cex=.75,
       bty='o')


# Do the prediction ellipses
# Define a function called plotEllipses: that will compute the ellipses and put them on the plot
# arguments are: 
# 	X.Var	-the variable to plot on the X axis
#	  Y.Var	-the varaible to plot on the Y axis
#	  Group.Var	-the variable (a factor) with the goups assigments on which the ellpises will be based
#	  PCT	  -The (1-alpha) level for the confidence ellipse, spcified as proportions (e.g. .90)
#		       Uses the ellipse function from the ellipse package.
# 	Method -Use "mcd" for robust covariance estaimtation, "classical" for regular covariance estimation
#		       This uses the cov.rob function from the MASS package


library(ellipse)
PlotEllipses <-function (X.Var, Y.Var, Group.Var, PCT, Method, lineColors){
  XY.Vars<- cbind(X.Var,Y.Var)
  NumGroups<-length(table(Group.Var))
  print(NumGroups)
  for (i in (1:NumGroups)) {
    IthGroupLabel<- names(table(Group.Var))[i]
     GroupIndices<- IthGroupLabel==Group.Var
    IthGroupData<-XY.Vars[GroupIndices,]
    # print(IthGroupData)
    CovRob <- cov.rob(IthGroupData, method=Method )
    CLs<-ellipse(CovRob$cov,centre=CovRob$center, level=PCT, n=100) 
    lines(CLs, col=lineColors[i], lwd=2)
  }
}

PlotEllipses (X.Var=LDAScores$x[,1],Y.Var=LDAScores$x[,2], 
              Group.Var= pam.K$cluster, 
              lineColors = fillCol1,
              PCT=.95, Method="classical")

# plot the outliers only - they seem include a comp.group
# and and a way out point:  JAM016 
plot(LDAPreds$x[,1],LDAPreds$x[,2],    
     pch=23, 
     bg=adjustcolor('yellow',alpha=.75) ,
     col="black",
     cex=2)
pointLabel(LDAPreds$x[,1],LDAPreds$x[,2],logJa.oo$ANID, allowSmallOverlap = F, cex=.5 )


# plot the variables on DF1 and DF2
plot(LDAfit$scaling[,1],LDAfit$scaling[,2],  type="n",
     cex=2, cex.lab= 1.5, xlab="Function 1", ylab="Function 2",
     xlim= c(-2,3), ylim= c(-2,2))
arrows(0, 0, LDAfit$scaling[,1],LDAfit$scaling[,2], len=.05, lwd=2, 
       col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.1* LDAfit$scaling[,1],1.1 * LDAfit$scaling[,2] , 
     rownames(LDAfit$scaling),
     col="black", cex=1)

 
# plot the scores on LDF1 and LDF3
par(mfrow=c(1,2))
plot(LDAScores$x[,1],LDAScores$x[,3],  	
     pch= 21,
     bg= adjustcolor(fillCol1[pam.K$cluster], alpha=.75) ,
     col="black", 
     cex=2, cex.lab= 1.5, xlab="Function 1", ylab="Function 3")

points(LDAPreds$x[,1],LDAPreds$x[,3],    
       pch=23, 
       bg=adjustcolor('yellow',alpha=.5) ,
       col="black",
       cex=2) 
pointLabel(LDAPreds$x[,1],LDAPreds$x[,3],logJa.oo$ANID, allowSmallOverlap = F, cex=.5 )

PlotEllipses (X.Var=LDAScores$x[,1],Y.Var=LDAScores$x[,3], 
              Group.Var= pam.K$cluster, 
              lineColors = fillCol1,
              PCT=.95, Method="classical")


# plot the variables on DF1 and DF3
plot(LDAfit$scaling[,1],LDAfit$scaling[,3],  type="n",
     cex=2, cex.lab= 1.5, xlab="Function 1", ylab="Function 3",
     xlim= c(-2,3), ylim= c(-3,4))
arrows(0, 0, LDAfit$scaling[,1],LDAfit$scaling[,3], len=.05, lwd=2, 
       col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.1* LDAfit$scaling[,1],1.1 * LDAfit$scaling[,3] , 
           rownames(LDAfit$scaling),
           col="black", cex=1)



points(LDAPreds$x[,1],LDAPreds$x[,3],    
       pch=23, 
       bg=adjustcolor('yellow',alpha=.75) ,
       col="black",
       cex=2)



# 3d plot of the LDA scores
library(rgl)
open3d()
plot3d(LDAScores$x[,1],
       LDAScores$x[,2],
       LDAScores$x[,3],
       xlab="DF 1", 
       ylab="DF 2", 
       zlab ="DF 3",
       cex.lab=1.5,
       col= adjustcolor(fillCol1[pam.K$cluster], alpha=.5),
       type="s",
       size=1.5)

spheres3d(LDAPreds$x[,1],
          LDAPreds$x[,2],
          LDAPreds$x[,3],    
          color='yellow', 
          radius=.5) 

# Assess the accuracy of the predictions with LOO CV
LDAfit.cv <- lda(pam.K$clustering ~ As + La + Lu + Nd + Sm + U +  Yb + Ce + Co + Cr + Cs + Eu + Fe + Hf + Rb + Sb + 
                Sc + Sr + Ta + Tb + Th + Zn + Zr + Al + Ba + Ca + Dy + K +  Mn + Na + Ti + V, 
              data=logJa.no, 
              na.action="na.omit",  CV=T, method="moment", prior= rep(1,max(pam.K$clustering))/max(pam.K$clustering) )

ct <- table(pam.K$clustering, LDAfit.cv$class)
ct
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))


# cross tabs for compsitional groups by sites 
table(pam.K$clustering)
groupByEstate <- table(pam.K$clustering, ja.no$Estate)
write.csv(groupByEstate, file ='groupByEstate')


# cross tabs for compsitional groups by Hauser et al's original 2 groups 
table(pam.K$clustering,ja.no$Compositional.Group)




