
# jaLA.R
# R code for the analysis of Bloch's LA samples from Jamaica and England
# includes detection of ouliers, PCA, clustering and LDA without them 
# Started:     3.26.2018 FDN
# Update:      4.4.2018 LB
# Update:      4.4.2018 FDN 


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


# 1. read in the data  ####

ja<-read.csv('C:/Users/fneiman/Desktop/JamaicaNAA/data/JamLAData.csv', sep=",", header=TRUE, stringsAsFactors=F)


ja<-subset(ja, !ja$Location=="MidAtlantic")

# check to see what it looks like
str(ja)
# what vars are in what rows?
cbind(1:ncol(ja), colnames(ja))



# 2. Identify missing values #### 
# In this dataset, MURR has coded missing values as 0's -- these are
# observation that are below detection limits.  
# Figure out the proportion observations for element that are below 
# detection limits 
Prop0s<- colSums(ja[,6:49]==0)/nrow(ja)
Prop0s[Prop0s>0]
colSums(ja[,6:49]==0)[Prop0s>0]

# Based on this, it looks like we need to get rid of Au 
ja$Au<-NULL
# Now what vars are in what rows?
cbind(1:ncol(ja), colnames(ja))

# Replace 0s with NAs
missingIndex <- ja[,6:48] == 0 
ja[,6:48][missingIndex] <- NA

# do the log transformation
logJa <- with(ja,cbind(Sample,log(ja[,6:48])))
                 
# imputation using the Amelia Package
set.seed(54897678)
a.out<-amelia(logJa, m=1,idvar="Sample",
            max.resample = 1000)

# check to see how the imputation looks
#overimpute(a.out, var="Ag")
#overimpute(a.out, var="U")

# Extract the imputed (and non-imputed) values
logJa <- a.out$imputations[[1]]


# convert of z-scores
logJa[,-1] <- scale(logJa[,-1]) 

# number of site (Estate) names
length(table(ja$Location))

# colors for Estates
fillCol <- c('white', brewer.pal (12, name='Paired'))
unique(ja$Location)

ja <- ja %>% mutate(estateCol = case_when(Location =='Drax Hall' ~ fillCol[1],
									Location =='Drax Hall Village' ~ fillCol[1],
                                    Location =='Good Hope' ~ fillCol[2],
                                    #Estate =='Juan de Bolas' ~ fillCol[3],
                                    Location =='Mona Village' ~ fillCol[4],
                                    Location =='Mona GH' ~ fillCol[4],
                                    #Estate =='Montpelier' ~ fillCol[5],
                                    #Estate =='Munchie' ~ fillCol[6],
                                    #Estate =='Old Kings House' ~ fillCol[7],
                                    #Estate =='Old Naval Dockyard' ~ fillCol[8],
                                    Location =='Papine' ~ fillCol[9],
                                    #Estate =='Saint Peters Church' ~ fillCol[10],
                                    Location =='Seville H15' ~ fillCol[11],
                                    Location =='Seville H16' ~ fillCol[11],
                                    Location =='Stewart Castle' ~ fillCol[12],
                                     Location =='SC Village' ~ fillCol[12],
                                    #Estate =='Thetford' ~ fillCol[13]),
                                    Location =='Buckley' ~ fillCol[3],
                                    Location =='Liverpool' ~ fillCol[5],
                                    #Location =='MidAtlantic' ~ fillCol[6],
                                     Location =='London' ~ fillCol[7],
                                     Location =='Staffordshire' ~ fillCol[8]),
                    coastPch = case_when(Location =='Drax Hall' ~ 21,
                    				Location =='Drax Hall Village' ~ 21,
                                     Location =='Good Hope' ~ 21,
                                     #Estate =='Juan de Bolas' ~ 22,
                                     Location =='Mona Village' ~ 22,
                                      Location =='Mona GH' ~ 22,
                                     #Estate =='Montpelier' ~ 21,
                                     #Estate =='Munchie' ~ 22,
                                     #Location =='Old Kings House' ~ 22,
                                     #Location =='Old Naval Dockyard' ~ 22,
                                     Location =='Papine' ~ 22,
                                     #Location =='Saint Peters Church' ~ 22,
                                     Location =='Seville H15' ~ 21,
                                     Location =='Seville H16' ~ 21,
                                      Location =='SC Village' ~ 21,
                                     Location =='Stewart Castle' ~ 21,
                                     #Location =='Thetford' ~ 22
                                      Location =='Buckley' ~ 23,
                                       Location =='Liverpool' ~ 23,
                                         #Location =='MidAtlantic' ~ 23,
                                          Location =='London' ~ 23,
                                          Location =='Staffordshire' ~ 23
                                     
                                     
                                     )
                    )
    

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

mddResults <- mdd(X= logJa[,-1], alpha= .90, pVal = (1-(1/295)))    
#355 with midatlantic
# fix the plot margins
par(mar=c(5,5,4,2))



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
# abline(h=15, lty=2)



labelIndex<- (mddResults$dists[,'rob.MD'] > mddResults$cutOff)

pointLabel(mddResults$dists[,'MD'][labelIndex], 
    mddResults$dists[,'rob.MD'][labelIndex], 
           logJa$Sample[labelIndex], cex=.75)


#pointLabel(mddResults$dists[,'MD'], mddResults$dists[,'rob.MD'], logJa$ANID)

Estates <- unique(cbind(ja$coastPch,ja$estateCol,ja$Location))
legend('topleft', pch = as.numeric(Estates[,1]), pt.bg=Estates[,2], 
       legend=Estates[,3],
       cex=.75,
       bty='o')



# Pull out the robust MD outliers: robust MD > 40
mdOutliers<- logJa$Sample[(mddResults$dists[,'rob.MD'] >= 40)]
mdOutliers <- cbind( rep('MD', length(mdOutliers)), as.character(mdOutliers))
mdOutliers

# 4. Outlier detection: Robust PCA al la Hubert ####

ctrl1<-CovControlMcd(nsamp=10000)


#  Hubert's RobPCA for 20 dimensions to decide on dimensionality
#all LA elements
pca<-PcaHubert(~ Li + Na + Al + Si + K + Ca +  Sc + Ti + V + Cr + Fe + Mn + Co + Ni + Cu + Zn + 
 Rb + Sr + Y + Zr + Nb + Mo + Ag + In + Sn + Sb + Cs + Ba +  La + Ce + Nd + Sm + Eu + Tb + Dy + Yb +  Lu + Hf + Ta + Pb + Bi + Th + U  ,
 #better subset of elements
 #pca<-PcaHubert(~Li + Na + Al + K + Ca + V + Fe + Mn + Co + Ni + Cu + Zn + Rb + Sr + Cs + Ba + Ce + Bi + U, 
               data=logJa,
               trace=T,
               maxdir=10000,
               scale=T, 
               alpha=.95, 
               mcd=T,
               k=20,
               kmax=50,
               control =ctrl1
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
bs<-broken.stick(length(attr(pca,"eigenvalues")))
barplot(ProportionVariance, names.arg=1:length(attr(pca,"eigenvalues")),
        lwd=1, xlab="PC Order",
        pch=21, bg="light grey", cex=1.5,
        ylab= "Proportion of Variance", cex.axis=1.5, cex.lab=1.5)
lines(bs[,1],bs[,2], col="black", lwd=2, lty=2)

#  Hubert's RobPCA for 4 dimensions
pca<-PcaHubert(~ Li + Na + Al + Si + K + Ca +  Sc + Ti + V + Cr + Fe + Mn + Co + Ni + Cu + Zn + 
 Rb + Sr + Y + Zr + Nb + Mo + Ag + In + Sn + Sb + Cs + Ba +  La + Ce + Nd + Sm + Eu + Tb + Dy + Yb +  Lu + Hf + Ta + Pb + Bi + Th + U  ,
#pca<-PcaHubert(~ Li + Na + Al + K + Ca + V + Fe + Mn + Co + Ni + Cu + Zn + Rb + Sr + Cs + Ba + Ce + Bi + U, 
               data=logJa,
               trace=T,
               maxdir=10000,
               scale=T, 
               alpha=.95, 
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
     xlab="Score Distances",
     ylab="Orthogonal Distances")
   
abline(v=sdCut,h=odCut, lty=2)

labelIndex<- (scoreDistances>sdCut)| (orthDistances>odCut)
pointLabel(scoreDistances[labelIndex],orthDistances[labelIndex],
     logJa$Sample[labelIndex],
     cex=.75,
     pos=4,
     adj=1)

Estates <- unique(cbind(ja$coastPch,ja$estateCol,ja$Location))
legend('bottomright', pch = as.numeric(Estates[,1]), pt.bg=Estates[,2], 
       legend=Estates[,3],
       cex=.5,
       bty='o')

# pull out ROBPCA outliers : orthogonal distances > 6
orthOutliers <- logJa$Sample[(orthDistances > 10)] 
orthOutliers <- cbind( rep('Orth', length(orthOutliers)), as.character(orthOutliers))


# 5. Remove Outliers from the data #####
# combine the MD and  ROBPCA ouliers -- the former are a subset of the latter
bothOutliers <- rbind(mdOutliers, orthOutliers)
table(bothOutliers[,1], bothOutliers[,2] )
outliers <- unique(bothOutliers[,2])


# create a data frame without outliers (.no)
logJa.no <- filter(logJa, ! Sample %in% outliers)
ja.no <- filter(ja, ! Sample %in% outliers)

# create data frame with only outliers (.oo)
logJa.oo <- filter(logJa,  Sample %in% outliers)
ja.oo <- filter(ja,  Sample %in% outliers)


# 6. Classical PCA w/o outliers  #####
pc1<-prcomp(~ Li + Na + Al + Si + K + Ca +  Sc + Ti + V + Cr + Fe + Mn + Co + Ni + Cu + Zn + 
 Rb + Sr + Y + Zr + Nb + Mo + Ag + In + Sn + Sb + Cs + Ba +  La + Ce + Nd + Sm + Eu + Tb + Dy + Yb +  Lu + Hf + Ta + Pb + Bi + Th + U  ,
#pc1<-prcomp(~ Li + Na + Al + K + Ca + V + Fe + Mn + Co + Ni + Cu + Zn + Rb + Sr + Cs + Ba + Ce + Bi + U, 
            data=logJa.no, retx = TRUE, center = TRUE, scale. = T)
# if scale.=TRUE is used, the PCs are extracted from the correlation matrix,
# vs. the covariance matrix 
 #get a summary
summary(pc1)


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
ProportionVariance<- pc1$sd^2/(sum(pc1$sd^2))
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


plot(pcaZScores[,1],pcaZScores[,2], 
     pch=ja.no$coastPch, 
     bg= adjustcolor(ja.no$estateCol, alpha=.5),
     col="black", 
     #xlim= c(-10,10), ylim=c(-6,10),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab="PC 2")


#plot the variables on PC1 and PC2
plot(pc1$rotation[,1],pc1$rotation[,2], asp=1,type="n",
     #xlim= c(-.3,.3), ylim=c(-.3,.3),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab='PC 2')
# Plot arrows: see ?arrows for the syntax 
arrows(0, 0, pc1$rotation[,1],pc1$rotation[,2], len=.2, lwd=2, col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.05* pc1$rotation[,1] , 1.05* pc1$rotation[,2] , names(pc1$rotation[,1]),
     col="black", cex=1.5)


plot(pcaZScores[,1],pcaZScores[,3], 
     pch=ja.no$coastPch, 
     bg= adjustcolor(ja.no$estateCol, alpha=.5),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab="PC 3")
#text(scores[,1]+1,scores[,3]+1, logJa.no$anid, cex=.5)


#Plots for the variables

plot(pc1$rotation[,1],pc1$rotation[,3], asp=1,type="n",
     xlim= c(-.3,.3), ylim=c(-.3,.3),
     cex=2, cex.lab= 1.5, xlab="PC 1", ylab='PC 3')
# Plot arrows: see ?arrows for the syntax 
arrows(0, 0, pc1$rotation[,1],pc1$rotation[,3], len=.2, lwd=2, col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
text(1.05* pc1$rotation[,1] , 1.05* pc1$rotation[,3] , names(pc1$rotation[,1]),
     col="black", cex=1.5)


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
# we first loook at average sillouette width

# write a function to do mean silhouette plots for values up to kMax
silAvgWidth <- function(data, maxK){
  avgWidth <-numeric(maxK-1)
  nClusters <- 2:maxK
  for (k in 2:maxK)
  {avgWidth[k-1] <- pam(data, k = k)$silinfo$avg.width}
  return(cbind(nClusters, avgWidth))
}

# for the raw data: call the function and make the plot
silPlot <- silAvgWidth(data = logJa.no[,-1] , maxK=30)

plot(silPlot, type='b', pch=21, cex= 2, bg='grey',
     xlab = 'Number of Clusters (k)',
     ylab = 'Mean Silouette Width',
     cex.lab= 1.5)
abline( v=10, col='red', lty=2)

# For the Rob PCA zscores on PC 1-4: call the function and make the plot
# first recompute without outliers
#  Hubert's RobPCA for 4 dimensions
pca <- PcaHubert(~ Li + Na + Al + Si + K + Ca +  Sc + Ti + V + Cr + Fe + Mn + Co + Ni + Cu + Zn + 
Rb + Sr + Y + Zr + Nb + Mo + Ag + In + Sn + Sb + Cs + Ba +  La + Ce + Nd + Sm + Eu + Tb + Dy + 
  Yb +  Lu + Hf + Ta + Pb + Bi + Th + U, 
               data=logJa.no,
               trace=T,
               maxdir=10000,
               scale=T, 
               alpha=.95, 
               mcd=T,
               k=4,
               kmax=50,
               control =ctrl1
)
robPCAZScores <- apply(getScores(pca), 2, scale) 
silPlot <- silAvgWidth(data = robPCAZScores, maxK=20)

plot(silPlot, type='b', pch=21, cex= 2, bg='grey',
     xlab = 'Number of Clusters (k)',
     ylab = 'Mean Silouette Width',
     cex.lab= 1.5)
abline( v=10, col='red', lty=2)


# for the Classical PCA zscores on PC 1-4: call the function and make the plot
silPlot <- silAvgWidth(data = pcaZScores[,1:4] , maxK=30)

plot(silPlot, type='b', pch=21, cex= 2, bg='grey',
     xlab = 'Number of Clusters (k)',
     ylab = 'Mean Silouette Width',
     cex.lab= 1.5)
abline( v=10, col='red', lty=2)


# for the Classical PCA zscores on PC 1-10: call the function and make the plot
silPlot <- silAvgWidth(data = pcaZScores[,1:10] , maxK=30)

plot(silPlot, type='b', pch=21, cex= 2, bg='grey',
     xlab = 'Number of Clusters (k)',
     ylab = 'Mean Silouette Width',
     cex.lab= 1.5)
abline( v=11, col='red', lty=2)


# now we look at  Tibshiani et al's gap statsitsic 


gapStatpam <- clusGap(logJa.no[,-1], FUN = pam, d.power=2,
                      K.max = 30, B = 50)
kHat<- maxSE(gapStatpam$Tab[, "gap"], gapStatpam$Tab[, "SE.sim"], method="Tibs2001SEmax")
kHat

plot(gapStatpam, pch=21, cex=2, bg='grey',
     xlab = 'Number of Clusters (k)',
     cex.lab= 1.5,
     main = 'Gap Statistics* for Data')
abline( v=6, col='red', lty=2)

apply(getScores(pca)[,1:4], 2, FUN=sd)


gapStatpam <- clusGap(robPCAZScores , FUN = pam, d.power=2,
                      K.max = 30, B = 50)

kHat<- maxSE(gapStatpam$Tab[, "gap"], gapStatpam$Tab[, "SE.sim"], method="Tibs2001SEmax")

plot(gapStatpam, pch=21, cex=2, bg='grey',
     xlab = 'Number of Clusters (k)',
     cex.lab= 1.5,
     main = 'Gap Statistics* for Robust PCs 1:4')
abline( v= 7, col='red', lty=2)


gapStatpam <- clusGap(pcaZScores[,1:4], FUN = pam, d.power=2,
                      K.max = 30, B = 50)
kHat<- maxSE(gapStatpam$Tab[, "gap"], gapStatpam$Tab[, "SE.sim"], method="Tibs2001SEmax")

plot(gapStatpam, pch=21, cex=2, bg='grey',
     xlab = 'Number of Clusters (k)',
     cex.lab= 1.5,
     main = 'Gap Statistics* for Classical PCs 1:4')
abline( v=7, col='red', lty=2)


gapStatpam <- clusGap(pcaZScores[,1:10], FUN = pam, d.power=2,
                      K.max = 30, B = 50)

kHat<- maxSE(gapStatpam$Tab[, "gap"], gapStatpam$Tab[, "SE.sim"],
             method="Tibs2001SEmax")

plot(gapStatpam, pch=21, cex=2, bg='grey',
     main= 'Gap Statistics* for Classical PCs 1:10',
     xlab = 'Number of Clusters (k)',
     cex.lab= 1.5)
abline( v=kHat, col='red', lty=2)




# is there a consensus on k=7. so...

pam.K <- pam(pcaZScores[,1:4], k=7)


# colors for Compositional Groups
fillCol1 <- brewer.pal (12, name='Paired')

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


plot3d(pcaZScores[,1],
       pcaZScores[,2],
       pcaZScores[,3],
       xlab="PC 1", 
       ylab="PC 2", 
       zlab ="PC 3",
       cex.lab=1.5,
       adjustcolor(fillCol1[pam.K$clustering], alpha=.5),
       type="s",
       size=1)


# 8. LDA w/o outliers  #####

library(MASS)
# LDA for ordination
LDAfit <- lda(pam.K$clustering ~ Li + Na + Al + Si + K + Ca +  Sc + Ti + V + Cr + Fe + Mn + Co + Ni + Cu + Zn + 
 Rb + Sr + Y + Zr + Nb + Mo + Ag + In + Sn + Sb + Cs + Ba +  La + Ce + Nd + Sm + Eu + Tb + Dy + Yb +  Lu + Hf + Ta + Pb + Bi + Th + U  ,
#LDAfit <- lda(pam.K$clustering ~ Li + Na + Al + K + Ca + V + Fe + Mn + Co + Ni + Cu + Zn + Rb + Sr + Cs + Ba + Ce + Bi + U, 
               data=logJa.no, 
               na.action="na.omit",  method="moment")


# summarize proportion of among-group variance accounted for

propVariance <- LDAfit$svd^2/ sum(LDAfit$svd^2)

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
plot(LDAScores$x[,1],LDAScores$x[,2],  	
     pch= ja.no$coastPch,
     bg= adjustcolor(fillCol1[pam.K$cluster], alpha=.75) ,
     col="black", 
     cex=2, cex.lab= 1.5, xlab="Function 1", ylab="Function 2")

points(LDAPreds$x[,1],LDAPreds$x[,2],    
       pch=4, 
       col="black",
       cex=2) 


#pointLabel(LDAScores$x[,1],LDAScores$x[,2], as.character(pam.K$clustering), cex=.5 )
pointLabel(LDAPreds$x[,1],LDAPreds$x[,2],logJa.oo$Sample, allowSmallOverlap = F, cex=.5 )
#pointLabel(LDAScores$x[,1],LDAScores$x[,2],logJa$Sample, allowSmallOverlap = F, cex=.5 )

legend('topright', 
       pch = c(23, 22, 21, rep(21, max(pam.K$clustering)), 4),
       pt.bg = c('white', 'white', 'white', fillCol1[1:max(pam.K$clustering)], 'grey'),
       legend = c('England','N.Coast', 'S.Coast', paste('Trial Group',1:max(pam.K$clustering)), 'Outliers'),
       pt.cex=2, ncol=2,
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
              PCT=.90, Method="classical")



# plot the outliers only 
plot(LDAPreds$x[,1],LDAPreds$x[,2],    
     pch=23, 
     bg=adjustcolor('yellow',alpha=.75) ,
     col="black",
     cex=2)

pointLabel(LDAPreds$x[,1],LDAPreds$x[,2], logJa.oo$Sample, allowSmallOverlap = F, cex=.5 )



#plot the variables on DF1 and DF2

plot(LDAfit$scaling[,1],LDAfit$scaling[,2],  type="n",
     cex=2, cex.lab= 1.5, xlab="Function 1", ylab="Function 2",
     xlim= c(-5,5), ylim= c(-4,4))
arrows(0, 0, LDAfit$scaling[,1],LDAfit$scaling[,2], len=.05, lwd=2, 
       col="grey")
abline(h=0,v=0,col="black", lwd=1, lty=2)
pointLabel(1.1* LDAfit$scaling[,1],1.1 * LDAfit$scaling[,2] , 
     rownames(LDAfit$scaling),
     col="black", cex=1)

 
#plot the cases on DF1 and DF3

plot(LDAScores$x[,1],LDAScores$x[,3],  	
     pch= ja.no$coastPch,
     bg= adjustcolor(fillCol1[pam.K$cluster], alpha=.75) ,
     col="black", 
     cex=2, cex.lab= 1.5, xlab="Function 1", ylab="Function 3")

# and any ouliers


points(LDAPreds$x[,1],LDAPreds$x[,3],    
       pch=4, 
       col="black",
       cex=2) 


PlotEllipses (X.Var=LDAScores$x[,1],Y.Var=LDAScores$x[,3], 
              Group.Var= pam.K$cluster, 
              lineColors = fillCol1,
              PCT=.90, Method="classical")



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
       size=1)


spheres3d(LDAPreds$x[,1],
          LDAPreds$x[,2],
          LDAPreds$x[,3],    
          color='black', 
          radius=.5) 


# Assess the accuracy of the prediction with LOO CV

LDAfit.cv <- lda(pam.K$clustering~ Li + Na + Al + Si + K + Ca +  Sc + Ti + V + Cr + Fe + Mn + Co + Ni + Cu + Zn + 
 Rb + Sr + Y + Zr + Nb + Mo + Ag + In + Sn + Sb + Cs + Ba +  La + Ce + Nd + Sm + Eu + Tb + Dy + Yb +  Lu + Hf + Ta + Pb + Bi + Th + U  ,
#LDAfit.cv <- lda(pam.K$clustering ~ Li + Na + Al + K + Ca + V + Fe + Mn + Co + Ni + Cu + Zn + Rb + Sr + Cs + Ba + Ce + Bi + U, 
              data=logJa.no, 
              na.action="na.omit",  CV=T, method="moment", prior= rep(1,max(pam.K$clustering))/max(pam.K$clustering) )

ct <- table(pam.K$clustering, LDAfit.cv$class)

ct
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))


table(pam.K$clustering, ja.no$JAMGROUP)
table(pam.K$clustering, ja.no$Location)

str(ja)




