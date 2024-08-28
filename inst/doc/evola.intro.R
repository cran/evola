## ----setup, include=FALSE-----------------------------------------------------
library(evola)

## -----------------------------------------------------------------------------
set.seed(1)
# Data
Gems <- data.frame(
  Color = c("Red", "Blue", "Purple", "Orange",
            "Green", "Pink", "White", "Black", 
            "Yellow"),
  Weight = round(runif(9,0.5,5),2),
  Value = round(abs(rnorm(9,0,5))+0.5,2)
)
head(Gems)

## -----------------------------------------------------------------------------
# Task: Gem selection. 
# Aim: Get highest combined value.
# Restriction: Max weight of the gem combined = 10. 
res0<-evolafit(cbind(Weight,Value)~Color, dt= Gems,
               # constraints: if greater than this ignore
               constraintsUB = c(10,Inf), 
               # constraints: if smaller than this ignore
               constraintsLB= c(-Inf,-Inf), 
               # weight the traits for the selection
               traitWeight = c(0,1), 
               # population parameters
               nCrosses = 100, nProgeny = 20, recombGens = 1, 
               # coancestry parameters
               A=NULL, lambda=c(0,0), nQTLperInd = 1, 
               # selection parameters
               propSelBetween = .9, propSelWithin =0.9, 
               nGenerations = 30, verbose = FALSE
) 

## -----------------------------------------------------------------------------
best=bestSol(res0)[2]; best # best solution for trait 1
res0$M[best,]
xa = res0$M[best,] %*% as.matrix(Gems[,c("Weight","Value")]); xa

## ----fig.show='hold'----------------------------------------------------------
data(DT_cpdata)
DT <- DT_cpdata
head(DT)

## ----fig.show='hold'----------------------------------------------------------
# get best 20 individuals weighting variance by 0.5
res<-evolafit(cbind(Yield, occ)~id, dt= DT, 
              # constraints: if sum is greater than this ignore 
              constraintsUB = c(Inf,20), 
              # constraints: if sum is smaller than this ignore
              constraintsLB= c(-Inf,-Inf), 
              # weight the traits for the selection
              traitWeight = c(1,0), 
              # population parameters
              nCrosses = 100, nProgeny = 10, 
              # coancestry parameters
              A=A, lambda=c(0.5,0), nQTLperInd = 2, 
              # selection parameters
              propSelBetween = 0.5, propSelWithin =0.5, 
              nGenerations = 20, verbose=FALSE) 

## ----fig.show='hold'----------------------------------------------------------
best = bestSol(res)[1];
sum(res$M[best,]) # total # of inds selected

## ----fig.show='hold'----------------------------------------------------------
pmonitor(res)
plot(DT$Yield, col=as.factor(res$M[best,]), 
     pch=(res$M[best,]*19)+1)


## -----------------------------------------------------------------------------
data(DT_technow)
DT <- DT_technow
DT$occ <- 1; DT$occ[1]=0
M <- M_technow
A <- A.mat(M)
head(DT)

## -----------------------------------------------------------------------------
# run the genetic algorithm
  res<-evolafit(formula = c(GY, occ)~hy, dt= DT, 
                # constraints: if sum is greater than this ignore
                constraintsUB = c(Inf,100), 
                # constraints: if sum is smaller than this ignore
                constraintsLB= c(-Inf,-Inf),
                # weight the traits for the selection
                traitWeight = c(1,0), 
                # population parameters
                nCrosses = 100, nProgeny = 10, 
                # coancestry parameters
                A=A, lambda=c(0.3,0), nQTLperInd = 100, 
                # selection parameters
                propSelBetween = 0.9, propSelWithin =0.9, 
                nGenerations = 10, verbose=FALSE) 
best = bestSol(res)[1]
sum(res$M[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
pmonitor(res)
plot(DT$GY, col=as.factor(res$M[best,]), 
       pch=(res$M[best,]*19)+1)

## -----------------------------------------------------------------------------
data(DT_wheat)
DT <- as.data.frame(DT_wheat)
DT$id <- rownames(DT) # IDs
DT$occ <- 1; DT$occ[1]=0 # to track occurrences
DT$dummy <- 1; DT$dummy[1]=0 # dummy trait
# if genomic
# GT <- GT_wheat + 1; rownames(GT) <- rownames(DT)
# A <-  GT%*%t(GT)
# A <- A/mean(diag(A))
# if pedigree
A <- A_wheat

## -----------------------------------------------------------------------------
##Perform eigenvalue decomposition for clustering
##And select cluster 5 as target set to predict
pcNum=25
svdWheat <- svd(A, nu = pcNum, nv = pcNum)
PCWheat <- A %*% svdWheat$v
rownames(PCWheat) <- rownames(A)
DistWheat <- dist(PCWheat)
TreeWheat <- cutree(hclust(DistWheat), k = 5 )
plot(PCWheat[,1], PCWheat[,2], col = TreeWheat, 
     pch = as.character(TreeWheat), xlab = "pc1", ylab = "pc2")
vp <- rownames(PCWheat)[TreeWheat == 3]; length(vp)
tp <- setdiff(rownames(PCWheat),vp)

## -----------------------------------------------------------------------------
As <- A[tp,tp]
DT2 <- DT[rownames(As),]

## -----------------------------------------------------------------------------
res<-evolafit(cbind(dummy, occ)~id, dt= DT2, 
                # constraints: if sum is greater than this ignore 
                constraintsUB = c(Inf, 100), 
                # constraints: if sum is smaller than this ignore
                constraintsLB= c(-Inf, -Inf), 
                # weight the traits for the selection
                traitWeight = c(1,0), 
                # population parameters
                nCrosses = 100, nProgeny = 10, 
                # coancestry parameters
                A=As,
                lambda=c(1,0), nQTLperInd = 80, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 15, verbose = FALSE)

best = bestSol(res)[1]
sum(res$M[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
cex <- rep(0.5,nrow(PCWheat))
names(cex) <- rownames(PCWheat)
cex[names(which(res$M[best,]==1))]=2
plot(PCWheat[,1], PCWheat[,2], col = TreeWheat, cex=cex,
     pch = TreeWheat, xlab = "pc1", ylab = "pc2")

## -----------------------------------------------------------------------------
DT2$cov <- apply(A[tp,vp],1,mean)

## -----------------------------------------------------------------------------
res<-evolafit(cbind(cov, occ)~id, dt= DT2, 
                # constraints: if sum is greater than this ignore 
                constraintsUB = c(Inf, 100), 
                # constraints: if sum is smaller than this ignore
                constraintsLB= c(-Inf, -Inf), 
                # weight the traits for the selection
                traitWeight = c(1,0), 
                # population parameters
                nCrosses = 100, nProgeny = 10, 
                # coancestry parameters
                A=As,
                lambda=c(1,0), nQTLperInd = 80, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 15, verbose = FALSE)
best = bestSol(res)[1]
sum(res$M[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
cex <- rep(0.5,nrow(PCWheat))
names(cex) <- rownames(PCWheat)
cex[names(which(res$M[best,]==1))]=2
plot(PCWheat[,1], PCWheat[,2], col = TreeWheat, cex=cex,
     pch = TreeWheat, xlab = "pc1", ylab = "pc2")

## -----------------------------------------------------------------------------
# data$motherN <- as.numeric(as.factor(data$mother))
# data$fatherN <- as.numeric(as.factor(data$father))
# fitnessf <- list(motherN= function(motherN,fatherN){
# v = table(c(motherN,fatherN))
# res <- ifelse(any(v > 4), 0, 1)
# return(res)
# })

