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
               A=NULL, lambda=0, nQTLperInd = 1, 
               # selection parameters
               propSelBetween = .9, propSelWithin =0.9, 
               nGenerations = 30, verbose = FALSE
) 
pmonitor(res0)

## -----------------------------------------------------------------------------
# index for the best solution for trait Value
best=bestSol(res0)["pop","Value"]; best 
# actual solution
res0$M[best,] 
# value and weight for the selected solution 
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
              A=A, lambda= (30*pi)/180 , nQTLperInd = 2, 
              # selection parameters
              propSelBetween = 0.5, propSelWithin =0.5, 
              nGenerations = 20, verbose=FALSE) 

## ----fig.show='hold'----------------------------------------------------------
best = bestSol(res)["pop","Yield"];
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
                A=A, lambda= (20*pi)/180 , nQTLperInd = 100, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 10, verbose=FALSE) 
best = bestSol(res)["pop","GY"]
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
                lambda=(60*pi)/180, nQTLperInd = 80, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 15, verbose = FALSE)

best = bestSol(res)["pop","dummy"]
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
                lambda=(60*pi)/180, nQTLperInd = 80, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 15, verbose = FALSE)
best = bestSol(res)["pop","cov"]
sum(res$M[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
cex <- rep(0.5,nrow(PCWheat))
names(cex) <- rownames(PCWheat)
cex[names(which(res$M[best,]==1))]=2
plot(PCWheat[,1], PCWheat[,2], col = TreeWheat, cex=cex,
     pch = TreeWheat, xlab = "pc1", ylab = "pc2")

## -----------------------------------------------------------------------------
data(DT_technow)
DT <- DT_technow
DT$occ <- 1; DT$occ[1]=0
M <- M_technow
A <- A.mat(M)

Z=with(DT,overlay(dent,flint) )#  Matrix::sparse.model.matrix(~dent-1, data=DT)
rownames(Z) <- DT$hy # needed to link to the QTL matrix

## -----------------------------------------------------------------------------
# regular fitness function
fitnessf <-function (Y, b, d, Q, Z) {
  fit <- Y %*% b - d
  return(fit)
}
# new fitness function with constraint
fitnessf <-function (Y, b, d, Q, Z) {
  X=Q%*%Z[colnames(Q),]
  bad <- as.vector( apply(X,1, function(x){length(which(x > 5))}) ) 
  bad <- which(bad > 0)
  fit <- Y %*% b - d
  if(length(bad) > 0){fit[bad,1]=min(fit[,1])}
  return(fit)
}

## -----------------------------------------------------------------------------
res<-evolafit(formula = c(GY, occ)~hy,
              dt= DT, 
              # constraints: if sum is greater than this ignore
              constraintsUB = c(Inf,50), 
              # constraints: if sum is smaller than this ignore
              constraintsLB= c(-Inf,-Inf),
              # weight the traits for the selection
              traitWeight = c(1,0), 
              # population parameters
              nCrosses = 100, nProgeny = 10, 
              # coancestry parameters
              A=A, lambda= (10*pi)/180 , nQTLperInd = 40, 
              # new fitness function and additional args
              fitnessf = fitnessf, Z=Z,
              # selection parameters
              propSelBetween = 0.5, propSelWithin =0.5, 
              nGenerations = 15, verbose=FALSE) 

best = bestSol(res)["pop","GY"]
xa = (res$M %*% DT$GY)[best,]; xa 
xAx = res$M[best,] %*% A %*% res$M[best,]; xAx 
sum(res$M[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
# check how many times an individual was used in the final crosses
crosses <- data.frame(cross=names(which( res$M[best,] == 1)))
table(unlist(strsplit(crosses$cross,":")))
# check performance of crosses selected
plot(DT$GY, col=as.factor(res$M[best,]), 
       pch=(res$M[best,]*19)+1)

## -----------------------------------------------------------------------------
data("mtcars")
mtcars <- as.data.frame(apply(mtcars,2,scale))
mtcars$inter <- 1
# head(mtcars)
# relationship between the 2 variables
# plot(mpg~hp, data=mtcars)
mod <- lm(mpg~hp, data=mtcars);mod

# create initial QTL effects
a <- seq(-1,1,.1);a
dt <- as.data.frame(expand.grid(a,a))
colnames(dt) <- paste0("alpha",1:ncol(dt))
dt$qtl=paste0("Q",1:nrow(dt))
dt$inter <- rnorm(nrow(dt))
head(dt)

# create n samples equivalent to the number of progeny
# you are planning to simulate (e.g., 1000)
sam <- sample(1:nrow(mtcars),500,replace = TRUE)
y <- mtcars$mpg[sam]
one <- rep(1,length(y))
x <- mtcars$hp[sam]
x2 <- mtcars$hp[sam]^2
X <- cbind(one,x)
plot(x,y)
# Task: linear regression
res0<-evolafit(formula=cbind(inter,alpha1)~qtl, dt= dt,
               # constraints: if greater than this ignore
               constraintsUB = c(Inf,Inf), 
               # constraints: if smaller than this ignore
               constraintsLB= c(-Inf,-Inf), 
               # weight the traits for the selection
               traitWeight = c(1,1), 
               # population parameters
               nCrosses = 50, nProgeny = 10, recombGens = 1, 
               # coancestry parameters
               A=NULL, lambda=0, nQTLperInd = 1, 
               # least MSE function (y - Xb)^2
               fitnessf=function(Y,b,d,Q,x,y){ apply(( (y%*%Jc(500)) - ( X%*%t(Y)) )^2,2,sum) },
               # selection parameters
               propSelBetween = 0.5, propSelWithin =0.5, selectTop=FALSE,
               nGenerations = 20, y=y, x=x, verbose = FALSE
) 

# develop a joint fitness function that uses all traits

pmonitor(res0)
# this time the best solution is the one that minimizes the error
error = ( stan(y) - apply( X*res0$pheno,1,sum ) )^2
best=which(error == min(error))[1]
xa=res0$M[best,] %*% as.matrix(dt[,c("inter","alpha1")]); xa

plot( as.matrix(mtcars[,c("inter","hp")]) %*% t(xa)  , mtcars$mpg,
      main="Correlation between GA-prediction and observed") # GA
plot( (mtcars$hp * mod$coefficients[2] ) + mod$coefficients[1] , mtcars$mpg,
      main="Correlation between lm-prediction and observed") # LM
# Correlation between GA-prediction and observed 
cor( as.matrix(mtcars[,c("inter","hp")]) %*% t(xa)  , mtcars$mpg) 
# Correlation between lm-prediction and observed
cor( (mtcars$hp * mod$coefficients[2] ) + mod$coefficients[1] , mtcars$mpg) # LM


