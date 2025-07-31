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
               b = c(0,1), 
               # population parameters
               nCrosses = 100, nProgeny = 20, recombGens = 1, 
               # coancestry parameters
               D=NULL, lambda=0, nQtlStart = 1, 
               # selection parameters
               propSelBetween = .9, propSelWithin =0.9, 
               nGenerations = 15, verbose = FALSE
) 
pmonitor(res0)

## -----------------------------------------------------------------------------
# index for the best solution for trait Value
best=bestSol(res0$pop)[,"Value"]; best 
# actual solution
Q <- pullQtlGeno(res0$pop, simParam = res0$simParam, trait=1); Q <- Q/2
Q[best,] 
# value and weight for the selected solution 
qa = Q[best,] %*% as.matrix(Gems[,c("Weight","Value")]); qa

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
              b = c(1,0), 
              # population parameters
              nCrosses = 100, nProgeny = 10, 
              # coancestry parameters
              D=A, lambda= (30*pi)/180 , nQtlStart = 2, 
              # selection parameters
              propSelBetween = 0.5, propSelWithin =0.5, 
              nGenerations = 15, verbose=FALSE) 

## ----fig.show='hold'----------------------------------------------------------
Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
best = bestSol(res$pop)[,"Yield"];
sum(Q[best,]) # total # of inds selected

## ----fig.show='hold'----------------------------------------------------------
pmonitor(res)
plot(DT$Yield, col=as.factor(Q[best,]), 
     pch=(Q[best,]*19)+1)


## -----------------------------------------------------------------------------
data(DT_technow)
DT <- DT_technow
DT$occ <- 1; DT$occ[1]=0
M <- M_technow
D <- A.mat(M)
head(DT)

## -----------------------------------------------------------------------------
# silent for CRAN checks restriction on vignettes time
# run the genetic algorithm
#   res<-evolafit(formula = c(GY, occ)~hy, dt= DT, 
#                 # constraints: if sum is greater than this ignore
#                 constraintsUB = c(Inf,100), 
#                 # constraints: if sum is smaller than this ignore
#                 constraintsLB= c(-Inf,-Inf),
#                 # weight the traits for the selection
#                 b = c(1,0), 
#                 # population parameters
#                 nCrosses = 100, nProgeny = 10, 
#                 # coancestry parameters
#                 D=D, lambda= (20*pi)/180 , nQtlStart = 100, 
#                 # selection parameters
#                 propSelBetween = 0.5, propSelWithin =0.5, 
#                 nGenerations = 15, verbose=FALSE) 
# 
# Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
# best = bestSol(res$pop)[1,1]
# sum(Q[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
# pmonitor(res)
# plot(DT$GY, col=as.factor(Q[best,]), 
#        pch=(Q[best,]*19)+1)

## -----------------------------------------------------------------------------
data(DT_wheat)
DT <- as.data.frame(DT_wheat)
DT$id <- rownames(DT) # IDs
DT$occ <- 1; DT$occ[1]=0 # to track occurrences
DT$dummy <- 1; DT$dummy[1]=0 # dummy trait
# if genomic
# GT <- GT_wheat + 1; rownames(GT) <- rownames(DT)
# D <-  GT%*%t(GT)
# D <- D/mean(diag(D))
# if pedigree
D <- A_wheat

## -----------------------------------------------------------------------------
##Perform eigenvalue decomposition for clustering
##And select cluster 5 as target set to predict
pcNum=25
svdWheat <- svd(D, nu = pcNum, nv = pcNum)
PCWheat <- D %*% svdWheat$v
rownames(PCWheat) <- rownames(D)
DistWheat <- dist(PCWheat)
TreeWheat <- cutree(hclust(DistWheat), k = 5 )
plot(PCWheat[,1], PCWheat[,2], col = TreeWheat, 
     pch = as.character(TreeWheat), xlab = "pc1", ylab = "pc2")
vp <- rownames(PCWheat)[TreeWheat == 3]; length(vp)
tp <- setdiff(rownames(PCWheat),vp)

## -----------------------------------------------------------------------------
As <- D[tp,tp]
DT2 <- DT[rownames(As),]

## -----------------------------------------------------------------------------
res<-evolafit(cbind(dummy, occ)~id, dt= DT2, 
                # constraints: if sum is greater than this ignore 
                constraintsUB = c(Inf, 100), 
                # constraints: if sum is smaller than this ignore
                constraintsLB= c(-Inf, -Inf), 
                # weight the traits for the selection
                b = c(1,0), 
                # population parameters
                nCrosses = 100, nProgeny = 10, 
                # coancestry parameters
                D=As,
                lambda=(60*pi)/180, nQtlStart = 80, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 15, verbose = FALSE)

Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
best = bestSol(res$pop)[,1]
sum(Q[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
cex <- rep(0.5,nrow(PCWheat))
names(cex) <- rownames(PCWheat)
cex[names(which(Q[best,]==1))]=2
plot(PCWheat[,1], PCWheat[,2], col = TreeWheat, cex=cex,
     pch = TreeWheat, xlab = "pc1", ylab = "pc2")

## -----------------------------------------------------------------------------
DT2$cov <- apply(D[tp,vp],1,mean)

## -----------------------------------------------------------------------------
res<-evolafit(cbind(cov, occ)~id, dt= DT2, 
                # constraints: if sum is greater than this ignore 
                constraintsUB = c(Inf, 100), 
                # constraints: if sum is smaller than this ignore
                constraintsLB= c(-Inf, -Inf), 
                # weight the traits for the selection
                b = c(1,0), 
                # population parameters
                nCrosses = 100, nProgeny = 10, 
                # coancestry parameters
                D=As,
                lambda=(60*pi)/180, nQtlStart = 80, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 15, verbose = FALSE)

Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
best = bestSol(res$pop)[,1]
sum(Q[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
cex <- rep(0.5,nrow(PCWheat))
names(cex) <- rownames(PCWheat)
cex[names(which(Q[best,]==1))]=2
plot(PCWheat[,1], PCWheat[,2], col = TreeWheat, cex=cex,
     pch = TreeWheat, xlab = "pc1", ylab = "pc2")

## -----------------------------------------------------------------------------
data(DT_technow)
DT <- DT_technow
DT$occ <- 1; DT$occ[1]=0
M <- M_technow
D <- A.mat(M)

Z=with(DT,overlay(dent,flint) )#  Matrix::sparse.model.matrix(~dent-1, data=DT)
rownames(Z) <- DT$hy # needed to link to the QTL matrix

## -----------------------------------------------------------------------------
# regular fitness function
fitnessf <-function (Y, b, d, Q, Z) {
  fit <- Y %*% b - d
  return(fit)
}
# new fitness function with constraint
fitnessf <-function (Y, b, Q, D, a, lambda, scale=TRUE, Z) {
  X=Q%*%Z[colnames(Q),]
  bad <- as.vector( apply(X,1, function(x){length(which(x > 5))}) )
  bad <- which(bad > 0)
  # (q'a)b - l(q'Dq)
  if(scale){
    fit <- stan( apply(Y,2,scale) %*% b) -  lambda*stan( Matrix::diag(Q %*% Matrix::tcrossprod(D, Q)) )
  }else{
    fit <- stan( Y %*% b) -  lambda*stan( Matrix::diag(Q %*% Matrix::tcrossprod(D, Q)) )
  }
  if(length(bad) > 0){fit[bad,1]=min(fit[,1])}
  return(fit)
}


## -----------------------------------------------------------------------------
# silent for CRAN checks restriction on time
# res<-evolafit(formula = c(GY, occ)~hy,
#               dt= DT, 
#               # constraints: if sum is greater than this ignore
#               constraintsUB = c(Inf,50), 
#               # constraints: if sum is smaller than this ignore
#               constraintsLB= c(-Inf,-Inf),
#               # weight the traits for the selection
#               b = c(1,0), 
#               # population parameters
#               nCrosses = 100, nProgeny = 10, 
#               # coancestry parameters
#               D=D, lambda= (10*pi)/180 , nQtlStart = 40, 
#               # new fitness function and additional args
#               fitnessf = fitnessf, Z=Z,
#               # selection parameters
#               propSelBetween = 0.5, propSelWithin =0.5, 
#               nGenerations = 15, verbose=FALSE) 
# 
# Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
# best = bestSol(res$pop)[,1]
# qa = (Q %*% DT$GY)[best,]; qa 
# qDq = Q[best,] %*% D %*% Q[best,]; qDq 
# sum(Q[best,]) # total # of inds selected

## -----------------------------------------------------------------------------
# # check how many times an individual was used in the final crosses
# crosses <- data.frame(cross=names(which( Q[best,] == 1)))
# table(unlist(strsplit(crosses$cross,":")))
# # check performance of crosses selected
# plot(DT$GY, col=as.factor(Q[best,]), 
#        pch=(Q[best,]*19)+1)

## -----------------------------------------------------------------------------

data("mtcars")
# we scale the variables
mtcars <- as.data.frame(apply(mtcars,2,scale))
mtcars$inter <- 1 # add an intercept if desired

# define the train and validation set
train <- sample(1:nrow(mtcars) , round((nrow(mtcars)*.4)))
validate <- setdiff(1:nrow(mtcars),train)
mtcarsT <- mtcars[train,]
mtcarsV <- mtcars[validate,]

##############################
# fit the regular linear model
head(mtcarsT)
mod <- lm(mpg~cyl+disp+hp+drat, data=mtcarsT);mod

##############################
# fit the genetic algorithm
# 1) create initial QTL effects to evolve
nqtls=100
dt <- data.frame(alpha=rnorm(nqtls,0,.3),qtl=paste0("Q",1:nqtls))
head(dt); nrow(dt)

# generate n samples equivalent to the number of progeny
# you are planning to start the simulation with (e.g., 500)
# these are fixed values
sam <- sample(1:nrow(mtcarsT),500,replace = TRUE)
y <- mtcarsT$mpg[sam]
X = mtcarsT[sam,c("cyl","disp","hp","drat")]

# Task: linear regression
res0<-evolafit(alpha~qtl, dt= dt,
               # constraints: if greater than this ignore
               constraintsUB = c(Inf), 
               # constraints: if smaller than this ignore
               constraintsLB= c(-Inf), 
               # weight the traits for the selection
               b = c(1), 
               # population parameters
               nCrosses = 50, nProgeny = 10, recombGens = 1, 
               # coancestry parameters
               D=NULL, lambda=0, nQtlStart = 4, fixNumQtlPerInd = TRUE,
               # least MSE function (y - Xb)^2; Y are betas; X*Y is X*beta; 
               # Y and X are fixed, we just evolve the betas
               fitnessf=regFun,
               # selection parameters
               propSelBetween = 0.65, propSelWithin =0.65, selectTop=FALSE,
               nGenerations = 10, y=y, X=X, verbose = FALSE
) 

# check how the fitness function changed across generations
pmonitor(res0, kind = 1)
# this time the best solution is the one that minimizes the error
Q <- pullQtlGeno(res0$pop, simParam = res0$simParam, trait=1); Q <- Q/2
bestid <- bestSol(res0$pop, selectTop = FALSE)[,"fitness"]
bestid
betas <- res0$simParam$traits[[1]]@addEff[which(Q[bestid,] > 0)]
betas

# plot predicted versus real values
plot( as.matrix(mtcarsV[,c("cyl","disp","hp","drat")]) %*% betas  , mtcarsV$mpg,
      xlab="predicted mpg value by GA", ylab="mpg",
      main="Correlation between GA-prediction and observed") # GA
plot( as.matrix(mtcarsV[,c("inter","cyl","disp","hp","drat")]) %*% mod$coefficients , mtcarsV$mpg,
      xlab="predicted mpg value by lm", ylab="mpg",
      main="Correlation between lm-prediction and observed") # LM
# Correlation between GA-prediction and observed 
cor( as.matrix(mtcarsV[,c("cyl","disp","hp","drat")]) %*% betas  , mtcarsV$mpg) 
# Correlation between lm-prediction and observed
cor( as.matrix(mtcarsV[,c("inter","cyl","disp","hp","drat")]) %*% mod$coefficients , mtcarsV$mpg) # LM



## -----------------------------------------------------------------------------
# "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# function to simulate cities
simCities <- function(n_cities = 5) {
  # extend "LETTERS" function to run from "A" to "ZZ"
  MORELETTERS <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  cities <- matrix(runif(2*n_cities,-1,1),ncol=2)
  rownames(cities) <- MORELETTERS[1:n_cities]
  colnames(cities) <- c("x","y")
  return(cities)
}

# function to plot cities
plotCities <- function(cities, route, dist=NULL, 
                      main="", bg="white", 
                      main_col = "black",
                      point_col = "deepskyblue",
                      start_col = "red",
                      line_col = "black") {
  # plot cities
  city_colors <- c(start_col, rep(point_col,nrow(cities)))
  
  par(bg=bg)
  plot(cities, 
       pch=16, cex=3,
       col=point_col, 
       ylim=c(-1,1.1), xlim=c(-1,1),
       # yaxt="n", xaxt="n",
       # ylab="", xlab="",
       bty="n",
       main=main, col.main=main_col)
  text(x=cities[,1], y=cities[,2], labels=rownames(cities))
  # plot route
  if(!missing(route)){
    for(i in 1:length(route)) {
      if(route[i]>0){
        nodes0 <- strsplit(names(route)[i], "-")[[1]]
        lines(x=cities[nodes0,"x"],
              y=cities[nodes0,"y"],
              col=line_col,
              lwd=1.5)
      }
    }
  }
  # add distance in legend
  if(!is.null(dist)) legend("topleft", 
                            bty="n", # no box around legend
                            legend=round(dist,4), 
                            bg="transparent", 
                            text.col="black")
  
}

# function to compute adjacency matrix
compAdjMat <- function(cities) return(as.matrix(dist(cities)))

## -----------------------------------------------------------------------------
nCities=5
cities <- simCities(nCities)
adjmat <- compAdjMat(cities)
# make the distance matrix a data frame
df2 <- as.data.frame(as.table(adjmat)) 
df2 <- df2[-which(df2$Var1 == df2$Var2),]
colnames(df2)[3] <- c("distances")
df2$route <- paste(df2$Var1, df2$Var2, sep="-")
df2


## -----------------------------------------------------------------------------
H <- with(df2, evola::overlay(Var1,Var2) )
rownames(H) <- df2$route
head(H)

## -----------------------------------------------------------------------------
salesf <- function(Y,b,Q,D,a,lambda ,H, nCities){
                # simple fitness function
                fitnessVal <-  (Y%*%b) 
                ###############
                # CONSTRAINTS
                ###############
                # calculate how many cities the solution has travelled to
                QH <- Q %*% H
                # ensure all cities have at least 2 edges
                edgeCheck <- apply(QH,1,function(x){if(all(x>1)){return(1)}else{return(0)} })
                #apply condition on arriving and leaving the city
                fitnessVal <- ifelse(edgeCheck == 0, Inf, fitnessVal) # if not touching at least 2 cities give inf distance
                # number of unique cities visited
                nCityCheck <- apply(QH,1,function(x){length(which(x > 0))}) 
                # apply condition on visiting all cities
                fitnessVal <- ifelse(nCityCheck < nCities, Inf, fitnessVal) # if not in all cities give an infinite distance
                return(fitnessVal)
              }

## -----------------------------------------------------------------------------
res<-evolafit(formula=distances~route, dt= df2,
              # constraints on traits: if greater than this ignore
              constraintsUB = c(Inf), 
              # constraints on traits: if smaller than this ignore
              constraintsLB= c(-Inf), 
              # weight the traits for the selection (fitness function)
              b = c(1), 
              # population parameters
              nCrosses = 50, nProgeny = 10, 
              # genome parameters
              recombGens = 1, nChr=1, mutRateAllele=0, 
              # start with at least n QTLs equivalent to n cities
              nQtlStart = nCities*2, 
              # coancestry parameters
              D=NULL, lambda=0, 
              fitnessf = salesf, 
              selectTop = FALSE, 
              # additional variables for the fitness function
              H=H, nCities=nCities,
              # selection parameters
              # propSelBetween = .8, propSelWithin =0.8, 
              nGenerations = 50, verbose=FALSE
) 

pmonitor(res, kind=1) # fitness should decrease
Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
best <- bestSol(res$pop, selectTop = FALSE)[,"fitness"]
Q[best,] # routes taken
Q[best,] %*% H # cities visited (should have a 2 so we arrived and left once)

plotCities(cities, route=Q[best,])


