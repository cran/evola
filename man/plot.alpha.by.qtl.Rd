\name{abq}
\alias{abq}
\title{Alpha by QTL plot}
\description{

A \code{plot} to observe the QTLs activated for the best solution(s) given a certain 
average effect of an allele substitution (alpha).

}
\usage{
abq(object, traitAlpha, selectTop=TRUE, n=1, \dots)
}
\arguments{
 \item{object}{model object returned by \code{evolafit}}
 \item{traitAlpha}{a character value to specify the trait alpha to be ploted in the y axis.}
 \item{selectTop}{a logical value to be passed internally to the \code{bestSol} function
 representing whether the best solutions selected should be the ones with the highest (TRUE)
 or lowest (FALSE) fitness value.}
 \item{n}{an integer value to be passed internally to the \code{bestSol} function 
 representing the number of best solutions.}
 \item{\dots}{Further arguments to be passed to the plot function.}
}
\value{
  vector of plot
 }
\author{Giovanny Covarrubias}
\examples{

set.seed(1)

# Data
Gems <- data.frame(
  Color = c("Red", "Blue", "Purple", "Orange",
            "Green", "Pink", "White", "Black", 
            "Yellow"),
  Weight = round(runif(9,0.5,5),2),
  Value = round(abs(rnorm(9,0,5))+0.5,2),
  Times=c(rep(1,8),0)
)

 \donttest{
# simple specification
res00<-evolafit(formula=cbind(Weight,Value)~Color, dt= Gems,
                # constraints on traits: if greater than this ignore
                constraintsUB = c(10,Inf), nGenerations = 10
)
abq(res00, traitAlpha = "Value")
}

}
\seealso{\code{\link{plot}}, \code{\link{evolafit}}}
\keyword{models}

