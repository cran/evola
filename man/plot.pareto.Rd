\name{pareto}
\alias{pareto}
\title{Pareto across iterations}
\description{\code{plot} for Pareto front.
}
\usage{
pareto(object, scaled=TRUE,trait=~fitness, pch=20, xlim, \dots)
}
\arguments{
 \item{object}{model object returned by \code{"evolafit"}}
 \item{scaled}{a logical value to specify the scale of the y-axis (gain in merit).}
 \item{trait}{a formula specifying the trait to be plotted. By default is the fitness.}
 \item{pch}{symbol for plotting points as desribed in par}
 \item{xlim}{upper and lower bound in the x-axis}
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
pareto(res00)
}

}
\seealso{\code{\link{plot}}, \code{\link{evolafit}}}
\keyword{models}

