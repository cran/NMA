\name{obj.forest}
\alias{obj.forest}
\title{Numerical objects of ranked forest plot for the synthesis results of network meta-analysis}
\usage{
obj.forest(x,method="NH",digits=3,ascending=TRUE)
}
\arguments{
\item{x}{Output object of \code{setup}}
\item{method}{Estimation and prediction method. \code{NH}: Noma-Hamura's improved REML-based methods (default). \code{REML}: The ordinary REML method. \code{fixed}: Fixed-effect model.}
\item{digits}{Number of decimal places}
\item{ascending}{Type of order. Default is ascending order, but it can be changed to descending order changing to \code{FALSE}.}
}
\value{
Numerical objects of ranked forest plot is generated. They can be used for \code{forestplot} function of \code{forestplot} package to make a customized ranked forest plot.
\itemize{
\item \code{labeltext}: A matrix that presents the label text table of the forestplot.
\item \code{coef}: A matrix that presents the point estimates and confidence limits.
\item \code{boxsize}: A vector that indicates the boxsizes.
}
}
\description{
Numerical objects of ranked forest plot for the synthesis results of network meta-analysis are generated. These objects can be used to make a customized forest plot using \code{forestplot} function of \code{forestplot} package.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

obj.forest(hf2)
}
