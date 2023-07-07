\name{nmarank}
\alias{nmarank}
\title{Calculating the ranking statistics of network meta-analysis}
\usage{
nmarank(x, B=20000, method="NH", ascending=TRUE)
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{B}{Number of parametric bootstrap resampling (default: 20000)}

\item{method}{Estimation and prediction method. \code{NH}: Noma-Hamura's improved REML-based methods (default). \code{REML}: The ordinary REML method. \code{fixed}: Fixed-effect model.}

\item{ascending}{A logical value that specify whether the ranking is defined by ascending or descending order}
}
\value{
Results of the ranking statistics of network meta-analysis are provided. Also, ranking probability plots are generated.
\itemize{
\item \code{SUCRA}: SUCRA estimates of individual treatment by parametric bootstrap.
\item \code{MEANRANK}: Mean rank estimates of individual treatment by parametric bootstrap.
\item \code{Probability of ranking}: Probability of ranking (best, 2nd, 3rd,..., worst) estimates of individual treatment by parametric bootstrap.
}
}
\description{
Ranking statistics of network meta-analysis such as SUCRA, MEANRANK, and probability of ranking are calculated by parametric bootstrap.
}
\references{
Salanti, G., Ades, A. E. and Ioannidis, J. P. (2011). Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: An overview and tutorial. \emph{Journal of Clinical Epidemiology} \strong{64}, 163–171.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

nmarank(hf2)
nmarank(hf2, ascending=FALSE)
}
