\name{nmarank}
\alias{nmarank}
\title{Calculating ranking statistics of network meta-analysis}
\usage{
nmarank(x, B=20000, method="NH", ascending=TRUE)
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{B}{Number of parametric bootstrap resampling (default: 20000)}

\item{method}{Estimation and prediction method. \code{NH}: Noma-Hamura's improved REML-based methods (default). \code{REML}: The ordinary REML method. \code{fixed}: Fixed-effect model.}

\item{ascending}{A logical value that specifies whether the ranking is defined by ascending or descending order. Set \code{ascending=TRUE} if a smaller value of the effect measure indicates a better outcome (e.g., mortality, frequency of adverse events, or other "harmful" outcomes). In this case, the treatment with the smallest estimate will be ranked first. Set \code{ascending=FALSE} if a larger value of the effect measure indicates a better outcome (e.g., cure rate, response rate, or other "beneficial" outcomes). In this case, the treatment with the largest estimate will be ranked first.}
}

\value{
The results of the ranking statistics of network meta-analysis are provided. Also, ranking probability plots are produced.
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
Chaimani, A., Higgins, J. P., Mavridis, D., Spyridonos, P., and Salanti, G. (2013). Graphical tools for network meta-analysis in STATA. \emph{PLoS One} \strong{8}, e76654.

Salanti, G., Ades, A. E. and Ioannidis, J. P. (2011). Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: An overview and tutorial. \emph{Journal of Clinical Epidemiology} \strong{64}, 163–171.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

nmarank(hf2)
nmarank(hf2, ascending=FALSE)
}
