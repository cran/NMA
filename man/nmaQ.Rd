\name{nmaQ}
\alias{nmaQ}
\title{Multivariate Q-statistic and its factorization}
\usage{
nmaQ(x)
}
\arguments{
\item{x}{Output object of \code{setup}}
}
\value{
Multivariate Q-statistic and its factorized ones (within and between designs) are provided.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{number of studies}: The number of synthesized studies.
\item \code{Within designs (individual designs)}: Q-statistics for individual designs and their P-values.
\item \code{Q-statistics}: Multivariate Q-statistics and its factorized ones (within and between designs), and their P-values.
}
}
\description{
Multivariate Q-statistic and its factorized versions (within and between designs) are provided. P-values of the corresponding Q-tests are also presented.
}
\references{
Jackson, D., White, I. R., and Riley, R. D. (2012). Quantifying the impact of between-study heterogeneity in multivariate meta-analyses. \emph{Statistics in Medicine} \strong{31}: 3805-3820. 

Krahn, U., Binder, H., and Konig, J. (2013). A graphical tool for locating inconsistency in network meta-analysis.  \emph{BMC Medical Research Methodology} \strong{13}, 35.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

nmaQ(hf2)
}
