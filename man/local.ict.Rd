\name{local.ict}
\alias{local.ict}
\title{Local inconsistency tests for all closed loops on the network}
\usage{
local.ict(x)
}
\arguments{
\item{x}{Output object of \code{setup}}
}
\value{
Results of the local inconsistency tests for all closed loops on the network are presented.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{reference}: Reference treatment category.
\item \code{N}: Number of studies.
\item \code{tau}: Between-studies SD estimate.
\item \code{X2-statistic}: Chi-squared statistics of the generalized Bucher's test.
\item \code{df}: Degree of freedom.
\item \code{P-value}: P-value of the generalized Bucher's test.
}
}
\description{
Local inconsistency tests for all closed loops on the network are performed. Higgins' inconsistency test (Generalized Bucher's test) that assesses the design-by-treatment interactions on the triangle loops are performed and their results are presented.
}
\references{
Bucher, H. C., Guyatt, G. H., Griffith, L. E., and Walter, S. D. (1997). The results of direct and indirect treatment comparisons in meta-analysis of randomized controlled trials. \emph{Journal of Clinical Epidemiology} \strong{50}, 683-691.

Higgins, J. P., Jackson, D., Barrett, J. K., Lu, G., Ades, A. E., and White, I. R. (2012). Consistency and inconsistency in network meta-analysis: concepts and models for multi-arm studies. \emph{Research Synthesis Methods} \strong{3}, 98-110.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

local.ict(hf2)
}
