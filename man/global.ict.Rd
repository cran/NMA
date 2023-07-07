\name{global.ict}
\alias{global.ict}
\title{Higgins' global inconsistency test}
\usage{
global.ict(x)
}
\arguments{
\item{x}{Output object of \code{setup}}
}
\value{
Results of the global inconsistency test are presented.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{reference}: Reference treatment category.
\item \code{number of studies}: Number of studies.
\item \code{designs}: Study designs (combinations of treatments of individual trials) on the network.
\item \code{Coefficients of the design-by-treatment interaction model}: Regression coefficients estimates and their SEs, 95\% confidence intervals and P-values.
\item \code{Between-studies_SD}: Between-studies SD estimate.
\item \code{Between-studies_COR}: Between-studies correlation coefficient estimate (=0.50).
\item \code{X2-statistic}: Chi-squared statistic of the global inconsistency test.
\item \code{df}: Degree of freedom.
\item \code{P-value}: P-value of the global inconsistency test.
}
}
\description{
Higgins' global inconsistency test based on the design-by-treatment interaction model. REML-based Wald test for the all possible design-by-treatment interactions on the network is performed.
}
\references{
Higgins, J. P., Jackson, D., Barrett, J. K., Lu, G., Ades, A. E., and White, I. R. (2012). Consistency and inconsistency in network meta-analysis: concepts and models for multi-arm studies. \emph{Research Synthesis Methods} \strong{3}, 98-110.

Jackson, D., Boddington, P., and White, I. R. (2016). The design-by-treatment interaction model: a unifying framework for modelling loop inconsistency in network meta-analysis. \emph{Research Synthesis Methods} \strong{7}, 329-332.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

global.ict(hf2)
}
