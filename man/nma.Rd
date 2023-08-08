\name{nma}
\alias{nma}
\title{Network meta-analysis based on contrast-based approach using the multivariate meta-analysis model}
\usage{
nma(x, eform=FALSE, method="NH")
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{eform}{A logical value that specify whether the outcome should be transformed by exponential function (default: FALSE)}

\item{method}{Estimation and prediction method. \code{NH}: Noma-Hamura's improved REML-based methods (default). \code{REML}: The ordinary REML method. \code{fixed}: Fixed-effect model.}
}
\value{
Results of the network meta-analysis using the multivariate meta-analysis model.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{reference}: Reference treatment category.
\item \code{number of studies}: The number of synthesized studies.
\item \code{method}: The estimation and prediction methods.
\item \code{Coef. (vs. treat1)}: Estimates, their SEs, Wald-type 95\% confidence intervals, and P-values for the grand mean parameter vector.
\item \code{tau (Between-studies_SD) estimate}: Between-studies SD (tau) estimate.
\item \code{tau2 (Between-studies_variance) estimate}: Between-studies variance (tau^2) estimate.
\item \code{Multivariate H2-statistic}: Jackson's multivariate H2-statistic.
\item \code{Multivariate I2-statistic}: Jackson's multivariate I2-statistic.
\item \code{Test for Heterogeneity}: Multivariate Q-statistic and P-value of the test for heterogeneity.
\item \code{95\%PI}: 95\% prediction intervals.
}
}
\description{
Network meta-analysis based on contrast-based approach using the multivariate random-effects meta-analysis model. The synthesis results and prediction intervals based on the consistency assumption are provided. The ordinary REML method and its improved higher order asymptotic methods (Noma-Hamura methods) are available.
}
\references{
Jackson, D., White, I. R., Riley, R. D. (2012). Quantifying the impact of between-study heterogeneity in multivariate meta-analyses. \emph{Statistics in Medicine} \strong{31}: 3805-3820. 

Nikolakopoulou, A., White, I. R., and Salanti, G. (2021). Network meta-analysis. In: Schmid, C. H., Stijnen, T., White, I. R., eds. \emph{Handbook of Meta-Analysis}. CRC Press; pp. 187-217.

Noma, H., Hamura, Y., Gosho, M., and Furukawa, T. A. (2023). Kenward-Roger-type corrections for inference methods of network meta-analysis and meta-regression. \emph{Research Synthesis Methods}, \doi{10.1002/jrsm.1652}.

Noma, H., Hamura, Y., Sugasawa, S., and Furukawa, T. A. (2023). Improved methods to construct prediction intervals for network meta-analysis. \emph{Research Synthesis Methods}, \doi{10.1002/jrsm.1651}.

White, I. R., Barrett, J. K., Jackson, D., and Higgins, J. P. (2012). Consistency and inconsistency in network meta-analysis: model estimation using multivariate meta-regression. \emph{Research Synthesis Methods} \strong{3}, 111-125.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)
hf3 <- setup(study=study,trt=trt,d=d,n=n,measure="RR",ref="Placebo",data=heartfailure)
hf4 <- setup(study=study,trt=trt,d=d,n=n,measure="RD",ref="Placebo",data=heartfailure)

nma(hf2, eform=TRUE)
nma(hf3, eform=TRUE)
nma(hf4)
}
