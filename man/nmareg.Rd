\name{nmareg}
\alias{nmareg}
\title{Network meta-regression based on contrast-based approach}
\usage{
nmareg(x, z, treats)
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{z}{Covariate name vector}

\item{treats}{A vector that specifies treatments to be assessed effect modifications that correspond to the elements of outcome vectors \code{y} in \code{x} (please specify the treatment numbers of \code{coding}; multiple outcomes can be specified jointly, as a vector).}
}
\value{
The results of the network meta-regression analysis are provided.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{Covariates}: Covariate that specified in \code{setup}.
\item \code{Outcome evaluated the effect modifications}: Treatment contrasts that the effect modifications are evaluated.
\item \code{Coefficients}: Estimates, their SEs, Wald-type 95\% confidence intervals, and P-values for the regression parameters (\code{cons}: intercept, \code{beta}: regression coefficient for the explanatory variable).
\item \code{Between-studies_SD}: Between-studies SD (tau) estimate.
\item \code{Between-studies_COR}: Between-studies correlation coefficient (ought to be 0.50).
}
}
\description{
Network meta-regression based on contrast-based approach using the multivariate meta-regression model. Effect modifications by study-level covariates (specified in the \code{setup} function) can be assessed. In many network meta-analysis, some treatment contrasts involve only 1 or 2 (or 0) direct comparisons, and the regression coefficients of the corresponding outcomes cannot be validly estimated (non-identifiable). Thus, the \code{nmareg} function can specify a subset of outcome variables to be modelled by the regression model (to be assessed the effect modifications) by \code{treats}. Currently, the parameter estimation is performed by the ordinary REML method.
}
\references{
Nikolakopoulou, A., White, I. R., Salanti, G. (2021). Network meta-analysis. In: Schmid, C. H., Stijnen, T., White, I. R., eds. \emph{Handbook of Meta-Analysis}. CRC Press; pp. 187-217.

Noma, H., Hamura, Y., Gosho, M., and Furukawa, T. A. (2023). Kenward-Roger-type corrections for inference methods of network meta-analysis and meta-regression. \emph{Research Synthesis Methods} \strong{14}, 731-741.

White, I. R., Barrett, J. K., Jackson, D., and Higgins, J. P. (2012). Consistency and inconsistency in network meta-analysis: model estimation using multivariate meta-regression. \emph{Research Synthesis Methods} \strong{3}, 111-125.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="OR",
ref="Placebo",data=heartfailure)

nmareg(hf2,z=SBP,treats=3)
nmareg(hf2,z=c(SBP,DBP),treats=c(3,4,6))
}
