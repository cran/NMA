\name{random.icm}
\alias{random.icm}
\title{Jackson's random inconsistency model}
\usage{
random.icm(x, digits=3)
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{digits}{Number of decimal places}
}
\value{
The results of the analysis of Jackson's random inconsistency model and I2-statistics are provided.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{reference}: Reference treatment category.
\item \code{number of studies}: Number of studies.
\item \code{number of designs}: Number of designs.
\item \code{designs}: Study designs (combinations of treatments of individual trials) on the network.
\item \code{Coef. (vs. treat 1)}: Regression coefficients estimates and their SEs, 95\% confidence intervals and P-values.
\item \code{Between-studies_SD}: Between-studies SD estimate.
\item \code{Between-designs_SD}: Between-designs SD estimate.
\item \code{Likelihood ratio tests for the variance components}: The results of the likelihood ratio tests for comparing (1) the fixed- and random-effects models without inconsistency effects (heterogeneity), (2) the random-effects models with and without inconsistency effects (inconsistency), and (3) the fixed-effect model without inconsistency effects and the random-effects model with inconsistency effects (heterogeneity + inconsistency).
\item \code{Heterogeneity and inconsistency statistics}: R-statistics and I2-statistics for comparing  (1) the fixed- and random-effects models without inconsistency effects (heterogeneity), (2) the random-effects models with and without inconsistency effects (inconsistency), and (3) the fixed-effect model without inconsistency effects and the random-effects model with inconsistency effects (heterogeneity + inconsistency).
}
}
\description{
Jackson's random inconsistency model for modelling the design-by-treatment interactions. Model-based testing results for heterogeneity and inconsistency (design-by-treatment interactions) and the I2-statistics are provided.
}
\references{
Jackson, D., White, I. R., and Riley, R. D. (2012). Quantifying the impact of between-study heterogeneity in multivariate meta-analyses. \emph{Statistics in Medicine} \strong{31}: 3805-3820. 

Jackson, D., Barrett, J. K., Rice, S., White, I. R., and Higgins, J. P. T. (2014). A design-by-treatment interaction model for network meta-analysis with random inconsistency effects. \emph{Statistics in Medicine} \strong{33}, 3639-3654.

Law, M., Jackson, D., Turner, R., Rhodes, K., and Viechtbauer, W. (2016). Two new methods to fit models for network meta-analysis with random inconsistency effects. \emph{BMC Medical Research Methodology} \strong{16}, 87.

Nikolakopoulou, A., White, I. R., and Salanti, G. (2021). Network meta-analysis. In: Schmid, C. H., Stijnen, T., White, I. R., eds. \emph{Handbook of Meta-Analysis}. CRC Press; pp. 187-217.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

random.icm(hf2)
}
