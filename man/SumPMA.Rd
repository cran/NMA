\name{SumPMA}
\alias{SumPMA}
\title{Pairwise meta-analyses for all treatment pairs with direct comparisons on the network}
\usage{
SumPMA(x,method="REML",test="z")
}
\arguments{
\item{x}{Output object of \code{setup}}
\item{method}{Method of the estimation of pairwise meta-analysis. All possible options of \code{rma} function in \code{metafor} package is available (default: the REML method (\code{REML})).}
\item{test}{Method of the statistical inference for pairwise meta-analysis. All possible options of \code{rma} function in \code{metafor} package is available (default: the standard normal approximation (\code{z})).}
}
\value{
The results of the meta-analyses for all possible treatment pairs are provided.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{measure}: Outcome measure.
\item \code{method}: Estimation method.
\item \code{test}: Inference method.
\item \code{Summary effect measures}: N (number of studies), summary estimates, 95\% confidence intervals, and P-values for all possible pairs.
\item \code{Heterogeneity measures}: N (number of studies), tau2 (heterogeneity variance) estimate, I2-statistic, and H2-statistic.
\item \code{Egger test}: N (number of studies), P-value of the Egger test for assessing publication bias.
}
}
\description{
Pairwise meta-analyses for all treatment pairs with direct comparisons on the network are performed. The synthesis analyses are performed by \code{rma} and \code{regtest} in \code{metafor} package.
}
\references{
DerSimonian, R., and Laird, N. M. (1986). Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} \strong{7}, 177-188.

Egger, M., Davey Smith, G., Schneider, M., and Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. \emph{BMJ} \strong{315}, 629-634.

Higgins, J. P. T., and Thompson, S. G. (2002). Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine} \strong{21}, 1539-1558.

IntHout, J., Ioannidis, J. P. A., and Borm, G. F. (2014). The Hartung–Knapp–Sidik–Jonkman method for random effects meta‐analysis is straightforward and considerably outperforms the standard DerSimonian–Laird method. \emph{BMC Medical Research Methodology} \strong{14}, 25.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software} \strong{36}, Issue 3.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

SumPMA(hf2,method="REML",test="z")
}
