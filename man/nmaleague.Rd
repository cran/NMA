\name{nmaleague}
\alias{nmaleague}
\title{Generating a league table}
\usage{
nmaleague(x, method="NH", eform=FALSE, digits=3, PI=FALSE, out.csv=NULL)
}
\arguments{
\item{x}{Output object of \code{setup}}
\item{method}{Estimation and prediction method. \code{NH}: Noma-Hamura's improved REML-based methods (default). \code{REML}: The ordinary REML method.}
\item{eform}{A logical value that specify whether the outcome should be transformed by exponential function (default: FALSE)}
\item{digits}{Number of decimal places}
\item{PI}{A logical value that specify whether the inference or prediction results are provided}
\item{out.csv}{A character object that specify a filename if the user wants to output the league table as a CSV file (e.g., \code{out.csv="out_league.csv"}).}
}
\value{
A league table is generated.
}
\description{
A league table is generated for all possible pairs of the treatments. The league table can be outputted as a CSV file through setting \code{out.csv="filename"}.
}
\references{
Nikolakopoulou, A., White, I. R., and Salanti, G. (2021). Network meta-analysis. In: Schmid, C. H., Stijnen, T., White, I. R., eds. \emph{Handbook of Meta-Analysis}. CRC Press; pp. 187-217.

Noma, H., Hamura, Y., Gosho, M., and Furukawa, T. A. (2023). Kenward-Roger-type corrections for inference methods of network meta-analysis and meta-regression. \emph{Research Synthesis Methods} \strong{14}, 731-741.

Noma, H., Hamura, Y., Sugasawa, S., and Furukawa, T. A. (2023). Improved methods to construct prediction intervals for network meta-analysis. \emph{Research Synthesis Methods} \strong{14}, 794-806.

White, I. R., Barrett, J. K., Jackson, D., and Higgins, J. P. (2012). Consistency and inconsistency in network meta-analysis: model estimation using multivariate meta-regression. \emph{Research Synthesis Methods} \strong{3}, 111-125.
}
\examples{
data(smoking)

smk2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="A",data=smoking)

nmaleague(smk2)                                        # default setting
nmaleague(smk2, eform=TRUE)                            # transformed to exponential-scale
nmaleague(smk2, eform=TRUE, digits=2)                  # digits can be changed
nmaleague(smk2, eform=TRUE, PI=TRUE)                   # prediction intervals
}
