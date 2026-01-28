\name{nmaweight}
\alias{nmaweight}
\title{Evaluating study weights and contribution matrix}
\usage{
nmaweight(x, digits=3)
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{digits}{Number of decimal places}
}
\value{
Contribution weight matrices for the consistency model are provided. Also, a heatmap for the contribution matrix of overall evidence is presented.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{Contribution of direct and indirect information}: Contribution rates of direct and indirect evidence for individual treatment pairs.
\item \code{Contribution weights: Direct comparison}: Contribution weight matrix for direct evidence.
\item \code{Contribution weights: Indirect comparison (BoS)}: Contribution weight matrix for indirect evidence (BoS; borrowing of strength of Jackson et al. (2017)).
\item \code{Contribution weights: Overall evidence}: Contribution weight matrix for overall evidence.
}
}
\description{
Contribution weight matrices to assess how individual studies influence the synthesized results are presented. Jackson et al. (2017) and Noma et al. (2017) showed the contribution rates are estimated by the factorized information, and the contribution weight matrices are calculated through the factorized information.
}
\references{
Jackson, D., White, I. R., Price, M., Copas, J., and Riley, R. D. (2017). Borrowing of strength and study weights in multivariate and network meta-analysis. \emph{Statistical Methods in Medical Research} \strong{26}, 2853-2868.

Noma, H., Tanaka, S., Matsui, S., Cipriani, A., and Furukawa, T. A. (2017). Quantifying indirect evidence in network meta-analysis. \emph{Statistics in Medicine} \strong{36}, 917-927.
}
\examples{
data(smoking)

smk2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="A",data=smoking)

nmaweight(smk2)
}
