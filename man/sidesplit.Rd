\name{sidesplit}
\alias{sidesplit}
\title{Sidesplitting for quantifying direct and indirect evidence for all possible treatment pairs and the inconsistency test}
\usage{
sidesplit(x)
}
\arguments{
\item{x}{Output object of \code{setup}}
}
\value{
The results of the sidesplitting for all possible treatment pairs are provided.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{reference}: Reference treatment category.
\item \code{Direct evidence}: Summary estimates, SEs, 95\% confidence intervals, and P-values for the direct evidence.
\item \code{Indirect evidence}: Summary estimates, SEs, 95\% confidence intervals, and P-values for the indirect evidence.
\item \code{Difference}: Differences of the summary estimates of direct and indirect evidence, and their inconsistency tests.
}
}
\description{
Noma's sidesplitting for quantifying direct and indirect evidence for all possible treatment pairs based on network meta-regression and the inconsistency tests are performed. For the bias correction that causes the involvement of multi-arm trials, we adopted the adjustment method of Noma et al. (2017) and Noma (2023).
}
\references{
Dias, S., Welton, N. J., Caldwell, D. M., and Ades, A. E. (2010). Checking consistency in mixed treatment comparison meta-analysis. \emph{Statistics in Medicine} \strong{29}, 932-944.

Noma, H. (2024). Sidesplitting using network meta-regression. \emph{Japanese Journal of Biometrics} \strong{44}, 107-118.

Noma, H., Tanaka, S., Matsui, S., Cipriani, A., and Furukawa, T. A. (2017). Quantifying indirect evidence in network meta-analysis. \emph{Statistics in Medicine} \strong{36}, 917-927.
}
\examples{
data(smoking)

smk2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="A",data=smoking)

sidesplit(smk2)
}
