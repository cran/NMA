\name{nmafunnel}
\alias{nmafunnel}
\title{Comparison-adjusted funnel plot}
\usage{
nmafunnel(x, method="NH", legends="topright")
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{method}{Estimation and prediction method. \code{NH}: Noma-Hamura's improved REML-based methods (default). \code{REML}: The ordinary REML method.}

\item{legends}{Location of the legend on the plot (default: \code{topright})}
}
\value{
Comparison-adjusted funnel plot for the studies involving treatment 1 (reference treatment specified in \code{setup}) is produced.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{summary}: \code{design}: design of studies, \code{N}: number of the corresponding studies, \code{n}: total sample size.
}
}
\description{
A comparison-adjusted funnel plot for the studies involving treatment 1 (reference treatment specified in \code{setup}) is produced.
}
\references{
Chaimani, A. and Salanti, G. (2012). Using network meta-analysis to evaluate the existence of small-study effects in a network of interventions. \emph{Research Synthesis Methods} \strong{3}, 161–176.

Chaimani, A., Higgins, J. P., Mavridis, D., Spyridonos, P., and Salanti, G. (2013). Graphical tools for network meta-analysis in Stata. \emph{PLoS One} \strong{8}, e76654.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)
hf4 <- setup(study=study,trt=trt,d=d,n=n,measure="RD",ref="Placebo",data=heartfailure)

nmafunnel(hf2,legends="bottomright")
nmafunnel(hf4)
}
