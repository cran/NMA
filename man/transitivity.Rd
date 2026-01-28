\name{transitivity}
\alias{transitivity}
\title{Checking transitivity}
\usage{
transitivity(x, z, gcol="blue", yrange, digits=3)
}
\arguments{
\item{x}{Output object of \code{setup}}

\item{z}{Covariate name for assessing transitivity (must be involved in \code{covariate} of the output object of \code{setup}}

\item{gcol}{Color of the plot}

\item{yrange}{Range of y-axis of the plot}

\item{digits}{Number of decimal places}
}
\value{
Summary statistics of the covariate among different study designs and its summary plot are presented.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{covariate}: Covariate that specified in \code{setup}.
\item \code{summary}: Summary of the covariate among different study designs. \code{N}: number of the corresponding studies, \code{n}: total sample size, \code{wt.mean}: weighted mean, \code{min}: minimum, \code{max}: maximum.
}
}
\description{
To check transitivity on the network, summary statistics of a certain covariate among different study designs are provided. Also, a summary plot for these statistics is presented.
}
\references{
Cipriani, A., Furukawa, T. A., Salanti, G., et al. (2018). Comparative efficacy and acceptability of 21 antidepressant drugs for the acute treatment of adults with major depressive disorder: a systematic review and network meta-analysis. \emph{Lancet} \strong{391}, 1357-1366.

Salanti, G. (2012). Indirect and mixed-treatment comparison, network, or multiple-treatments meta-analysis: many names, many benefits, many concerns for the next generation evidence synthesis tool. \emph{Research Synthesis Methods} \strong{3}, 80-97. 
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="OR",
ref="Placebo",data=heartfailure)

transitivity(hf2, SBP)
transitivity(hf2, DBP)
transitivity(hf2, pubyear)
}
