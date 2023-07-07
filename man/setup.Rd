\name{setup}
\alias{setup}
\title{Transforming arm-level data to contrast-based summary statistics and making objects for the network meta-analysis}
\usage{
setup(study,trt,d,n,m,s,measure,ref,data)
}
\arguments{
\item{study}{Study ID}
\item{trt}{Treatment variable. It can be formed as both of numbered treatment (=1,2,...) and characters (e.g., "Placebo", "ARB", "Beta blocker").}
\item{d}{Number of events (for dichotomous outcome).}
\item{n}{Sample size (for dichotomous and continuous outcome).}
\item{m}{Mean of the outcome variable (for continuous outcome).}
\item{s}{Standard deviation of the outcome variable (for continuous outcome).}
\item{measure}{Outcome measure (can be \code{OR} (odds ratio), \code{RR} (risk ratio), and \code{RD} (risk difference) for dichotomous outcome, and \code{MD} (mean difference) and \code{SMD} (standardized mean difference) for continuous outcome.}
\item{ref}{Reference treatment category that should be involved in \code{trt}.}
\item{data}{A data frame that involves the arm-based data.}
}
\value{
Contrast-based summary statistics are generated.
\itemize{
\item \code{coding}: A table that presents the correspondence between the numerical code and treatment categories (the reference category is coded as 1).
\item \code{reference}: Reference treatment category.
\item \code{measure}: Outcome measure.
\item \code{N}: The number of study.
\item \code{p}: The dimension of the contrast-based statistics.
\item \code{df}: The degree of freedom.
\item \code{study}: The ID variable that specifies studies.
\item \code{trt}: The original vector that specifies treatment categories.
\item \code{treat}: A numerical vector that specifies treatment categories based on the \code{coding} table.
\item \code{d}: The original vector that specifies number of events.
\item \code{n}: The original vector that specifies sample sizes.
\item \code{m}: The original vector that specifies means.
\item \code{s}: The original vector that specifies standard deviations.
\item \code{y}: Contrast-based summary estimates.
\item \code{S}: Vectored within-study covariance matrix.
}
}
\description{
A setup function to generate R objects that can be used for various network meta-analysis functions in this package. Users should prepare arm-level datasets from literatures, and this function can transform them to the contrast-based summary statistics.
Both of dichotomous and continuous outcomes can be handled. The type of outcome variable can be specified by the \code{measure}. If the \code{measure} is specified as \code{OR}, \code{RR} or \code{RD}, the outcome should be dichotomous, and \code{d} and \code{n} are needed to compute the summary statistics.
Besides, if the \code{measure} is specified as \code{MD} or \code{SMD}, the outcome should be continuous, and \code{m}, \code{s} and \code{n} are needed to compute the summary statistics.
The output objects of \code{setup} function is used for the other computational or graphical tools in this package.
}
\references{
Noma, H. (2023). Within-study covariance estimators for network meta-analysis with contrast-based approach. Forthcoming.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)
hf3 <- setup(study=study,trt=trt,d=d,n=n,measure="RR",ref="Placebo",data=heartfailure)
hf4 <- setup(study=study,trt=trt,d=d,n=n,measure="RD",ref="Placebo",data=heartfailure)

data(antidiabetic)

ad2 <- setup(study=id,trt=t,m=y,s=sd,n=n,measure="MD",ref="Placebo",data=antidiabetic)
ad3 <- setup(study=id,trt=t,m=y,s=sd,n=n,measure="SMD",ref="Placebo",data=antidiabetic)
}
