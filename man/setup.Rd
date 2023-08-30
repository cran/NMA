\name{setup}
\alias{setup}
\title{Transforming arm-level data to contrast-based summary statistics and making objects for network meta-analysis}
\usage{
setup(study,trt,d,n,m,s,z,measure,ref,data)
}
\arguments{
\item{study}{Study ID}
\item{trt}{Treatment variable. It can be formed as both of numbered treatment (=1,2,...) and characters (e.g., "Placebo", "ARB", "Beta blocker").}
\item{d}{Number of events (for dichotomous outcome).}
\item{n}{Sample size (for dichotomous and continuous outcome).}
\item{m}{Mean of the outcome variable (for continuous outcome).}
\item{s}{Standard deviation of the outcome variable (for continuous outcome).}
\item{z}{Covariate name vector to be used for network meta-regression analysis or transitivity analysis (optional).}
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
\item \code{covariate}: Covariate name(s).
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
\item \code{Z}: The data frame that specifies covariates matrix (design matrix).
\item \code{y}: Contrast-based summary estimates.
\item \code{S}: Vectored within-study covariance matrix.
}
}
\description{
A setup function to generate R objects that can be used for network meta-analysis. 
Users should prepare arm-level datasets, and the \code{setup} function transforms the arm-level data to the contrast-based summary statistics. 
Both of dichotomous and continuous outcomes can be treated. The type of outcome variable can be specified by the \code{measure}. If the \code{measure} is specified as \code{OR}, \code{RR} or \code{RD}, the outcome should be dichotomous, and \code{d} and \code{n} are needed to compute the summary statistics.
Besides, if the \code{measure} is specified as \code{MD} or \code{SMD}, the outcome should be continuous, and \code{m}, \code{s} and \code{n} are needed to compute the summary statistics. 
Several covariates can be involved as \code{z} for network meta-regression analysis (\code{nmareg}) and transitivity analysis (\code{transitivity}). 
}
\references{
Noma, H. (2023). Within-study covariance estimators for network meta-analysis with contrast-based approach. Jxiv, 490. \doi{10.51094/jxiv.490}.
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)
hf3 <- setup(study=study,trt=trt,d=d,n=n,measure="RR",ref="Placebo",data=heartfailure)
hf4 <- setup(study=study,trt=trt,d=d,n=n,measure="RD",ref="Placebo",data=heartfailure)

hf5 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="OR",
ref="Placebo",data=heartfailure)

data(antidiabetic)

ad2 <- setup(study=id,trt=t,m=y,s=sd,n=n,measure="MD",ref="Placebo",data=antidiabetic)
ad3 <- setup(study=id,trt=t,m=y,s=sd,n=n,measure="SMD",ref="Placebo",data=antidiabetic)
}
