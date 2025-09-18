\docType{data}
\name{heartfailure}
\alias{heartfailure}
\title{Sciarretta et al. (2011)'s network meta-analysis data}
\format{An arm-based dataset with 26 studies}
\usage{
data(heartfailure)
}
\description{
A network meta-analysis data from Sciarretta et al. (2011) that compared 7 antihypertensive drug classes and placebo for occurrence of heart failure.
\itemize{
  \item \code{study}: Study ID
  \item \code{trial}: Trial name
  \item \code{trt}: Treatment (AB (Alpha blocker), ACE (ACE inhibitor), ARB, BB (Beta blocker), CCB (Calcium-channel blocker), CT (conventional treatments), Diuretic (DD), Placebo)
  \item \code{n}: Sample size
  \item \code{d}: Number of events (occurrence of heart failure)
  \item \code{pubyear}: Publication year
  \item \code{SBP}: Mean of baseline systolic blood pressure (mmHg)
  \item \code{DBP}: Mean of diastolic systolic blood pressure (mmHg)
}
}
\references{
Sciarretta, S., Palano, F., Tocci, G., Baldini, R., and Volpe, M. (2011). 
Antihypertensive treatment and development of heart failure in hypertension: a Bayesian network meta-analysis of studies in patients with hypertension and high cardiovascular risk. 
\emph{Archives of Internal Medicine} 
\strong{171}: 384-394.
}
\keyword{datasets}
