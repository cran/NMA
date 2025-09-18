\docType{data}
\name{smoking}
\alias{smoking}
\title{Smoking cessation data}
\format{An arm-based dataset with 24 studies.}
\usage{
data(smoking)
}
\description{
A network meta-analysis data for smoking cessation from Lu and Ades (2006) and Higgins et al. (2012).
\itemize{
  \item \code{study}: Study ID.
  \item \code{trt}: A character variable that indicates the type of intervention, \code{A}: No contact, \code{B}: Self help, \code{C}: Individual counselling, \code{D}: Group counselling.
  \item \code{n}: Number of participants of the intervention.
  \item \code{d}: Number of successes of the intervention.
}
}
\references{
Lu, G., Ades, A. E. (2006).
Assessing evidence inconsistency in mixed treatment comparisons.
\emph{Journal of the American Statistical Association} \strong{101}:447-459.

Higgins, J. P. T., Jackson, D., Barrett, J. K. et al (2012).
Consistency and inconsistency in network meta-analysis: concepts and models for multi-arm studies.
\emph{Research Synthesis Methods} \strong{3}:98-110.
}
\keyword{datasets}
