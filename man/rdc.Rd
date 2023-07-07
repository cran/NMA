\name{rdc}
\alias{rdc}
\title{Rounding a numerical value}
\usage{
rdc(a,digits)
}
\arguments{
\item{a}{A numerical value to be rounded}
\item{digits}{Number of decimal places}
}
\value{
\itemize{
The rounded value is returned as a character.
}
}
\description{
A function that returns a rounded value as a character.
}
\examples{
rdc(2.412, 3)
rdc(2.41, 3)
rdc(2.4, 3)
rdc(2, 3)

rdc(-2.41, 3)
rdc(-2.4, 3)
rdc(-2, 3)

rdc(0, 3)
}
