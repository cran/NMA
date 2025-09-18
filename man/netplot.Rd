\name{netplot}
\alias{netplot}
\title{Generating a networkplot}
\usage{
netplot(x,text=TRUE,col="black",bg="blue",base.lwd=1,base.cex=1)
}
\arguments{
\item{x}{Output object of \code{setup}}
\item{text}{A logical value that specifies whether the treatment labels are added}
\item{col}{Outer circumferential color of the nodes (default: black)}
\item{bg}{Color of the node (default: blue)}
\item{base.lwd}{A parameter adjusting edge widths (default: 1)}
\item{base.cex}{A parameter adjusting node sizes (default: 1)}
}
\value{
A networkplot is produced.
}
\description{
Generating a networkplot. The sizes of the nodes and edges are proportional to the corresponding sample sizes of direct comparisons. 
}
\examples{
data(heartfailure)

hf2 <- setup(study=study,trt=trt,d=d,n=n,measure="OR",ref="Placebo",data=heartfailure)

netplot(hf2)                                      # default color and sizes
netplot(hf2,base.lwd=1.5,base.cex=1.5)            # change the sizes
netplot(hf2,col="red",bg="red")                   # change the color
netplot(hf2,text=FALSE)                           # without texts
}
