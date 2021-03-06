% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plotIsingGraph}
\alias{plotIsingGraph}
\title{plotIsingGraph}
\usage{
plotIsingGraph(x, weight = 0.6, layout = "kk")
}
\arguments{
\item{x}{\code{flowReMix} object.}
}
\value{
\code{ggplot} ggraph object
}
\description{
Plot an ising graph.
}
\examples{
data(fit505)
plotIsingGraph(fit505,weight=0.9) +
guides(shape=guide_legend(nrow=1),
  size=guide_legend(nrow=1),
  color=guide_legend(nrow=4))
}
