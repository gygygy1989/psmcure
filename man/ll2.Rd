% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ll2.R
\name{ll2}
\alias{ll2}
\title{Likelihood function for latency survival function in EM algorithm}
\usage{
ll2(trybeta, Time, Status, age, X, w, mort.s, mort.h)
}
\arguments{
\item{trybeta}{The estimated parameters in latency survival function.}

\item{Time}{Survival time from data set.}

\item{Status}{Right censoring status: 1 is death, 0 is censored.}

\item{age}{the variable of age}

\item{X}{The covaraites for latency survival function.}

\item{w}{unobserved cure status indicator.}

\item{mort.s}{background survival rate from life table.}

\item{mort.h}{background hazard rate from life table.}
}
\value{
None
}
\description{
Likelihood function for latency survival function in EM algorithm
}
