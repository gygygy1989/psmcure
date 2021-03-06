% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/age_base.R
\name{split_base}
\alias{split_base}
\title{Semiparametric baseline estimate stratified by age.}
\usage{
split_base(Time, Status, X, beta, w, mort.s, mort.h)
}
\arguments{
\item{Time}{observed survival time}

\item{Status}{right censored status (1=death, 0=censored)}

\item{X}{covariates in latency survival function}

\item{beta}{parameters for X}

\item{w}{unobserved cure indicator}

\item{mort.s}{background survival rate}

\item{mort.h}{background hazard rate}
}
\value{
None
}
\description{
Semiparametric baseline estimate stratified by age.
}
