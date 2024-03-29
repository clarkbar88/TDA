% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TDA_Power.R
\name{TDA_Power}
\alias{TDA_Power}
\title{Complete TDA power analyis}
\usage{
TDA_Power(value, year, site, trend.magnitude.threshold, power = 0.8,
  power.simulations.n = 1000)
}
\arguments{
\item{value}{vector of values (concentration, BIBI, etc.)}

\item{year}{vector of years - if multiple values for same year, values are averaged (arithmetically)}

\item{site}{vector of sites}

\item{trend.magnitude.threshold}{absolute value of ecologically meaningful trend magnitude. if not provided, defaults to 10 percent of average value}

\item{power}{1-B - defaults to 0.8}

\item{power.simulation}{how many permutations to run to calculate power - default to 1000}
}
\value{
A tibble of trends and power.
}
\description{
Complete TDA power analyis
}
\examples{
data(all_bibi_data)
tmp<-subset(all_bibi_data,Site \%in\% c('08BEA3312','08BEA3321'))
TDA_Kendall(tmp$`Overall Score`,tmp$Year,tmp$Sites,0.5)
}
