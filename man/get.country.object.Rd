\name{get.country.object}
\Rdversion{1.1}
\alias{get.country.object}

\title{
Accessing Country Information
}
\description{
The function returns an object containing country name, code and index. 
}
\usage{
get.country.object(country, meta = NULL, country.table = NULL, 
    index = FALSE)
}
\arguments{
  \item{country}{Country name, code or index. If it is an index, the argument \code{index} must be set to \code{TRUE}.}
  \item{meta}{Object of class \code{\link{bayesTFR.mcmc.meta}}. If it is not given, the argument \code{country.table} must be given.}
  \item{country.table}{A table containing columns \dQuote{name} and \dQuote{code} from which the country info can be extracted. Only relevant, if \code{meta} is \code{NULL}. }
  \item{index}{Logical determining if the argument \code{country} is an index.}
}
\details{
Given partial information about a country (i.e. having either name or code or index), the function returns  an object containing all three pieces of information. Only countries from the \code{\link{UN2008}} dataset are considered that are included in the simulations and predictions. Note that country index is an internal index used in various components of a \code{\link{bayesTFR.mcmc.meta}} object. Thus, it does not relate to the column \dQuote{Index} in \code{\link{UN2008}}.
}
\value{
List with components:
\item{name}{Country name}
\item{code}{Country code}
\item{index}{Country index}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{country.names}}, \code{\link{UN2008}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
m <- get.tfr.mcmc(sim.dir)
# all three calls should give the same answer
get.country.object('China', m$meta)
get.country.object(156, m$meta)
get.country.object(56, m$meta, index=TRUE)
}

\keyword{ attribute }

