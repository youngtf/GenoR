\name{Mono2bialle}
\alias{Mono2bialle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mono2bialle
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Mono2bialle(alleledata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alleledata}{
%%     ~~Describe \code{alleledata} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (alleledata) 
{
    nmar = ncol(alleledata)
    if (nrow(alleledata)\%\%2 == 1) 
        stop("wrong nrow!")
    nind = nrow(alleledata)/2
    indnames = rownames(alleledata)
    nchar.indnames = nchar(indnames)
    if (any(nchar.indnames[1:nind] != nchar.indnames[(nind + 
        1):(2 * nind)])) {
        stop("wrong indnames")
    }
    indnames1 = substr(indnames[1:nind], 1, nchar.indnames[1:nind] - 
        1)
    indnames2 = substr(indnames[(nind + 1):(2 * nind)], 1, nchar.indnames[1:nind] - 
        1)
    if (any(indnames1 != indnames2)) 
        stop("wrong indnames")
    Allele1 = alleledata[1:nind]
    Allele2 = alleledata[(nind + 1):(2 * nind)]
    res.bia = paste0(Allele1, Allele2)
    dim(res.bia) = c(nind, nmar)
    res.bia
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
