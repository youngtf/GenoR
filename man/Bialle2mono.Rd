\name{Bialle2mono}
\alias{Bialle2mono}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bialle2mono
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Bialle2mono(genodata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genodata}{
%%     ~~Describe \code{genodata} here~~
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
function (genodata) 
{
    nind = nrow(genodata)
    nmar = ncol(genodata)
    if (is.null(rownames(genodata))) {
        warning("no valid rownames for the input data")
        rownames(genodata) = paste0("R", 1:nind)
    }
    matAlle = matrix(" ", nind * 2, nmar)
    matAlle[1:nind, ] = substr(genodata, 1, 1)
    matAlle[(nind + 1):(2 * nind), ] = substr(genodata, 2, 2)
    rownames(matAlle) = c(paste0(rownames(genodata), ".1"), paste0(rownames(genodata), 
        ".2"))
    colnames(matAlle) = colnames(genodata)
    matAlle
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
