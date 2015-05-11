\name{mapReorder}
\alias{mapReorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mapReorder
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mapReorder(mapinfo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mapinfo}{
%%     ~~Describe \code{mapinfo} here~~
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
function (mapinfo) 
{
    chromname = as.character(unique(mapinfo[, 2]))
    isnum = chromname \%in\% as.character(0:100)
    if (any(!isnum)) {
        chromsort = c(sort(as.integer(chromname[isnum])), chromname[!isnum])
        mapinfo = mapinfo[order(match(mapinfo[, 2], chromsort), 
            mapinfo[, 3]), ]
    }
    else {
        chromsort = as.character(sort(as.integer(chromname)))
        mapinfo = mapinfo[order(mapinfo[, 2], mapinfo[, 3]), 
            ]
    }
    list(mapinfo, chromsort)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
