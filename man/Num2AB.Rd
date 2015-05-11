\name{Num2AB}
\alias{Num2AB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Num2AB
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Num2AB(genodata, code = c(1, 0, -1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genodata}{
%%     ~~Describe \code{genodata} here~~
}
  \item{code}{
%%     ~~Describe \code{code} here~~
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
function (genodata, code = c(1, 0, -1)) 
{
    genodata[(genodata == code[3])] = "BB"
    genodata[(genodata == code[2])] = "AB"
    genodata[(genodata == code[1])] = "AA"
    genodata
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
