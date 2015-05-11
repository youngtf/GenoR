\name{strsplit.mat}
\alias{strsplit.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{strsplit.mat
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
strsplit.mat(vec.char, sep = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec.char}{
%%     ~~Describe \code{vec.char} here~~
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
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
function (vec.char, sep = " ") 
{
    res.list = strsplit(vec.char, sep)
    len.char = lapply(res.list, length)
    if (length(unique(len.char)) > 1) 
        stop("not the same pattern")
    res.vec.char = unlist(res.list)
    matrix(res.vec.char, ncol = unique(len.char), by.row = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
