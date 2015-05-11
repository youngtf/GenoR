\name{blank.remover}
\alias{blank.remover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{blank.remover
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
blank.remover(vec.char)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec.char}{
%%     ~~Describe \code{vec.char} here~~
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
function (vec.char) 
{
    blank.rm.ele = function(char) {
        n = nchar(char)
        while (substr(char, n, n) == " ") {
            char = substr(char, 1, (n - 1))
            n = nchar(char)
        }
        char
    }
    res.char = vec.char
    for (i in 1:length(vec.char)) {
        res.char[i] = blank.rm.ele(vec.char[i])
    }
    res.char
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
