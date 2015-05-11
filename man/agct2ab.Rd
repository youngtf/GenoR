\name{agct2ab}
\alias{agct2ab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{agct2ab
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
agct2ab(genodata, char.na = "-", allAlleles = c("A", "G", "C", "T"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genodata}{
%%     ~~Describe \code{genodata} here~~
}
  \item{char.na}{
%%     ~~Describe \code{char.na} here~~
}
  \item{allAlleles}{
%%     ~~Describe \code{allAlleles} here~~
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
function (genodata, char.na = "-", allAlleles = c("A", "G", "C", 
    "T")) 
{
    nmar = ncol(genodata)
    nind = nrow(genodata)
    cat("Fragmentation...\n")
    genodata.mono = Bialle2mono(genodata)
    genodpata.mono[genodata.mono == char.na] == NA
    cat("Transformation...\n")
    res.mat = genodata.mono
    codingSheet = matrix("0", 2, nmar)
    pb <- txtProgressBar(min = 0, max = nmar, style = 3)
    for (i in 1:nmar) {
        uniAllele = unique(genodata.mono[, i])
        isvalid = which(!is.na(uniAllele))
        uniAllele.nona = uniAllele[isvalid]
        npoly = length(isvalid)
        if (npoly == 1) {
            res.mat[, i] = sub(uniAllele.nona, "A", genodata.mono[, 
                i])
            codingSheet[, i] = c(uniAllele.nona, uniAllele.nona)
        }
        else if (npoly == 2) {
            res.mat[, i] = sub(uniAllele.nona[1], "A", genodata.mono[, 
                i])
            res.mat[, i] = sub(uniAllele.nona[2], "B", genodata.mono[, 
                i])
            codingSheet[, i] = uniAllele.nona
        }
        else {
            res.mat[, i] = rep(NA, nind * 2)
        }
        setTxtProgressBar(pb, i)
    }
    res.mat
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
