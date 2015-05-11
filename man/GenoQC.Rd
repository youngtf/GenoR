\name{GenoQC}
\alias{GenoQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GenoQC
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GenoQC(genodata)
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
    nmar = ncol(genodata)
    nind = nrow(genodata)
    res = matrix(0, 3, nmar)
    rownames(res) = c("Missing", "MAF", "HWtest")
    res[1, ] = colSums(is.na(genodata))/nind
    idx.allmissing = which(res[1, ] == 1)
    genodata.info = genodata[, -idx.allmissing]
    AlleAMat = (genodata.info == 1) * 2 + (genodata.info == 0)
    callMat = (!is.na(genodata.info)) * 2
    AAF = colSums(AlleAMat, na.rm = T)/colSums(callMat)
    BAF = 1 - AAF
    res[2, -idx.allmissing] = ifelse((AAF <= 0.5), AAF, BAF)
    nValid = rowSums((!is.na(genodata.info)))
    if (any(nValid == 0)) 
        stop("allmissing unremoved, debug the program.")
    nAA_AB_BB = cbind(colSums((genodata.info == 1), na.rm = T), 
        colSums((genodata.info == 0), na.rm = T), colSums((genodata.info == 
            -1), na.rm = T))
    EnAA_AB_BB = cbind((AAF^2) * nValid, (AAF * BAF * 2) * nValid, 
        (BAF^2) * nValid)
    res[3, -idx.allmissing] = colSums(((nAA_AB_BB - EnAA_AB_BB)^2)/EnAA_AB_BB)
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
