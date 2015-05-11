\name{mono2PED}
\alias{mono2PED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mono2PED
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mono2PED(monodata, FamilyID = NULL, IndID = NULL, PID = NULL, MID = NULL, Sex = NULL, Pheno = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{monodata}{
%%     ~~Describe \code{monodata} here~~
}
  \item{FamilyID}{
%%     ~~Describe \code{FamilyID} here~~
}
  \item{IndID}{
%%     ~~Describe \code{IndID} here~~
}
  \item{PID}{
%%     ~~Describe \code{PID} here~~
}
  \item{MID}{
%%     ~~Describe \code{MID} here~~
}
  \item{Sex}{
%%     ~~Describe \code{Sex} here~~
}
  \item{Pheno}{
%%     ~~Describe \code{Pheno} here~~
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
function (monodata, FamilyID = NULL, IndID = NULL, PID = NULL, 
    MID = NULL, Sex = NULL, Pheno = NULL) 
{
    if (nrow(monodata)\%\%2 == 1) 
        stop("wrong nrow!")
    nind = nrow(monodata)/2
    if (is.null(FamilyID)) 
        FamilyID = paste0("Family_", seq(nind))
    if (is.null(PID)) 
        PID = paste0("Sire_", seq(nind))
    if (is.null(MID)) 
        MID = paste0("Dam_", seq(nind))
    if (is.null(Sex)) 
        Sex = rep(1, nind)
    if (is.null(Pheno)) 
        Pheno = rep(1, nind)
    if (is.null(IndID)) {
        IndID.Raw = rownames(monodata)[1:nind]
        len.names = nchar(IndID.Raw)
        IndID = substr(IndID.Raw, 1, (len.names - 2))
    }
    res.PED = data.frame(FamilyID = FamilyID, IndID = IndID, 
        PID = PID, MID = MID, Sex = Sex, Pheno = Pheno, matrix(monodata, 
            nrow = nind))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
