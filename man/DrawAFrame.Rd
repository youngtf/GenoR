\name{DrawAFrame}
\alias{DrawAFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{DrawAFrame
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DrawAFrame(res.reorder, type = c("pos", "idx"), ylim = NULL, ylab = NULL, title = NULL, gap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.reorder}{
%%     ~~Describe \code{res.reorder} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
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
function (res.reorder, type = c("pos", "idx"), ylim = NULL, ylab = NULL, 
    title = NULL, gap) 
{
    match.arg(arg = type, choices = c("pos", "idx"))
    mapinfo = res.reorder[[1]]
    chromsort = res.reorder[[2]]
    n.chrom = length(chromsort)
    if (type == "pos") {
        chromlen0 = tapply(mapinfo[, 3], mapinfo[, 2], max) + 
            1
    }
    else if (type == "idx") {
        chromlen0 = tapply(mapinfo[, 3], mapinfo[, 2], length)
    }
    chromlen = chromlen0[match(chromsort, names(chromlen0))]
    chrom.start = cumsum(c(1, (chromlen))) + gap * (0:(n.chrom))
    x.max = chrom.start[n.chrom + 1] + gap
    chrom.mid = chrom.start[-(n.chrom + 1)] + chromlen/2
    if (is.null(ylab)) 
        ylab = "-log(p-val)"
    if (is.null(title)) 
        title = "Model"
    plot(0, 0, type = "n", xlim = c(0, x.max), ylim = ylim, ylab = ylab, 
        xlab = "Chromosome", xaxt = "n", main = title)
    res = list(n.chrom = n.chrom, chrom.start = chrom.start, 
        chrom.mid = chrom.mid)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
