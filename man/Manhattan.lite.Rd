\name{Manhattan.lite}
\alias{Manhattan.lite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Manhattan.lite
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Manhattan.lite(res, mapinfo, gap = 5e+07, ylim = NULL, ylab = NULL, title = NULL, ltype = "p", cols = c("dark blue", "cornflowerblue"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{mapinfo}{
%%     ~~Describe \code{mapinfo} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
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
  \item{ltype}{
%%     ~~Describe \code{ltype} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
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
function (res, mapinfo, gap = 5e+07, ylim = NULL, ylab = NULL, 
    title = NULL, ltype = "p", cols = c("dark blue", "cornflowerblue")) 
{
    if (is.null(ylim)) 
        ylim = c(0, max(res[, 2]) * 1.1)
    res.reorder = mapReorder(mapinfo)
    mapinfo = res.reorder[[1]]
    chromsort = res.reorder[[2]]
    res.frame = DrawAFrame(res.reorder, type = "pos", ylim = ylim, 
        ylab = ylab, title = title, gap = gap)
    idx.res = match(res[, 1], mapinfo[, 1])
    idx.res.left = which(!is.na(idx.res))
    idx.res.right = idx.res[idx.res.left]
    res.all = rep(0, nrow(mapinfo))
    res.all[idx.res.right] = res[idx.res.left, 2]
    data.plot = cbind(mapinfo, res.all)
    for (chr in 1:res.frame$n.chrom) {
        idx.chrom = which(data.plot[, 2] == chromsort[chr])
        xs = res.frame$chrom.start[chr] + data.plot[idx.chrom, 
            3]
        ys = data.plot[idx.chrom, 4]
        points(xs, ys, col = cols[(chr\%\%2) + 1], pch = 16, type = ltype)
    }
    axis(side = 1, at = res.frame$chrom.mid, labels = chromsort)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
