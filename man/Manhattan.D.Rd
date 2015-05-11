\name{Manhattan.D}
\alias{Manhattan.D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Manhattan.D
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Manhattan.D(res, mapinfo, gap = 5e+07, ylim = NULL, ylab = NULL, title = NULL, ltype = "p", cols = c("dark blue", "cornflowerblue"), cols2 = c("darkred", "red"))
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
  \item{cols2}{
%%     ~~Describe \code{cols2} here~~
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
    title = NULL, ltype = "p", cols = c("dark blue", "cornflowerblue"), 
    cols2 = c("darkred", "red")) 
{
    if (is.null(ylim)) 
        ylim = c(0, max(res[, (2:3)]) * 1.1)
    res.reorder = mapReorder(mapinfo)
    mapinfo = res.reorder[[1]]
    chromsort = res.reorder[[2]]
    DrawAFrame = function(res.reorder, typre = "pos", ylim = ylim, 
        ylab = ylab, title = title, gap = gap) idx.res = match(res[, 
        1], mapinfo[, 1])
    idx.res.left = which(!is.na(idx.res))
    idx.res.right = idx.res[idx.res.left]
    res.all = matrix(0, nrow(mapinfo), 2)
    res.all[idx.res.right, 1] = res[idx.res.left, 2]
    res.all[idx.res.right, 2] = res[idx.res.left, 3]
    data.plot = cbind(mapinfo, res.all)
    for (chr in 1:n.chrom) {
        idx.chrom = which(data.plot[, 2] == chromsort[chr])
        xs = chrom.start[chr] + data.plot[idx.chrom, 3]
        ys = data.plot[idx.chrom, 4]
        ysd = data.plot[idx.chrom, 5]
        points(xs, ys, col = cols[(chr\%\%2) + 1], pch = 16, type = ltype)
        points(xs, ysd, col = cols2[(chr\%\%2) + 1], pch = 16, 
            type = ltype)
    }
    axis(side = 1, at = chrom.mid, labels = chromsort)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
