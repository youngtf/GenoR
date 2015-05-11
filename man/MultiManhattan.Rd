\name{MultiManhattan}
\alias{MultiManhattan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MultiManhattan
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MultiManhattan(res, mapinfo, ylim0 = NULL, ylab = NULL, gap = NULL, ltype = "h", title = NULL, sub = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{mapinfo}{
%%     ~~Describe \code{mapinfo} here~~
}
  \item{ylim0}{
%%     ~~Describe \code{ylim0} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
}
  \item{ltype}{
%%     ~~Describe \code{ltype} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
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
function (res, mapinfo, ylim0 = NULL, ylab = NULL, gap = NULL, 
    ltype = "h", title = NULL, sub = NULL) 
{
    res.reorder = mapReorder(mapinfo)
    mapinfo = res.reorder[[1]]
    chromsort = res.reorder[[2]]
    if (is.null(gap)) {
        if (ltype == "h") 
            gap = 2e+07
        else gap = 5e+07
    }
    n.chrom = length(chromsort)
    chromlen0 = tapply(mapinfo[, 3], mapinfo[, 2], max) + 1
    chromlen = chromlen0[match(chromsort, names(chromlen0))]
    chrom.start = cumsum(c(1, (chromlen))) + gap * (0:(n.chrom))
    x.max = chrom.start[n.chrom + 1] + gap
    chrom.mid = chrom.start[-(n.chrom + 1)] + chromlen/2
    nres = ncol(res) - 1
    idx.res = match(res[, 1], mapinfo[, 1])
    idx.res.left = which(!is.na(idx.res))
    idx.res.right = idx.res[idx.res.left]
    res.all = matrix(0, nrow(mapinfo), nres)
    for (icol in 1:nres) {
        res.all[idx.res.right, icol] = res[idx.res.left, (icol + 
            1)]
    }
    data.plot = cbind(mapinfo, res.all)
    parBackup = par()
    par(mfrow = c(1, nres), mar = c(4, 2, 4, 1))
    for (i in 1:nres) {
        res = data.plot[, (3 + i)]
        if (is.null(ylim0)) {
            ylim = c(0, max(res) * 1.1)
        }
        else ylim = ylim0
        if (is.null(ylab)) 
            ylabi = " "
        else ylabi = ylab[i]
        if (is.null(title)) 
            titlei = " "
        else titlei = title[i]
        if (is.null(sub)) 
            subi = " "
        else subi = sub[i]
        plot(0, 0, frame.plot = FALSE, type = "n", ylim = c(0, 
            x.max), xlim = ylim, xlab = ylabi, ylab = "Chromosome", 
            xaxt = "n", yaxt = "n", main = titlei, sub = subi)
        axis(side = 2, at = chrom.mid, labels = chromsort, cex.lab = 0.5)
        axis(side = 1)
        for (chr in 1:n.chrom) {
            idx.chrom = which(data.plot[, 2] == chromsort[chr])
            ys = chrom.start[chr] + data.plot[idx.chrom, 3]
            xs = data.plot[idx.chrom, 3 + i]
            colchr = c("dark blue", "cornflowerblue")[(chr\%\%2) + 
                1]
            if (ltype == "h") {
                lines(x = c(0, 0), y = c(min(ys), max(ys)), lwd = 2, 
                  col = colchr, )
                for (p in which(xs > 0)) {
                  lines(x = c(0, xs[p]), y = c(ys[p], ys[p]), 
                    col = colchr)
                }
            }
            else points(xs, ys, col = colchr, pch = 16)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
