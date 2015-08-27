#------------------------------------------------------------------------------
# R file  --Tianfu Yang
# Type:              Functions
# Subtype/Project:   Ploting
# Descripetions:     functions used in ploting
# Last Update:       2014-05-20
# Contents:
# 1
#   Function:    Manhattan.lite(res,mapinfo,ylim=NULL,xlab=NULL
#                             ,gap=50000000,ltype="p",title= NULL)
#   Description: Manhattan plot for a single result of GWAS
# 2
#   Function:    MultiManhattan(res,mapinfo,ylim0 = NULL,ylab = NULL
#                               ,gap  = NULL,ltype= "h",title= NULL,sub = NULL)
#   Description: Manhattan plot for comparison of results of GWAS
# 3
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Analysis Block
# Type:              Functions
# Descripetions:     
# Last Update:       Nov 10, 2014
#------------------------------------------------------------------------------
DrawAFrame = function(res.reorder
                         ,type=c("pos","idx")
                         ,ylim=NULL,ylab=NULL,title=NULL,gap){
  match.arg(arg = type,choices = c("pos","idx"))
  ### order mapinfo
  mapinfo   = res.reorder[[1]]
  chromsort = res.reorder[[2]]
  ### build the frame    
  n.chrom   = length(chromsort)
  if (type == "pos"){
    chromlen0 = tapply(mapinfo[,3],mapinfo[,2],max)+1
  } else if (type == "idx"){
    chromlen0 = tapply(mapinfo[,3],mapinfo[,2],length)
  }
  chromlen  = chromlen0[match(chromsort,names(chromlen0))]
  chrom.start = cumsum(c(1,(chromlen))) + gap * (0:(n.chrom))
  x.max       = chrom.start[n.chrom+1] + gap
  chrom.mid   = chrom.start[-(n.chrom+1)] + chromlen/2
  
  if(is.null(ylab))  ylab = "-log(p-val)"
  if(is.null(title)) title = "Model"
  
  plot(0,0,
       type="n",xlim=c(0,x.max),ylim=ylim,
       ylab=ylab,xlab="Chromosome",xaxt="n",
       main = title)
  res = list(n.chrom = n.chrom, chrom.start = chrom.start,chrom.mid=chrom.mid)
}

#------------------------------------------------------------------------------
# Function:    Manhattan.lite(res,mapinfo,ylim=NULL,xlab=NULL
#                             ,gap=50000000,ltype="p",title= NULL)
# Description: Manhattan plot for a single result of GWAS
#------------------------------------------------------------------------------
Manhattan.lite = function(res
                          , mapinfo
                          , gap  = 50000000   # frame
                          , ylim = NULL       # frame
                          , ylab = NULL       # frame
                          , title= NULL       # frame
                          , ltype= "p"        # points
                          , cols = c("dark blue","cornflowerblue") # points
                          , axes = TRUE
){    
  if(is.null(ylim))  ylim = c(0,max(res[,2])*1.1)
  
  res.reorder = mapReorder(mapinfo)
  mapinfo   = res.reorder[[1]]
  chromsort = res.reorder[[2]]
  
  res.frame = DrawAFrame(res.reorder,type="pos",
                         ylim=ylim,ylab=ylab,title=title,gap=gap)
  
  ### prepare plotting data
  idx.res = match(res[,1],mapinfo[,1])
  idx.res.left  = which(!is.na(idx.res))
  idx.res.right = idx.res[idx.res.left]
  res.all = rep(0,nrow(mapinfo))
  res.all[idx.res.right] = res[idx.res.left,2]
  data.plot = cbind(mapinfo,res.all)
  
  for (chr in 1:res.frame$n.chrom){                     ### points
    idx.chrom  = which(data.plot[,2] == chromsort[chr])
    xs = res.frame$chrom.start[chr]+data.plot[idx.chrom,3]
    ys = data.plot[idx.chrom,4]
    points(xs,ys,col=cols[(chr%%2)+1],pch=16,type=ltype)
  }
  if (axes)
    axis(side=1,at=res.frame$chrom.mid,labels=chromsort) ### left axis
  
  res.frame
}


Manhattan.win = function(res, mapinfo, gap  = 50
                         , ylim = NULL, ylab = NULL, title= NULL
                         , ltype= "p", cols = NULL
                     ){    
  # res is a list with multiple "res"
  
  if(is.null(ylim)){
    ylim = c(0,max(unlist(lapply(res,function(x) max(x[,2]))))*1.1)
  }
  res.reorder = mapReorder(mapinfo)
  mapinfo   = res.reorder[[1]]
  chromsort = res.reorder[[2]]
  nvalues = length(res)
  res.frame = DrawAFrame(res.reorder,type="pos",
                         ylim=ylim,ylab=ylab,title=title,gap=gap)
    
  if (is.null(cols) || length(cols) < nvalues) cols = seq(nvalues)
  
  ### prepare plotting data
  
  for (ivalue in 1: nvalues){
    idx.res = match(res[[ivalue]][,1],mapinfo[,1])
    idx.res.left  = which(!is.na(idx.res))
    idx.res.right = idx.res[idx.res.left]
  
    res.all = rep(0,nrow(mapinfo))
    res.all[idx.res.right] = as.matrix(res[[ivalue]][idx.res.left,2])
    data.plot = cbind(mapinfo,res.all)
  
    for (chr in 1:res.frame$n.chrom){                     ### points
      idx.chrom  = which(data.plot[,2] == chromsort[chr])
      xs = res.frame$chrom.start[chr]+data.plot[idx.chrom,3]
      ys = data.plot[idx.chrom,4]
      points(xs,ys,col=cols[ivalue],pch=16,type=ltype)
    }
  }
  axis(side=1,at=res.frame$chrom.mid,labels=chromsort) ### left axis
}


# -----------------------------------------------------------------------------
# Analysis Block
# Type:              Analysis/Visualization/Functions/UnitTest
# Descriptions:      A function to plot lsmeans results
# Last Update:       Jun 29, 2015
# -----------------------------------------------------------------------------
plot.lsmeans = function(res.lm,SNPnames
                        ,nSNPs = length(SNPnames)
                        ,mfrow = c(1,nSNPs)
                        ,mar   = c(4.1,5.1,2.1,1.1)
                        ,ylim  = c(2,6)
                        ,ylab = "lsmean"
                        ,xlab = NULL){
  # check parameters
  if (!is.null(xlab)){
    if(length(xlab) != nSNPs){
      stop("The length of xlab is different from that of SNPs.")  
    }
  }
  plot.new()
  par(mfrow = mfrow)
  par(mar = mar)
  
  for(i in seq(nSNPs)){
    res.lsm = summary(lsmeans(res.lm,SNPnames[i]))
    v.lsm = res.lsm$lsmean
    v.lower = res.lsm$lower.CL
    v.upper = res.lsm$upper.CL
    v.name = names(res.lsm)[1]
    if (is.null(xlab)){
      x_lab = v.name
    } else {
      x_lab = xlab[i]
    }
    n.level = length(v.lower)
    name.levels = levels(res.lsm[[1]])
    if (i == 1){
      plot(0,col="white",ylim = ylim, xlim =c(0.5,n.level + 0.5),
           xlab = "", ylab = ylab,axes = F,frame.plot =T)
    } else {
      plot(0,col="white",ylim = ylim, xlim =c(0.5,n.level + 0.5),
           xlab = "", ylab = "",axes = F,frame.plot =T)
    }
    title(xlab = x_lab,line = 4)
    axis(2)
    axis(1,labels = c(name.levels),at = 1:n.level)
    
    points(x = 1:n.level,y=v.lsm)
    for (j in 1:n.level){
      lines(x=c(j,j),y=c(v.lower[j],v.upper[j]))
      lines(x=c(j-0.1,j+0.1),y=c(v.lower[j],v.lower[j]))
      lines(x=c(j-0.1,j+0.1),y=c(v.upper[j],v.upper[j]))
    }
    
    # table
    genotype.snp = as.character(res.lm$model[,v.name])
    for (j in 1:n.level){
      text(j,ylim[1],sum(genotype.snp == name.levels[j]))
    }
  }
}
