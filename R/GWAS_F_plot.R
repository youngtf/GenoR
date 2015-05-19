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
#------------------------------------------------------------------------------
# Function:    MultiManhattan(res,mapinfo,ylim0 = NULL,ylab = NULL
#                            ,gap  = NULL,ltype= "h",title= NULL,sub = NULL)
# Description: Manhattan plot for comparison of results of GWAS
#------------------------------------------------------------------------------
# MultiManhattan = function(res
#                           , mapinfo
#                           , ylim0 = NULL
#                           , ylab = NULL
#                           , gap  = NULL
#                           , ltype= "h"
#                           , title= NULL
#                           , sub = NULL
#                           ){
#     ### order mapinfo
#     res.reorder = mapReorder(mapinfo)
#     mapinfo   = res.reorder[[1]]
#     chromsort = res.reorder[[2]]
#     ### X axis  
#     if(is.null(gap)){
#         if (ltype == "h") gap = 20000000 else gap = 50000000
#     }
#     n.chrom   = length(chromsort)
#     chromlen0 = tapply(mapinfo[,3],mapinfo[,2],max)+1
#     chromlen  = chromlen0[match(chromsort,names(chromlen0))]
#     chrom.start = cumsum(c(1,(chromlen))) + gap * (0:(n.chrom))
#     x.max     	= chrom.start[n.chrom+1] + gap
#     chrom.mid   = chrom.start[-(n.chrom+1)] + chromlen/2
#     
#     nres = ncol(res) - 1
#     
#     ### prepare plotting data
#     idx.res = match(res[,1],mapinfo[,1])
#     idx.res.left  = which(!is.na(idx.res))
#     idx.res.right = idx.res[idx.res.left]
#     res.all = matrix(0,nrow(mapinfo),nres)
#     for (icol in 1: nres){
#         res.all[idx.res.right,icol] = res[idx.res.left,(icol+1)]
#     }
#     data.plot = cbind(mapinfo,res.all)
#     
#     ### Big frame
#     parBackup = par()
#     par(mfrow = c(1,nres)
#         , mar = c(4,2,4,1))
#     
#     for(i in 1:nres){
#         res = data.plot[,(3+i)]
#         if(is.null(ylim0)){
#             ylim = c(0,max(res)*1.1)
#         }else ylim = ylim0
#         if(is.null(ylab))  ylabi = " "  else ylabi = ylab[i]
#         if(is.null(title)) titlei = " " else titlei = title[i]
#         if(is.null(sub))   subi = " "   else subi = sub[i]
#         plot(0,0,frame.plot=FALSE,
#              type="n",ylim=c(0,x.max),xlim=ylim,
#              xlab=ylabi,ylab="Chromosome",xaxt="n",yaxt="n",
#              main = titlei,sub=subi)
#         axis(side=2,at=chrom.mid,labels=chromsort,cex.lab=0.5)
#         axis(side=1)
#         for (chr in 1:n.chrom){
#             idx.chrom  = which(data.plot[,2] == chromsort[chr])
#             ys = chrom.start[chr]+data.plot[idx.chrom,3]
#             xs = data.plot[idx.chrom,3+i]
#             colchr = c("dark blue","cornflowerblue")[(chr%%2)+1]
#             
#             if(ltype=="h"){
#                 lines(x=c(0,0),y=c(min(ys),max(ys)),lwd=2
#                       ,col=colchr,)
#                 for(p in which(xs>0)){
#                     lines(x=c(0,xs[p]),y=c(ys[p],ys[p]),col=colchr)
#                 }
#             }else points(xs,ys,col=colchr,pch=16)
#         }
#     }
# }
#------------------------------------------------------------------------------
# Updated Mar 17, 2015 16:41
# Function:     
# Description:  
# input:        
# ouput:        
#------------------------------------------------------------------------------
# Manhattan.D = function(res, mapinfo, gap  = 50000000
#                        , ylim = NULL, ylab = NULL, title= NULL
#                        , ltype= "p", cols = c("dark blue","cornflowerblue")
#                        , cols2 = c("darkred","red")
# ){
#   
#   if(is.null(ylim))  ylim = c(0,max(res[,(2:3)])*1.1)
#   
#   res.reorder = mapReorder(mapinfo)
#   mapinfo   = res.reorder[[1]]
#   chromsort = res.reorder[[2]]
#   
#   DrawAFrame = function(res.reorder,typre="pos",
#                         ylim=ylim,ylab=ylab,title=title,gap=gap)
#     
#     ### prepare plotting data
#     idx.res = match(res[,1],mapinfo[,1])
#   idx.res.left  = which(!is.na(idx.res))
#   idx.res.right = idx.res[idx.res.left]
#   res.all = matrix(0,nrow(mapinfo),2)
#   res.all[idx.res.right,1] = res[idx.res.left,2]
#   res.all[idx.res.right,2] = res[idx.res.left,3]
#   data.plot = cbind(mapinfo,res.all)
#   
#   ### points
#   for (chr in 1:n.chrom){
#     idx.chrom  = which(data.plot[,2] == chromsort[chr])
#     xs = chrom.start[chr]+data.plot[idx.chrom,3]
#     ys = data.plot[idx.chrom,4]
#     ysd = data.plot[idx.chrom,5]
#     points(xs,ys
#            ,col=cols[(chr%%2)+1]
#            ,pch=16,type=ltype)
#     points(xs,ysd
#            ,col=cols2[(chr%%2)+1]
#            ,pch=16,type=ltype)
#   }
#   ### axis
#   # left axis
#   axis(side=1,at=chrom.mid,labels=chromsort)
#   
# }
#------------------------------------------------------------------------------
# Function:    
# Description: 
#------------------------------------------------------------------------------
mapReorder = function(mapinfo){
    # first column is marker name, second is chromosome, third is map position.
    chromname = as.character(unique(mapinfo[,2]))
    isnum = chromname %in% as.character(0:100)
    if (any(!isnum)){
        chromsort = c(sort(as.integer(chromname[isnum]))
                      ,chromname[!isnum])
        mapinfo   = mapinfo[order(match(mapinfo[,2],chromsort),mapinfo[,3]),]
    }else{
        chromsort = as.character(sort(as.integer(chromname)))
        mapinfo   = mapinfo[order(mapinfo[,2],mapinfo[,3]),]
    }
    list(mapinfo,chromsort)
}


# -----------------------------------------------------------------------------
# Updated May 18, 2015 9:50 PM
# Function:     
# Description:  
# input:        
# ouput:        
# -----------------------------------------------------------------------------
CreateWindowMap = function(Map.single, distance){
  
  ## Reorder the map
  Map.single = mapReorder(Map.single)[[1]]
  
  n.chrom = length(unique(Map.single[,2]))
  pos.start = rep(0,n.chrom)
  pos.end = tapply(Map.single[,3],Map.single[,2],max)
  
  n.windows = ceiling(pos.end / distance)
  
  info.chr = rep(1:n.chrom,times =n.windows)
  window.name = unlist(lapply(n.windows,function(x) seq(1:x)),use.names = F)
  info.pos = window.name - 0.5
  info.ID = paste0(info.chr,"_",window.name)
  
  res = data.frame(ID = info.ID, chr = info.chr, pos = info.pos,stringsAsFactors = FALSE)
}
  
  
