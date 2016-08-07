#------------------------------------------------------------------------------
# R file  --Tianfu Yang
# Type:              Functions
# Subtype/Project:   Ploting
# Descripetions:     functions used in ploting
# Last Update:       Aug 27, 2015 5:10 PM
# Contents:
# 1
#   Function:    DrawAFrame(map.sorted, gap, ylim, axes = FALSE,
#                           xlab = "", ylab = "", title = "")
# 2
#   Function:    Manhattan(res.gwas.ID, res.gwas.value, map.sorted,
#                          gap  = 50000000, ylim = NULL, xlab = "Chromosome", 
#                          ylab = "", title= "", axes = TRUE, ltype = "p", 
#                          pch = 16, cols = c("dark blue","cornflowerblue"))
# 3 
#   Function:    PlotLSmeans(res.lm,SNPnames,    nSNPs = length(SNPnames),
#                            mfrow = c(1,nSNPs), mar   = c(4.1,5.1,2.1,1.1)
#                            ylim  = c(2,6),     ylab  = "lsmean",
#                            xlab  = NULL)
# -----------------------------------------------------------------------------
# To-do list:
# 1 suppress message in Manhattan                               Done!
# 2 clear possible warning in the join of tbls in Manhattan     Done!
# 3 Transfer the plotting system to ggplot2
# 4 Annotation tool for manhattan plot
# -----------------------------------------------------------------------------
# UPDATED Aug 27, 2015 5:35 PM
# FUNCTION:     DrawAFrame(map.sorted, gap, ylim, axes = FALSE,
#'                        xlab = "", ylab = "", title = "")
#' @title       Position calculation and frame plotting
#' @param       map.sorted A sorted map
#' @param       gap        The gap between chromosomes in the plot
#' @param       ylim       The y-axis range
#' @param       axes       Logical, draw the axes or not
#' @param       xlab,ylab  Parameters passed to plot.default()       
#' @param       title      Character, the title of the plot 
#' @return      A data frame of chromosome information
# -----------------------------------------------------------------------------
#' @note        This is a private function that calculate the genome positions
#'              and draw a basic frame for manhattan plots
#' @examples    
#' SID = paste("SNP",1:20,sep="_")
#' chromosome = c(rep("0",3),rep("1",5),rep("2",3),rep("X",6),rep("Y",3))
#' position = c(31, 71, 45, 60, 34, 41, 19, 23, 63, 43,  
#'               2, 79,  6 ,95, 39, 97, 15, 32, 56, 88)
#'               
#' (map.sorted = SortMap(SID, chromosome, position))
#' # $Chromosomes
#' #   chr.id chr.is.allosome chr.length
#' # 0      0           FALSE         71
#' # 1      1           FALSE         60
#' # 2      2           FALSE         63
#' # X      X            TRUE         97
#' # Y      Y            TRUE         88
#' # 
#' # $SNPs
#' # $SNPs$`0`
#' # ID chromosome position
#' # 1 SNP_1          0       31
#' # 3 SNP_3          0       45
#' # 2 SNP_2          0       71
#' # 
#' # $SNPs$`1`
#' # ID chromosome position
#' # 4 SNP_7          1       19
#' # 5 SNP_8          1       23
#' # 2 SNP_5          1       34
#' # 3 SNP_6          1       41
#' # 1 SNP_4          1       60
#' # 
#' # $SNPs$`2`
#' # ID chromosome position
#' # 3 SNP_11          2        2
#' # 2 SNP_10          2       43
#' # 1  SNP_9          2       63
#' # 
#' # $SNPs$X
#' # ID chromosome position
#' # 2 SNP_13          X        6
#' # 6 SNP_17          X       15
#' # 4 SNP_15          X       39
#' # 1 SNP_12          X       79
#' # 3 SNP_14          X       95
#' # 5 SNP_16          X       97
#' # 
#' # $SNPs$Y
#' # ID chromosome position
#' # 1 SNP_18          Y       32
#' # 2 SNP_19          Y       56
#' # 3 SNP_20          Y       88
#' # 
#' # 
#' # attr(,"class")
#' # [1] "geno.map"
#' DrawAFrame(map.sorted, gap = 30, ylim = c(1,3))
#' #   chr.id chr.is.allosome chr.length chr.start chr.end chr.mid
#' # 1      1           FALSE         60         1      61    31.0
#' # 2      2           FALSE         63        91     154   122.5
#' # X      X            TRUE         97       184     281   232.5
#' # Y      Y            TRUE         88       311     399   355.0
# -----------------------------------------------------------------------------

DrawAFrame = function(map.sorted, gap, ylim, axes = FALSE,
                      xlab = "", ylab = "", title = ""){
  # horizontal axis
    chr.info = map.sorted$Chromosomes
    if (chr.info$chr.id[1] == "0"){
      chr.info = chr.info[-1,]
    }
    n.chr    = nrow(chr.info)
    chr.info$chr.start = cumsum(c(1,(chr.info$chr.length)))[1:(n.chr)] + 
                         gap * (0:(n.chr-1))
    chr.info$chr.end   = chr.info$chr.start + chr.info$chr.length
    chr.info$chr.mid   = (chr.info$chr.start + chr.info$chr.end) / 2
    
    xlim = c(0,chr.info$chr.end[n.chr] + gap)
  
  # plot the frame
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, main = title, axes = axes)
  
  # results
  return(chr.info)
}

# -----------------------------------------------------------------------------
# UPDATED Aug 28, 2015 12:26 PM
# FUNCTION:     Manhattan(res.gwas.ID, res.gwas.value, map.sorted,
#                         gap  = 50000000, ylim = NULL, xlab = "Chromosome", 
#                         ylab = "", title= "", axes = TRUE, ltype = "p", 
#                         pch = 16, cols = c("dark blue","cornflowerblue"))
#' @title       Draw Manhattan plots
#' @param       res.gwas.ID     The ID of markers from GWAS result
#' @param       res.gwas.value  The value of markers from GWAS result
#' @param       map.sorted      A sorted map
#' @param       gap             The gap between chromosomes in the plot
#' @param       ylim            The y-axis range
#' @param       xlab,ylab       Parameters passed to plot.default() 
#' @param       title           Character, the title of the plot
#' @param       axes            Logical, draw the axes or not
#' @param       ltype           "lty" passed to plot.default()
#' @param       cols            A Vector with 2 elements of color. can be used
#'                              to distinguish two adjacent chromosomes.
#' @param       pch             Type of points
#' @return      A data frame of chromosome information (from DrawAFrame())
# -----------------------------------------------------------------------------
#' @export      
#' @note        This is a function that draws Manhattan plots
# -----------------------------------------------------------------------------

Manhattan = function(res.gwas.ID,
                     res.gwas.value,
                     map.sorted,
                     gap  = 50000000, ylim = NULL,              # frame
                     xlab = "Chromosome", ylab = "", title= "", # frame
                     axes = TRUE,                               # frame
                     ltype = "p", pch = 20,                     # points
                     cols = c("dark blue","cornflowerblue")     # points
                    ){
  # calculate ylim if needed
    if(is.null(ylim)){  
      ylim = c(0,max(res.gwas.value)*1.1)
    }
  # draw the frame
    res.frame = DrawAFrame(map.sorted, gap = gap, ylim = ylim, axes = FALSE,
                           xlab = xlab, ylab = ylab, title = title)
  
  # result data frame
    res.df = data.frame(ID = as.character(res.gwas.ID), 
                        value = res.gwas.value,
                        stringsAsFactors = FALSE)
    
  # loop across chromosomes 
    for (i.chr in 1:nrow(res.frame)){    
      name.chr = res.frame$chr.id[i.chr]
      point.data = map.sorted$SNPs[[name.chr]]
      point.data = suppressMessages(left_join(point.data,res.df))
      point.data$value[is.na(point.data$value)] = ylim[1] - 1  # hide 
      point.data$x = res.frame$chr.start[i.chr] + point.data$position
      points(point.data$x,
             point.data$value,
             col = cols[(i.chr %% 2) + 1],
             pch = pch, type = ltype)
    }
  # axes
    if (axes){
      axis(side = 2)                                         ### left axis
      axis(side = 1, 
           at     = res.frame$chr.mid,
           labels = res.frame$chr.id)                        ### bottom axis
    }

    return(res.frame)
}

