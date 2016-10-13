#------------------------------------------------------------------------------
# R file  --Tianfu Yang
# Type:              Functions
# Subtype/Project:   Ploting
# Descripetions:     functions used in ploting
# Last Update:       Aug 27, 2015 5:10 PM
# Contents:

# -----------------------------------------------------------------------------
# Note
# The basic logic in this module:
# 1. the mapinfo is pre-processed for further use
# 2. based on a given "plotting map" (the genomic regions of interests), a
#    genetic position can be converted into a x-coordinate
# -----------------------------------------------------------------------------
# Note
# Standardization of Chr ID: Natural sort
# -----------------------------------------------------------------------------
# Note
# The x-axis can be treated as a set of bins, where a bin can be a chromosome, 
# or a specific genomic region on a chromosome. map2x() maps a genomic 
# map and a "bin list" into a simple list of x-coordinate, 
# -----------------------------------------------------------------------------
# sample bin_list: three column data.frame
# Chr   start_pos end_pos
# [str] [int]     [int]
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Block
# Descriptions:      new sortmap
# -----------------------------------------------------------------------------
load("~/Dropbox/P_GWAS/3_beefMQ/DarkCut/F2_TidyData/map80Raw.Rdata")

sorted_map = sort_map(MarID = Map80k$MarID, 
                      chromosome = Map80k$Chr,
                      position = Map80k$Pos,
                      allosome = c("X","Y"),
                      unmapped = "0",
                      mitochondria = "MT"
                      )
sorted_map_ext = subset_sorted_map(sorted_map)
bin_info_df    = create_plot_map(sorted_map_ext, 50000000)

sort_map = function(MarID, 
                    chromosome, 
                    position, 
                    allosome = NULL,     ## default: no allosome
                    unmapped = NULL,     ## default: no unmapped
                    mitochondria = NULL  ## default: no MT genome
                   ){
  
  chr_factor = naturalfactor(chromosome)
  res_map = as.tbl(data.frame(MarID            = MarID, 
                              Chr              = chr_factor, 
                              Pos              = position, 
                              stringsAsFactors = F))
  
  ## add attributes
  attributes(res_map)$chrom_id     = levels(chr_factor)
  attributes(res_map)$n_chromosome = length(attributes(res_map)$chrom_id)
  attributes(res_map)$allosome     = allosome
  attributes(res_map)$unmapped     = unmapped
  attributes(res_map)$MT           = mitochondria
  attributes(res_map)$chr_len      = tapply(res_map$Pos, 
                                            res_map$Chr, max)
  
  res_map_sorted = res_map[order(res_map$Chr, res_map$Pos),]

  return(res_map)
}

subset_sorted_map = function(sorted_map, 
                             bin_list = NULL, 
                             plot_unmapped = F,
                             plot_allosome = F,
                             plot_mt = F){
  if_all_chr = is.null(bin_list)
  ## bin_list pre-process
  if(is.null(bin_list)){
    bin_list = data.frame(Chr   = attributes(sorted_map)$chrom_id,
                          Start = 0,
                          End   = attributes(sorted_map)$chr_len)
    sorted_map_ext = mutate(sorted_map, bin = as.numeric(Chr))
  } else {
    sorted_map_ext = mutate(sorted_map, bin = NA)
    for (i in 1:nrow(bin_list)){
      if_selected = (as.character(sorted_map_ext$Chr) == as.character(bin_list$Chr[i])) &
                    (sorted_map_ext$Pos <= bin_list$End[i]) &
                    (sorted_map_ext$Pos >= bin_list$Start[i])
      sorted_map_ext$bin[if_selected] = i
    } 
  }
  ## attributes
  kept_map_info = attributes(sorted_map)
  kept_map_info$row.names = NULL
  kept_map_info$names = NULL
  kept_map_info$class = NULL
  attributes(sorted_map_ext) = c(attributes(sorted_map_ext), kept_map_info)
  ## remove useless bins
  if(!plot_unmapped & !is.null(attributes(sorted_map_ext)$unmapped)){
    is_unmapped = as.character(sorted_map_ext$Chr) %in% attributes(sorted_map_ext)$unmapped
    sorted_map_ext$bin[is_unmapped] = NA
  }
  if(!plot_allosome & !is.null(attributes(sorted_map_ext)$allosome)){
    is_allosome = as.character(sorted_map_ext$Chr) %in% attributes(sorted_map_ext)$allosome
    sorted_map_ext$bin[is_allosome] = NA
  }
  if(!plot_mt & !is.null(attributes(sorted_map_ext)$MT)){
    is_mt = as.character(sorted_map_ext$Chr) %in% attributes(sorted_map_ext)$MT
    sorted_map_ext$bin[is_mt] = NA
  }
  attributes(sorted_map_ext)$if_all_chr = if_all_chr
  return(sorted_map_ext)
}

create_plot_map = function(sorted_map_ext, gap){
  bin_info = by(sorted_map_ext, sorted_map_ext$bin, 
                function(x) data.frame(chr = unique(as.character(x$Chr)), 
                                       start = min(x$Pos), 
                                       end = max(x$Pos)), 
                simplify = T)
  bin_info_df = do.call("rbind", bin_info)
  bin_info_df$bin = as.integer(names(bin_info))
  if (attributes(sorted_map_ext)$if_all_chr){
    bin_info_df$start = 1
  }
  bin_info_df$chr_length = bin_info_df$end - bin_info_df$start + 1
  # new map
  
  # |========|----|======|----|....|========|
  #   bin1    int1  bin 2 int2      last_bin
  #
  n_chr   = nrow(bin_info_df)
  gap_dev = gap * (0:(n_chr-1)) ## position deviation because of gap
  cumu_pos = cumsum(c(1,(bin_info_df$chr_length)))[1:(n_chr)]
  
  bin_info_df$x_chr_start =  cumu_pos + gap_dev
    
  bin_info_df$x_chr_end   = bin_info_df$x_chr_start + bin_info_df$chr_length - 1
  bin_info_df$x_chr_mid   = (bin_info_df$x_chr_start + bin_info_df$x_chr_end) / 2
  
  return(bin_info_df)
}

# -----------------------------------------------------------------------------
# UPDATED Aug 27, 2015 5:35 PM
# FUNCTION:     DrawAFrame(map.sorted, gap, 
#                          ylim, axes = FALSE,
#'                         xlab = "", ylab = "", title = "")
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
#' # DrawAFrame(map.sorted, gap = 30, ylim = c(1,3))        ## not exported
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


#   
# # -----------------------------------------------------------------------------
# # Aug 27, 2015 2:59 PM
# # geno.map object
# # This is a class for genomic maps of SNPs. It is essentially a list:
# # List of 2
# # $ Chromosomes:'data.frame':	5 obs. of  3 variables:
# #   $ chr.id         : chr [1:5] "0" "1" "2" "X" ...
# #   $ chr.is.allosome: logi [1:5] FALSE FALSE FALSE TRUE TRUE
# #   $ chr.length     : int [1:5] 71 60 63 97 88
# # $ SNPs       :List of 5
# #   $ 0:'data.frame':	3 obs. of  3 variables:
# #      $ SID       : chr [1:3] "SNP_1" "SNP_3" "SNP_2"
# #      $ chromosome: chr [1:3] "0" "0" "0"
# #      $ position  : int [1:3] 31 45 71
# #   $ 1:'data.frame':	5 obs. of  3 variables:
# #      $ SID       : chr [1:5] "SNP_7" "SNP_8" "SNP_5" "SNP_6" ...
# #      $ chromosome: chr [1:5] "1" "1" "1" "1" ...
# #      $ position  : int [1:5] 19 23 34 41 60
# #   $ 2:'data.frame':	3 obs. of  3 variables:
# #      $ SID       : chr [1:3] "SNP_11" "SNP_10" "SNP_9"
# #      $ chromosome: chr [1:3] "2" "2" "2"
# #      $ position  : int [1:3] 2 43 63
# #   $ X:'data.frame':	6 obs. of  3 variables:
# #      $ SID       : chr [1:6] "SNP_13" "SNP_17" "SNP_15" "SNP_12" ...
# #      $ chromosome: chr [1:6] "X" "X" "X" "X" ...
# #      $ position  : int [1:6] 6 15 39 79 95 97
# #   $ Y:'data.frame':	3 obs. of  3 variables:
# #      $ SID       : chr [1:3] "SNP_18" "SNP_19" "SNP_20"
# #      $ chromosome: chr [1:3] "Y" "Y" "Y"
# #      $ position  : int [1:3] 32 56 88
# # - attr(*, "class")= chr "geno.map"
# # -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# UPDATED Aug 27, 2015 1:11 PM
# FUNCTION:     SortMap(SID,chromosome,position)
# @title       build a geno_map object with sorted information
# @param       SID        A vector of SNP IDs (character)
# @param       chromosome A vector of chromosome (character) where the SNP 
#                         is located. Support [1-9]{1,2}|[A-Z] ie. 0-99 or 
#                         a single upper letter. 0 for unmapped SNPs
# @param       position   A vector of position (integer) where the SNP 
#                         is located.
# @param       allosome   A vector of the names of allosomes. By default all
#                         chromosomes natated with upper letters will be 
#                         treated as allosomes
# @return      A geno.map object
# -----------------------------------------------------------------------------
# @export      
# @note        This is a function for sortting map information of SNPs.
# @examples    
# SID = paste("SNP",1:20,sep="_")
# chromosome = c(rep("0",3),rep("1",5),rep("2",3),rep("X",6),rep("Y",3))
# position = c(31, 71, 45, 60, 34, 41, 19, 23, 63, 43,  
#               2, 79,  6 ,95, 39, 97, 15, 32, 56, 88)
# 
# data.frame(SID,chromosome,position)
# #       SID chromosome position
# # 1   SNP_1          0       31
# # 2   SNP_2          0       71
# # 3   SNP_3          0       45
# # 4   SNP_4          1       60
# # 5   SNP_5          1       34
# # 6   SNP_6          1       41
# # 7   SNP_7          1       19
# # 8   SNP_8          1       23
# # 9   SNP_9          2       63
# # 10 SNP_10          2       43
# # 11 SNP_11          2        2
# # 12 SNP_12          X       79
# # 13 SNP_13          X        6
# # 14 SNP_14          X       95
# # 15 SNP_15          X       39
# # 16 SNP_16          X       97
# # 17 SNP_17          X       15
# # 18 SNP_18          Y       32
# # 19 SNP_19          Y       56
# # 20 SNP_20          Y       88
# 
# SortMap(SID, chromosome, position)
# # $Chromosomes
# #   chr.id chr.is.allosome chr.length
# # 0      0           FALSE         71
# # 1      1           FALSE         60
# # 2      2           FALSE         63
# # X      X            TRUE         97
# # Y      Y            TRUE         88
# # 
# # $SNPs
# # $SNPs$`0`
# #    ID    chromosome position
# # 1 SNP_1           0       31
# # 3 SNP_3           0       45
# # 2 SNP_2           0       71
# # 
# # $SNPs$`1`
# #    ID    chromosome position
# # 4 SNP_7           1       19
# # 5 SNP_8           1       23
# # 2 SNP_5           1       34
# # 3 SNP_6           1       41
# # 1 SNP_4           1       60
# # 
# # $SNPs$`2`
# #    ID    chromosome position
# # 3 SNP_11          2        2
# # 2 SNP_10          2       43
# # 1  SNP_9          2       63
# # 
# # $SNPs$X
# #    ID    chromosome position
# # 2 SNP_13          X        6
# # 6 SNP_17          X       15
# # 4 SNP_15          X       39
# # 1 SNP_12          X       79
# # 3 SNP_14          X       95
# # 5 SNP_16          X       97
# # 
# # $SNPs$Y
# #    ID    chromosome position
# # 1 SNP_18          Y       32
# # 2 SNP_19          Y       56
# # 3 SNP_20          Y       88
# # 
# # 
# # attr(,"class")
# # [1] "geno.map"
# -----------------------------------------------------------------------------
# SortMap_0 = function(SID, chromosome, position, allosome = NULL){
#   # checking chromosome information:
#   chromname = as.character(unique(chromosome))
#   n.chromosome = length(chromname)
#   idx.num      = grep("[1-9]{1,2}",chromname)
#   idx.char     = grep("[A-Z]",     chromname)
#   idx.unmapped = which(chromname == "0")
#   if (any(!(seq(n.chromosome) %in% c(idx.num,idx.char,idx.unmapped)))){
#     stop("invalid chromosome name(s):",
#          chromname[!(seq(n.chromosome) %in% c(idx.num,idx.char,idx.unmapped))]
#     )
#   }
#   # chromosome ID
#   chr.ID = character(0)
#   if (length(idx.unmapped) == 1){
#     chr.ID = c(chr.ID, "0")
#   }
#   if (length(idx.num) > 0){
#     chr.ID = c(chr.ID, as.character(sort(as.integer(chromname[idx.num]))))
#   }
#   if (length(idx.num) > 0){
#     chr.ID = c(chr.ID, sort(chromname[idx.char]))
#   }
#   # isAutosome
#   if (is.null(allosome)){
#     chr.is.allosome = (chr.ID %in% chromname[idx.char])
#   } else {
#     chr.is.allosome = (chr.ID %in% allosome)
#   }
#   
#   # Maps
#   map.info = vector("list",n.chromosome)
#   names(map.info) = chr.ID
#   for (i in 1:n.chromosome){
#     chr.name = chr.ID[i]
#     idx.this.chr = which(chromosome == chr.name)
#     chr.snp = SID[idx.this.chr]
#     chr.pos = position[idx.this.chr]
#     map.temp = data.frame(ID         = chr.snp,
#                           chromosome = chr.name,
#                           position   = chr.pos,
#                           stringsAsFactors = FALSE)
#     map.info[[i]] = map.temp[order(chr.pos),]
#   }
#   # chromosome length
#   chr.length = unlist(lapply(map.info,function(x) max(x$position)))
#   
#   res.map = list(Chromosomes = data.frame(chr.id          = chr.ID,
#                                           chr.is.allosome = chr.is.allosome,
#                                           chr.length      = chr.length,
#                                           stringsAsFactors = FALSE),
#                  SNPs        = map.info)
#   attributes(res.map)$class = c("geno.map",attributes(res.map)$class)
#   return(res.map)
# }
