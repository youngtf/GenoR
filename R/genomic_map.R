# -----------------------------------------------------------------------------
# R file             --Tianfu Yang
# Type:              Functions/Analysis/Visualization
# Subtype/Project:   plot module
# Descriptions:    
# Last Update:       Aug 25, 2015  
# -----------------------------------------------------------------------------
# To-do
# 1 A function that help to group SNPs into windows
# 2 reconsider: why do i need a exported sortmap function?
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Aug 27, 2015 2:59 PM
# geno.map object
# This is a class for genomic maps of SNPs. It is essentially a list: 
# List of 2
# $ Chromosomes:'data.frame':	5 obs. of  3 variables:
#   $ chr.id         : chr [1:5] "0" "1" "2" "X" ...
#   $ chr.is.allosome: logi [1:5] FALSE FALSE FALSE TRUE TRUE
#   $ chr.length     : int [1:5] 71 60 63 97 88
# $ SNPs       :List of 5
#   $ 0:'data.frame':	3 obs. of  3 variables:
#      $ SID       : chr [1:3] "SNP_1" "SNP_3" "SNP_2"
#      $ chromosome: chr [1:3] "0" "0" "0"
#      $ position  : int [1:3] 31 45 71
#   $ 1:'data.frame':	5 obs. of  3 variables:
#      $ SID       : chr [1:5] "SNP_7" "SNP_8" "SNP_5" "SNP_6" ...
#      $ chromosome: chr [1:5] "1" "1" "1" "1" ...
#      $ position  : int [1:5] 19 23 34 41 60
#   $ 2:'data.frame':	3 obs. of  3 variables:
#      $ SID       : chr [1:3] "SNP_11" "SNP_10" "SNP_9"
#      $ chromosome: chr [1:3] "2" "2" "2"
#      $ position  : int [1:3] 2 43 63
#   $ X:'data.frame':	6 obs. of  3 variables:
#      $ SID       : chr [1:6] "SNP_13" "SNP_17" "SNP_15" "SNP_12" ...
#      $ chromosome: chr [1:6] "X" "X" "X" "X" ...
#      $ position  : int [1:6] 6 15 39 79 95 97
#   $ Y:'data.frame':	3 obs. of  3 variables:
#      $ SID       : chr [1:3] "SNP_18" "SNP_19" "SNP_20"
#      $ chromosome: chr [1:3] "Y" "Y" "Y"
#      $ position  : int [1:3] 32 56 88
# - attr(*, "class")= chr "geno.map"
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# UPDATED Aug 27, 2015 4:18 PM
# FUNCTION:     CreateWindowMap(SID, chromosome, position, 
#                               allosome = NULL, distance = 10^6)
#' @title       Create a windows map
#' @param       SID        A vector of SNP IDs (character)
#' @param       chromosome A vector of chromosome (character) where the SNP 
#'                         is located. Support [1-9]{1,2}|[A-Z] ie. 0-99 or 
#'                         a single upper letter. 0 for unmapped SNPs
#' @param       position   A vector of position (integer) where the SNP 
#'                         is located.
#' @param       allosome   A vector of the names of allosomes. By default all
#'                         chromosomes natated with upper letters will be 
#'                         treated as allosomes
#' @param       distance   A integer of window size
#' @return      A geno.map.window object
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function create a windows map for "windows result" based 
#'              original map.
#' @examples    
#' SID = paste("SNP",1:20,sep="_")
#' chromosome = c(rep("0",3),rep("1",5),rep("2",3),rep("X",6),rep("Y",3))
#' position = sample(1:100,20)
#' 
#' data.frame(SID,chromosome,position)
#' #       SID chromosome position
#' # 1   SNP_1          0       31
#' # 2   SNP_2          0       71
#' # 3   SNP_3          0       45
#' # 4   SNP_4          1       60
#' # 5   SNP_5          1       34
#' # 6   SNP_6          1       41
#' # 7   SNP_7          1       19
#' # 8   SNP_8          1       23
#' # 9   SNP_9          2       63
#' # 10 SNP_10          2       43
#' # 11 SNP_11          2        2
#' # 12 SNP_12          X       79
#' # 13 SNP_13          X        6
#' # 14 SNP_14          X       95
#' # 15 SNP_15          X       39
#' # 16 SNP_16          X       97
#' # 17 SNP_17          X       15
#' # 18 SNP_18          Y       32
#' # 19 SNP_19          Y       56
#' # 20 SNP_20          Y       88
#' 
#' CreateWindowMap(SID, chromosome, position, distance = 30)
#' # $Chromosomes
#' #   chr.id chr.is.allosome chr.length chr.n.window
#' # 0      0           FALSE         71            3
#' # 1      1           FALSE         60            2
#' # 2      2           FALSE         63            3
#' # X      X            TRUE         97            4
#' # Y      Y            TRUE         88            3
#' # 
#' # $Windows
#' # $Windows$`0`
#' #    ID chromosome position
#' # 1 0_1          0        0
#' # 2 0_2          0        1
#' # 3 0_3          0        2
#' # 
#' # $Windows$`1`
#' #    ID chromosome position
#' # 1 1_1          1        0
#' # 2 1_2          1        1
#' # 
#' # $Windows$`2`
#' #    ID chromosome position
#' # 1 2_1          2        0
#' # 2 2_2          2        1
#' # 3 2_3          2        2
#' # 
#' # $Windows$X
#' #    ID chromosome position
#' # 1 X_1          X        0
#' # 2 X_2          X        1
#' # 3 X_3          X        2
#' # 4 X_4          X        3
#' # 
#' # $Windows$Y
#' #    ID chromosome position
#' # 1 Y_1          Y        0
#' # 2 Y_2          Y        1
#' # 3 Y_3          Y        2
#' # 
#' # 
#' # attr(,"class")
#' # [1] "geno.map.window"
# -----------------------------------------------------------------------------

CreateWindowMap = function(SID, chromosome, position, 
                           allosome = NULL, 
                           distance = 10^6   ## size of windows
                           ){
  
  ## Reorder the map
    map.sorted = SortMap(SID, chromosome, position, allosome = NULL)
  ## chromosome info
    chr.new = map.sorted$Chromosomes
    chr.new$chr.n.window = ceiling(chr.new$chr.length / distance)
    n.chromosome = length(chr.new$chr.id)
  ## windows info
    map.info = vector("list",n.chromosome)
    names(map.info) = chr.new$chr.id
    for (i in 1:n.chromosome){
      chr.name = chr.new$chr.id[i]
      chr.n.window = chr.new$chr.n.window[i]
      chr.win  = paste0(chr.name,"_", 1:chr.n.window)
      chr.pos  = (1:chr.n.window) - 1
      map.temp = data.frame(ID         = chr.win,
                            chromosome = chr.name,
                            position   = chr.pos,
                            stringsAsFactors = FALSE)
      map.info[[i]] = map.temp[order(chr.pos),]
    }
  ## result 
    map.win = list(Chromosomes = chr.new,
                   Windows     = map.info)
    attributes(map.win)$class = c("geno.map.window",attributes(map.win)$class)
  return(map.win)
}


# -----------------------------------------------------------------------------
# Block
# Descriptions:      
# -----------------------------------------------------------------------------
#' 
#' # SNP2Window(SID, chromosome, position, distance = 30)
#' # $Chromosomes
#' # chr.id chr.is.allosome chr.length chr.n.window
#' # 0      0           FALSE         77            3
#' # 1      1           FALSE         98            4
#' # 2      2           FALSE         88            3
#' # X      X            TRUE         94            4
#' # Y      Y            TRUE         68            3
#' # 
#' # $SNPs
#' # $SNPs$`0`
#' # ID chromosome position window
#' # 2 SNP_2          0       41    0_2
#' # 3 SNP_3          0       72    0_3
#' # 1 SNP_1          0       77    0_3
#' # 
#' # $SNPs$`1`
#' # ID chromosome position window
#' # 3 SNP_6          1       18    1_1
#' # 4 SNP_7          1       22    1_1
#' # 5 SNP_8          1       23    1_1
#' # 1 SNP_4          1       51    1_2
#' # 2 SNP_5          1       98    1_4
#' # 
#' # $SNPs$`2`
#' # ID chromosome position window
#' # 1  SNP_9          2       27    2_1
#' # 2 SNP_10          2       52    2_2
#' # 3 SNP_11          2       88    2_3
#' # 
#' # $SNPs$X
#' # ID chromosome position window
#' # 6 SNP_17          X       11    X_1
#' # 5 SNP_16          X       34    X_2
#' # 3 SNP_14          X       37    X_2
#' # 2 SNP_13          X       60    X_2
#' # 1 SNP_12          X       65    X_3
#' # 4 SNP_15          X       94    X_4
#' # 
#' # $SNPs$Y
#' # ID chromosome position window
#' # 1 SNP_18          Y        2    Y_1
#' # 3 SNP_20          Y       44    Y_2
#' # 2 SNP_19          Y       68    Y_3
#' # 
#' # attr(,"class")
#' # [1] "geno.snp2window" "geno.map"  
# -----------------------------------------------------------------------------
# SNP2Window = function(SID, chromosome, position, 
#                       allosome = NULL, 
#                       distance = 10^6    ## size of windows
#                       ){
#   
#   ## Reorder the map
#   map.sorted = SortMap(SID, chromosome, position, allosome = NULL)
#   ## chromosome info
#   map.sorted$Chromosomes$chr.n.window = ceiling(map.sorted$Chromosomes$chr.length / distance)
#   n.chromosome = length(map.sorted$Chromosomes$chr.id)
#   ## add windows info
#   for (i in 1:n.chromosome){
#     window.idx = ceiling(map.sorted$SNPs[[i]]$position / distance)
#     map.sorted$SNPs[[i]]$window = paste0(map.sorted$SNPs[[i]]$chromosome, "_", window.idx)
#   }
#   ## result 
#   attributes(map.sorted)$class = c("geno.snp2window",attributes(map.sorted)$class)
#   return(map.sorted)
# }
