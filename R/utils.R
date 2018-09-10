#------------------------------------------------------------------------------
# R file  --Tianfu Yang
# Type:              Functions
# Subtype/Project:   Workflow
# Descriptions:      Functions about the workflow and efficiency of analysis
# Last Update:       Aug 26, 2015 2:59 PM
# Contents:
# 1
# FUNCTION:     clc.functions()
# @title        Clean all (customed) functions in present environment
# 2
# FUNCTION:     repmat(a,n,m)
# @title        Replication of matrix (mat) n x m times 
# 3
# FUNCTION:     strsplit.mat(vec.char,sep=" ")
# @title        A matrix version of strsplit()
# 4
# FUNCTION:     blank.remover(vec.char)
# @title        Remove blanks in the end of characters
# 5
# FUNCTION:     ScaleMatrics(mat,scale.coef,MARGIN)
# @title        Matrics scaling
#------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# UPDATED Aug 26, 2015 2:47 PM
# FUNCTION:     clc.functions()
#' @title       clean all (customed) functions in present environment
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function looks for all functions defined in current 
#'              session and cleaning them. It may be useful when a 
#'              environment image is loaded and functions may be masked.
# -----------------------------------------------------------------------------

clc.functions=function(){
    ls.objtype = function(){
        ls.obj = ls(envir = globalenv())
        ls.type = ls.obj
        for (i in 1:length(ls.obj)) ls.type[i] = class(get(ls.obj[i]))
        cbind(ls.obj,ls.type)
    }
    allnames  = ls.objtype()
    funcnames = allnames[which(allnames[,2] == "function"),1]
    rm(list=funcnames,envir=globalenv())
}

# -----------------------------------------------------------------------------
# UPDATED Aug 26, 2015 2:50 PM
# FUNCTION:     repmat(a,n,m)
#' @title       Replication of matrix (mat) n x m times 
#' @param       mat A matrix
#' @param       nrow Number of rows
#' @param       ncol Number of cols
#' @return      A matrix
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function is an application of kronecker multiplication.
#' @examples    
#' mat1 = matrix(1:4,2,2)
#' repmat(mat1,2,3)
#' #      [,1] [,2] [,3] [,4] [,5] [,6]
#' # [1,]    1    3    1    3    1    3
#' # [2,]    2    4    2    4    2    4
#' # [3,]    1    3    1    3    1    3
#' # [4,]    2    4    2    4    2    4
# -----------------------------------------------------------------------------

repmat = function(mat,nrow,ncol){
  kronecker(matrix(1,nrow,ncol),mat)
}

# -----------------------------------------------------------------------------
# UPDATED Aug 26, 2015 2:54 PM
# FUNCTION:     strsplit.mat(vec.char,sep=" ")
#' @title       A matrix version of strsplit()
#' @param       vec.char A vector of characters which will be split.
#' @param       sep      The separator passed to strsplit()
#' @return      A matrix
# -----------------------------------------------------------------------------
#' @export      
#' @note        Instead of outputing a list, this func try to return a matrix.
#' @examples    
#' vec.char = c("QWERTYUIOP",
#'              "ASDFGHJKL:",
#'              "ZXCVBNM<>?")
#' strsplit.mat(vec.char, "")
#' 
#' #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#' # [1,] "Q"  "W"  "E"  "R"  "T"  "Y"  "U"  "I"  "O"  "P"  
#' # [2,] "A"  "S"  "D"  "F"  "G"  "H"  "J"  "K"  "L"  ":"  
#' # [3,] "Z"  "X"  "C"  "V"  "B"  "N"  "M"  "<"  ">"  "?"  
# -----------------------------------------------------------------------------

strsplit.mat = function(vec.char,sep=" "){
    res.list = strsplit(vec.char,sep)
    len.char = unlist(lapply(res.list,length))
    if (length(unique(len.char)) > 1) {
      stop("not the same number of columns")
    }
    res.vec.char = unlist(res.list)
    matrix(res.vec.char,ncol=unique(len.char), byrow=TRUE)
}

# -----------------------------------------------------------------------------
# UPDATED Aug 26, 2015 2:59 PM
# FUNCTION:     blank.remover(vec.char)
# @title       Remove blanks in the end of characters
# @param       vec.char A vector of characters for testing
# @return      A vector of characters
# -----------------------------------------------------------------------------
#  /@export     ## replaced by built-in function R/base::trimws()
# @note        This function detects blanks in the end of characters and 
#              delete them
# @examples    
# vec.char = c("QWERTYUIOP   ",
#              "ASDFGHJKL:  ",
#              "ZXCVBNM<>? ")
# blank.remover(vec.char)
# # Blank found! 
# # Blank found! 
# # Blank found! 
# # [1] "QWERTYUIOP" "ASDFGHJKL:" "ZXCVBNM<>?"
# 
# -----------------------------------------------------------------------------
# 
# blank.remover = function(vec.char){
#     blank.rm.ele = function(char){
#         blank_found = FALSE
#         n = nchar(char)
#         while(substr(char,n,n) == " "){
#             blank_found = TRUE
#             char = substr(char,1,(n-1))
#             n = nchar(char)
#         }
#         if (blank_found) cat("Blank found! \n")
#         char
#     }
#     res.char = vec.char
#     for (i in 1:length(vec.char)){
#       if (!is.na(res.char[i]))
#         res.char[i] = blank.rm.ele(vec.char[i])
#     } 
#     res.char
# }

# -----------------------------------------------------------------------------
# UPDATED Aug 26, 2015 3:22 PM
# FUNCTION:     ScaleMatrics(mat,scale.coef,MARGIN)
#' @title       Matrics scaling
#' @param       mat        A numeric matrix (m * n)
#' @param       scale.coef A numeric factor (m or n)
#' @param       MARGIN     An integer:
#'                         1 for row manipulation and 2 for column manipulation
#' @return      An sparse numeric matrix m * n
# -----------------------------------------------------------------------------
#' @export      
#' @note        Each vector in the matrix is scaled using the corresponding
#'              element in the scale.coef.
#'              Left multiplication manipulate the rows, while right 
#'              multiplication manipulate the cols. This function uses diagonal 
#'              sparse functions to multiply each row/column by a nonzero 
#'              scalar (can be different across rows/columns)
#'              Thanks for Dr. Zhiquan Wang's contribution.
#' @examples 
#' (mat = matrix(1,3,5))
#' # [,1] [,2] [,3] [,4] [,5]
#' # [1,]    1    1    1    1    1
#' # [2,]    1    1    1    1    1
#' # [3,]    1    1    1    1    1
#' 
#' ScaleMatrics(mat,c(2,3,5),1)
#' # 3 x 5 Matrix of class "dgeMatrix"
#' #      [,1] [,2] [,3] [,4] [,5]
#' # [1,]    2    2    2    2    2
#' # [2,]    3    3    3    3    3
#' # [3,]    5    5    5    5    5
#' 
#' ScaleMatrics(mat,c(2,3,5,7,8),2)
#' # 3 x 3 Matrix of class "dgeMatrix"
#' # [,1] [,2] [,3]
#' # [1,]    2    3    5
#' # [2,]    2    3    5
#' # [3,]    2    3    5
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# May 8, 2015 10:09 AM
# 1 May not be able to handle missing value 
# 2 variable m is not correctly assign when MARGIN = 2        check!
# -----------------------------------------------------------------------------

ScaleMatrics = function(mat,scale.coef,MARGIN){
  n.coef = length(scale.coef)
  ## a sparse diagnal matrix
  scale.matrix = sparseMatrix(i=seq(n.coef),j=seq(n.coef),x=scale.coef)
  if (MARGIN == 1){
    m = nrow(mat)
    if (m != n.coef) {
      stop("Wrong dimension of input matrices")
    }
    res = scale.matrix %*% mat
  } else if (MARGIN == 2){
    m = ncol(mat)
    if (m != n.coef) {
      stop("Wrong dimension of input matrices")
    }
    res = mat %*% scale.matrix
  } else {
    stop("Should not get here!")
  }
  return(res)
}
