#------------------------------------------------------------------------------
# R file  --Tianfu Yang
# Type:              Functions
# Subtype/Project:   Workflow
# Descriptions:      Functions about the workflow and efficiency of analysis
# Last Update:       2014-05-19
# Contents:
# 1
#   Function:    sv(data,rows=1:3)
#   Description: An alternative of head(), with column numbers
# 2
#   Function:    clc.func()
#   Description: clean all functions in present environment
# 3
#   Function:    repmat(a,n,m)
#   Description: replication of matrix (a) n x m times 
# 4
#   Function:    strsplit.mat(vec.char,sep=" ")
#   Description: instead of a list, this func try to return a matrix.
# 5
#   Function:    blank.remover(vec.char)
#   Description: to remove the blank in the end of characters
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Function:    sv(data,rows=1:3)
# Description: An alternative of head(), with column numbers
#------------------------------------------------------------------------------
sv = function(data, rows=1:3){
    ncol.d = ncol(data)
    dataToPrint = rbind(colnames(data),as.matrix(data[rows,]))
    rownames(dataToPrint)[1] = "COLNAME"
    colnames(dataToPrint) = 1:ncol.d
    print(dataToPrint,quote=FALSE)
}

#------------------------------------------------------------------------------
# Function:    clc.func()
# Description: clean all functions in present environment
#------------------------------------------------------------------------------
clc.func=function(){
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

#------------------------------------------------------------------------------
# Function:    repmat(a,n,m)
# Description: replication of matrix (a) n x m times 
#------------------------------------------------------------------------------
repmat <- function(a,nrow,ncol) {kronecker(matrix(1,nrow,ncol),a)}

#------------------------------------------------------------------------------
# Function:    strsplit.mat(vec.char,sep=" ")
# Description: instead of a list, this func try to return a matrix.
#------------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Unit test 
# vec.char = c("QWERTYUIOP",
#              "ASDFGHJKL:",
#              "ZXCVBNM<>?")
# strsplit.mat(vec.char, "")
# 
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,] "Q"  "W"  "E"  "R"  "T"  "Y"  "U"  "I"  "O"  "P"  
# [2,] "A"  "S"  "D"  "F"  "G"  "H"  "J"  "K"  "L"  ":"  
# [3,] "Z"  "X"  "C"  "V"  "B"  "N"  "M"  "<"  ">"  "?"  

# -----------------------------------------------------------------------------

strsplit.mat = function(vec.char,sep=" "){
    res.list = strsplit(vec.char,sep)
    len.char = unlist(lapply(res.list,length))
    if (length(unique(len.char)) > 1) stop("not the same number of columns")
    res.vec.char = unlist(res.list)
    matrix(res.vec.char,ncol=unique(len.char), byrow=TRUE)
}

#------------------------------------------------------------------------------
# Function:    blank.remover(vec.char)
# Description: to remove the blank in the end of characters
#------------------------------------------------------------------------------
blank.remover = function(vec.char){
    blank.rm.ele = function(char){
        n = nchar(char)
        while(substr(char,n,n) == " "){
            char = substr(char,1,(n-1))
            n = nchar(char)
        }
        char
    }
    res.char = vec.char
    for (i in 1:length(vec.char)){
        res.char[i] = blank.rm.ele(vec.char[i])
    } 
    res.char
}

# -----------------------------------------------------------------------------
# Updated May 8, 2015 9:58 AM
# Function:     rowScaling(InputMatrix,scaleFactor)
#               colScaling(InputMatrix,scaleFactor)
# Description:  Each vector in the matrix is scaled using the corresponding
#               element in the scaleFactor
# input:        InputMatrix: a numeric matrix (m * n)
#               scaleFactor: a numeric factor (m / n)
# ouput:        a numeric matrix m * n
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# May 8, 2015 10:09 AM
# 1 May not be able to handle missing value 
# -----------------------------------------------------------------------------

rowScaling = function(InputMatrix,scaleFactor){
  n1 = nrow(InputMatrix)
  n2 = length(scaleFactor)
  if (n1 != n2) {
    stop("Wrong dimension of input matrices")
  }
  scaleMatrix = sparseMatrix(i=seq(n2),j=seq(n2),x=scaleFactor)
  res = scaleMatrix %*% InputMatrix
  res
}

colScaling = function(InputMatrix,scaleFactor){
  n1 = ncol(InputMatrix)
  n2 = length(scaleFactor)
  if (n1 != n2) {
    stop("Wrong dimension of input matrices")
  }
  scaleMatrix = sparseMatrix(i=seq(n2),j=seq(n2),x=scaleFactor)
  res = InputMatrix %*% scaleMatrix
  res
}