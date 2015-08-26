#------------------------------------------------------------------------------
# R file  --Tianfu Yang
# Type:              Functions
# Subtype/Project:   GenoCleaning
# Descripetions:     functions used to clean the genotype data
# Last Update:       Aug 26, 2015 1:08 PM
# Contents:
#
#                ______        ______  <5>   ______        __        _________ 
#      RAW ===> | "TG" |      | "AB" | ---> |-1/0/1| ---> |QC| ---> | imputed |
#               |______|      |______| <--- |______| <7>  |__| <8>  |_________|
#                  |^            |^     <6>    ^                         |
#              <1> || <2>    <1> ||<2>         |                         |
#                  v|            v|            | <4>                     v
#  ______   <9>  ______        ______          |                    __________     
# | PED  | <--- |"T""G"| ---> |"A""B"| ---------                   | analysis |
# |______|      |______| <3>  |______|                             |__________|
# 
# 1                                                       
# FUNCTION:     geno.biAL2monoAL()
# @title        divide bi-allelic genotype data into mono-allelic
# 2                                                 
# FUNCTION:     geno.monoAL2biAL(monoAlle,na.geno,na.output)
# @title        combine the allele data into Bi-allele data
# 3                                                 
# FUNCTION:     geno.AGCT2AB(monoAlle,na.output)
# @title        Transform the coding of monoAL data to (A/B/NA).output
# 4                                                 
# FUNCTION:     geno.monoAL2Num(monoAlle_AB,na.geno,na.output)
# @title        Transform the coding from mono-allele to integer (-1/0/1/-5)
# 5
# FUNCTION:     geno.AB2Num(genodata,na.geno = NULL)
# @title        Transform the coding between BB/AB/AA and -1/0/1
# 6
# FUNCTION:     geno.Num2AB(genodata.int,code=c(1,0,-1),na.geno = NULL)
# @title        Transform the coding between BB/AB/AA and -1/0/1
# 7
# FUNCTION:     geno.qc(genodata,code=c(1,0,-1),na.string = NULL)
# @title        Quality control for genotype data
# 8 
# FUNCTION:     geno.impu.ave(genodata)
# @title        Impute the missing cells with average
# 9
# FUNCTION:     geno.mono2PED(monodata,IndID,lite,FamilyID,PID,MID,Sex,Pheno)
# @title        Create PED files
#------------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Notes
# Descriptions:      To-do List
# To-do              1 construction of AB matrix          CHECK!      
#                    2 unit tests with "CC/TT"            CHECK! 
#                    3 comments for PED                   CHECK!
#                    4 re-write QC function               CHECK!
#                    5 combine PedR and PED function
# Last Update:       Aug 26, 2015
# -----------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Feb 16, 2015 14:12
# monoAlle: A class for the mono-allelic genotypic data
# this object is a list 
# monoAlle
#  - dim : a vector with 2 elements (nind and nmar)
#  - names.ind
#  - names.mar
#  - names.base : a 2-row matrix of the bases (A/B or A/G/C/T)
#  - Allele_1   : matrix with code A/B
#  - Allele_2   : matrix with code A/B
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# UPDATED Aug 25, 2015 2:56 PM
# FUNCTION:     geno.biAL2monoAL()
#' @title       Divide bi-allelic genotype data into mono-allelic
#' @param       genodata  A matrix of bi-allelic genotype (n(nind) x m(nmar))
#' @param       na.geno   A 2-digit character for missing genotype in genedata
#' @param       na.output A 2-digit character for out put data 
#' @return      A monoAlle onject
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function divides bi-allelic genotype matrix into a
#'              mono-allelic object
#' @examples    
#' genodata = rbind(c("AA","AB","--","CC")
#'                 ,c("AA","--","BB","TT")
#'                 ,c("AB","AB","AB","--"))
#' (monoAlle = geno.biAL2monoAL(genodata))
#' #
#' # $dim
#' # Row Col 
#' # 3   4 
#' # 
#' # $names.ind
#' # [1] "R1" "R2" "R3"
#' # 
#' # $names.mar
#' # [1] "C1" "C2" "C3" "C4"
#' # 
#' # $names.base
#' # $names.base$C1
#' # [1] "A" "B"
#' # 
#' # $names.base$C2
#' # [1] "A" "B"
#' # 
#' # $names.base$C3
#' # [1] "B" "A"
#' # 
#' # $names.base$C4
#' # [1] "C" "T"
#' # 
#' # 
#' # $Allele_1
#' # C1  C2  C3  C4 
#' # R1 "A" "A" "0" "C"
#' # R2 "A" "0" "B" "T"
#' # R3 "A" "A" "A" "0"
#' # 
#' # $Allele_2
#' # C1  C2  C3  C4 
#' # R1 "A" "B" "0" "C"
#' # R2 "A" "0" "B" "T"
#' # R3 "B" "B" "B" "0"
#' # 
#' # attr(,"class")
#' # [1] "monoAL"
#' # Warning messages:
#' # 1: In geno.biAL2monoAL(genodata) :
#' #   no valid rownames for the input data, simple names will be given.
#' # 2: In geno.biAL2monoAL(genodata) :
#' #   no valid colnames for the input data, simple names will be given.
# -----------------------------------------------------------------------------

geno.biAL2monoAL = function(genodata
                            ,na.geno = "--"   # all genotype should be biallelic
                            ,na.output = "00" # widely used in Plink
                            ){
  
  ## function for allele list
    GenoList2AlleleList = function(genolist){
      allele_1 = substr(genolist,1,1)
      allele_2 = substr(genolist,2,2)
      alleleList = unique(c(allele_1,allele_2))
    }
  ## output res
    res.data = vector("list",6)
    names(res.data) = c("dim",
                        "names.ind","names.mar","names.base",
                        "Allele_1","Allele_2")
  ## result: dim
    res.data$dim = dim(genodata)
    names(res.data$dim) = c("Row","Col")
  ## check the data: dim
    if (res.data$dim[1] >= res.data$dim[2]){
      warning("More rows than columns for the input data.")
    }  
  ## check the data: dimnames
    if (is.null(rownames(genodata))){
      warning("no valid rownames for the input data, simple names will be given.")
      rownames(genodata) = paste0("R",1:res.data$dim[1])
    }
    if (is.null(colnames(genodata))){
      warning("no valid colnames for the input data, simple names will be given.")
      colnames(genodata) = paste0("C",1:res.data$dim[2])
    }
  ## result: names
    res.data$names.ind = rownames(genodata)
    res.data$names.mar = colnames(genodata)
  ## Genotype data
  ## na.string: na.output
    if(any(genodata == na.geno)){
      genodata[genodata == na.geno] = na.output
    }
  ## split the genotype 
    # substr() can be used in matries
    res.data$Allele_1  = substr(genodata,1,1)
    res.data$Allele_2  = substr(genodata,2,2)
  ## allele list
    genotype.list = apply(genodata,2,unique)
    genotype.list.nona = lapply(genotype.list
                                ,function(x) x[x != na.output])
    allele.list = lapply(genotype.list.nona
                        ,GenoList2AlleleList)
    n.allele = unlist(lapply(allele.list,length))
  ## check the allele number (should be biploid)
    if (any(n.allele > 2)){
    idx.innormal = which(n.allele > 2)
    warning("more than 2 possible alleles: ",idx.innormal)
  }
  ## result: allele list
    res.data$names.base = allele.list
  ## Attributes
    attributes(res.data)$class = c("monoAL",attributes(res.data)$class)
  ## Return  
    res.data
}
# -----------------------------------------------------------------------------
# UPDATED Aug 25, 2015 3:23 PM
# FUNCTION:     geno.monoAL2biAL(monoAlle,na.geno,na.output)
#' @title       combine the allele data into Bi-allele data
#' @param       monoAlle  A monoAL object
#' @param       na.geno   A 2-digit character for missing genotype in monoAlle
#' @param       na.output A 2-digit character for out put data 
#' @return      A matrix of bi-allelic genotype
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function combines the allele data into Bi-allele data
#' @examples    
#' genodata = rbind(c("AA","AB","--","CC")
#'                 ,c("AA","--","BB","TT")
#'                 ,c("AB","AB","AB","--"))
#' monoAlle = geno.biAL2monoAL(genodata)
#' geno.monoAL2biAL(monoAlle)
#' #    C1   C2   C3   C4  
#' # R1 "AA" "AB" "--" "CC"
#' # R2 "AA" "--" "BB" "TT"
#' # R3 "AB" "AB" "AB" "--"
# -----------------------------------------------------------------------------

geno.monoAL2biAL = function(monoAlle
                            ,na.geno   = "00" # widely used in Plink
                            ,na.output = "--" # 
                            ){
  ## check the type of biallele data
    if (!("monoAL" %in% attributes(monoAlle)$class)){
      stop("The input data should be a monoAL data")
    }
  ## Result
    res.bia = paste0(monoAlle$Allele_1,monoAlle$Allele_2)
    dim(res.bia) = monoAlle$dim
    rownames(res.bia) = monoAlle$names.ind
    colnames(res.bia) = monoAlle$names.mar
    res.bia[res.bia == na.geno] = na.output
    res.bia
}

# -----------------------------------------------------------------------------
# UPDATED Aug 25, 2015 3:27 PM
# FUNCTION:     geno.AGCT2AB(monoAlle,na.output)
#' @title       Transform the coding of monoAL data to (A/B/NA).output
#' @param       monoAlle  A monoAL object
#' @param       na.output A integer for out put data 
#' @return      A monoAL_AB object
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function transforms the coding of monoAL data to (A/B/NA) 
#'              output, which is the index in the allele list / missing.
#               The A/B code is only determined by the allele list in the 
#               raw monoAL data (which is also kept in the result)
#' @examples    
#' genodata = rbind(c("AA","AB","--","CC")
#'                 ,c("AA","--","BB","TT")
#'                 ,c("AB","AB","AB","--"))
#' monoAlle = geno.biAL2monoAL(genodata)
#' (monoAlle.AB = geno.AGCT2AB(monoAlle))
#' # $dim
#' # Row Col 
#' # 3   4 
#' # 
#' # $names.ind
#' # [1] "R1" "R2" "R3"
#' # 
#' # $names.mar
#' # [1] "C1" "C2" "C3" "C4"
#' # 
#' # $names.base
#' # $names.base$C1
#' # [1] "A" "B"
#' # 
#' # $names.base$C2
#' # [1] "A" "B"
#' # 
#' # $names.base$C3
#' # [1] "B" "A"
#' # 
#' # $names.base$C4
#' # [1] "C" "T"
#' # 
#' # 
#' # $Allele_1
#' # [,1] [,2] [,3] [,4]
#' # [1,] "A"  "A"  "0"  "A" 
#' # [2,] "A"  "0"  "A"  "B" 
#' # [3,] "A"  "A"  "B"  "0" 
#' # 
#' # $Allele_2
#' # [,1] [,2] [,3] [,4]
#' # [1,] "A"  "B"  "0"  "A" 
#' # [2,] "A"  "0"  "A"  "B" 
#' # [3,] "B"  "B"  "A"  "0" 
#' # 
#' # attr(,"class")
#' # [1] "monoAL"    "monoAL_AB"
# -----------------------------------------------------------------------------

geno.AGCT2AB = function(monoAlle,na.output = 0){
  
  ## check the type of biallele data
    if (!("monoAL" %in% attributes(monoAlle)$class)){
      stop("The input data should be a monoAL data")
    }
  ## dim
    nmar = monoAlle$dim[2]
    nind = monoAlle$dim[1]
  ## initial result list     
    res.list = vector("list",0)
    res.list$dim        = monoAlle$dim
    res.list$names.ind  = monoAlle$names.ind
    res.list$names.mar  = monoAlle$names.mar
    res.list$names.base = monoAlle$names.base
  ## getting allele info  
    allele_List = monoAlle$names.base
  ## transfer
    res.list$Allele_1 = apply(t(seq(nmar)), 2
                             ,function(x) 
                              c("A","B")[match(monoAlle$Allele_1[,x]
                                        ,allele_List[[x]])]
                             )
    res.list$Allele_2 = apply(t(seq(nmar)), 2
                             ,function(x) 
                              c("A","B")[match(monoAlle$Allele_2[,x]
                                              ,allele_List[[x]])]
                             )
  ## All NA is missing
    res.list$Allele_1[is.na(res.list$Allele_1)] = na.output
    res.list$Allele_2[is.na(res.list$Allele_2)] = na.output
  ## Attributes
    attributes(res.list)$class = c("monoAL","monoAL_AB",
                                   attributes(res.list)$class)
  ## return the result
    res.list
}

# -----------------------------------------------------------------------------
# UPDATED Aug 26, 2015 1:08 PM
# FUNCTION:     geno.monoAL2Num(monoAlle_AB,na.geno,na.output)
#' @title       transform the coding from mono-allele to integer (-1/0/1/-5)
#' @param       monoAlle_AB A monoAL_AB object
#' @param       na.geno   An integer for missing genotype in monoAlle
#' @param       na.output An integer for out put data 
#' @return      A matrix of integer
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function transforms the coding from mono-allele to 
#'              integer (-1/0/1/-5)
#' @examples    
#' genodata = rbind(c("AA","AB","--","CC")
#'                 ,c("AA","--","BB","TT")
#'                 ,c("AB","AB","AB","--"))
#' monoAlle = geno.biAL2monoAL(genodata)
#' monoAlle.AB = geno.AGCT2AB(monoAlle)
#' (mat.int = geno.monoAL2Num(monoAlle.AB))
#' #    C1 C2 C3 C4
#' # R1  1  0 -5  1
#' # R2  1 -5  1 -1
#' # R3  0  0  0 -5
# -----------------------------------------------------------------------------

geno.monoAL2Num = function(monoAlle_AB
                           ,na.geno   = 0
                           ,na.output = -5
                           ){
  ## check the type of biallele data
    if (!("monoAL_AB" %in% attributes(monoAlle_AB)$class)){
      stop("The input data should be a monoAL_AB data")
    }
    
  ## Mark those NAs 
    Allele_1 = monoAlle_AB$Allele_1
    Allele_1[Allele_1 == na.geno] = NA
    Allele_2 = monoAlle_AB$Allele_2
    Allele_2[Allele_2 == na.geno] = NA
  ## change the coding  
    res.num = (Allele_1 == Allele_2) *        # homo or not
              sign((Allele_1 == 1) - 0.5 )    # first AL or not
  ## other info
    dim(res.num) = monoAlle_AB$dim
    rownames(res.num) = monoAlle_AB$names.ind
    colnames(res.num) = monoAlle_AB$names.mar
    res.num[is.na(res.num)] = na.output
    res.num
}

# -----------------------------------------------------------------------------
# UPDATED Aug 25, 2015 3:58 PM
# FUNCTION:     geno.AB2Num(genodata, na.geno = NULL)
#' @rdname      geno.AB2Num
#' @title       transform the coding between BB/AB/AA and -1/0/1
#' @param       genodata  A matrix of bi-allelic genotype (n(nind) x m(nmar))
#' @param       na.geno   A 2-digit character for missing genotype in monoAlle
#' @return      A matrix of integer
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function transforms the coding from AA/AB/BB to -1/0/1
#' @examples    
#' genodata = rbind(c("AA","AB","--","CC")
#'                 ,c("AA","--","BB","TT")
#'                 ,c("AB","AB","AB","--"))
#' monoAlle = geno.biAL2monoAL(genodata)
#' monoAlle.AB = geno.AGCT2AB(monoAlle)
#' geno.AB = geno.monoAL2biAL(monoAlle.AB)
#' (mat.AB2Num = geno.AB2Num(geno.AB,"--")) 
#' #    C1 C2 C3 C4
#' # R1  1  0 NA  1
#' # R2  1 NA  1 -1
#' # R3  0  0  0 NA
#' 
# -----------------------------------------------------------------------------

geno.AB2Num = function(genodata,na.geno = NULL){
  if ((!is.null(na.geno)) & (any(genodata == na.geno))){
    genodata[(genodata == na.geno)] = NA    
  }
  return((genodata == "AA") - (genodata == "BB"))
}

# -----------------------------------------------------------------------------
# UPDATED Aug 25, 2015 3:58 PM
# FUNCTION:     geno.Num2AB(genodata.int,code=c(1,0,-1),na.geno = NULL)
#' @rdname      geno.AB2Num
#' @param       genodata.int  A matrix of integer for genotype
#' @param       code          A vector of integer for AA/AB/BB
#' @param       na.geno.int       A integer for missing genotype
#' @return      A matrix of bi-allelic genotype
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function transforms the coding from -1/0/1 to BB/AB/AA
#' @examples    
#' genodata = rbind(c("AA","AB","--","CC")
#'                 ,c("AA","--","BB","TT")
#'                 ,c("AB","AB","AB","--"))
#' monoAlle = geno.biAL2monoAL(genodata)
#' monoAlle.AB = geno.AGCT2AB(monoAlle)
#' geno.AB = geno.monoAL2biAL(monoAlle.AB)
#' (mat.AB2Num = geno.AB2Num(geno.AB,"--"))
#' #    C1 C2 C3 C4
#' # R1  1  0 NA  1
#' # R2  1 NA  1 -1
#' # R3  0  0  0 NA
#' 
#' geno.Num2AB(geno.AB2Num(geno.AB,"--"))
#' # C1   C2   C3   C4  
#' # R1 "AA" "AB" NA   "AA"
#' # R2 "AA" NA   "AA" "BB"
#' # R3 "AB" "AB" "AB" NA 
#' 
#' a = geno.AB2Num(geno.AB,"--")
#' a[is.na(a)] = -9
#' geno.Num2AB(a,na.geno = -9)
#' # C1   C2   C3   C4  
#' # R1 "AA" "AB" NA   "AA"
#' # R2 "AA" NA   "AA" "BB"
#' # R3 "AB" "AB" "AB" NA  
# -----------------------------------------------------------------------------
geno.Num2AB = function(genodata.int,
                       code=c(1,0,-1),
                       na.geno.int = NULL){
  if (!is.null(na.geno.int) && na.geno.int %in% code){
    stop("Invalid missing code (ambiguous code)")
  }
  if ((!is.null(na.geno.int)) & (any(genodata.int == na.geno.int))){
    genodata.int[(genodata.int == na.geno.int)] = NA
  }
  if (any(genodata.int == code[3])){
    genodata.int[(genodata.int == code[3])] = "BB"
  }
  if (any(genodata.int == code[2])){
    genodata.int[(genodata.int == code[2])] = "AB"
  }
  if (any(genodata.int == code[1])){
    genodata.int[(genodata.int == code[1])] = "AA"
  }
  genodata.int
}
# -----------------------------------------------------------------------------
# UPDATED Aug 26, 2015 12:55 PM
# FUNCTION:     geno.qc(genodata,code=c(1,0,-1),na.string = NULL)
#' @title       Quality control for genotype data
#' @param       genodata   A matrix of genotype
#' @param       code       A vector of codes used in the genotype matrix
#' @param       na.string  The character used to code the missing value.
#' @return      A data.frame including basic statistics for SNPs.
# -----------------------------------------------------------------------------
#' @export      
#' @note        Quality control for genotype data
#' @examples  
#' genodata_2 = rbind(c("AA","AB","--","AA","--")
#'                   ,c("AA","--","BB","BB","--")
#'                   ,c("AB","AB","AB","--","--"))
#'                   
#' geno.qc(genodata_2,c("AA","AB","BB"),"--") 
#' #   ngeno nmiss     rmiss NumAA NumAB NumBB       AAF       BAF       MAF
#' # 1     3     0 0.0000000     2     1     0 0.8333333 0.1666667 0.1666667
#' # 2     2     1 0.3333333     0     2     0 0.5000000 0.5000000 0.5000000
#' # 3     2     1 0.3333333     0     1     1 0.2500000 0.7500000 0.2500000
#' # 4     2     1 0.3333333     1     0     1 0.5000000 0.5000000 0.5000000
#' # 5     0     3 1.0000000     0     0     0 0.0000000 0.0000000 0.0000000
# -----------------------------------------------------------------------------

geno.qc = function(genodata,code=c(1,0,-1),na.string = NULL){
  # basic variables
  nind = nrow(genodata)
  # missing rate
  if (!is.null(na.string)){
    genodata[genodata == na.string] = NA
  }
  num.missing = colSums(is.na(genodata))
  num.calling = nind - num.missing
  rate.missing = num.missing / nind
  # genotype frequency
  num.AA = colSums(genodata == code[1], na.rm=T)
  num.AB = colSums(genodata == code[2], na.rm=T)
  num.BB = colSums(genodata == code[3], na.rm=T)
  # allele frequency
  num.A.allele = num.AA * 2 + num.AB
  num.B.allele = num.BB * 2 + num.AB
  num.calling.nozero = num.calling
  num.calling.nozero[num.calling.nozero == 0] = 10^14
  freq.A.allele = num.A.allele / (num.calling.nozero * 2)
  freq.B.allele = num.B.allele / (num.calling.nozero * 2)
  minor.allele.freq = pmin(freq.A.allele,freq.B.allele)
  # output
  res.freq = data.frame(ngeno  = num.calling,
                        nmiss  = num.missing,
                        rmiss  = rate.missing,
                        NumAA  = num.AA,
                        NumAB  = num.AB,
                        NumBB  = num.BB,
                        AAF    = freq.A.allele,
                        BAF    = freq.B.allele,
                        MAF    = minor.allele.freq,
                        row.names        = colnames(genodata),
                        stringsAsFactors = FALSE)
  return(res.freq)
}

# -----------------------------------------------------------------------------
# UPDATED Aug 25, 2015 3:19 PM
# FUNCTION:     geno.impu.ave(genodata)
#' @title       Impute the missing cells with average
#' @param       genodata A numeric matrix of genotype (missing coded as NA)
#' @return      A numeric matrix of genotype
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function imputes missing value in genotype data using
#'              the average value of that column. This method can only be
#'              used for genotype data coded as integer (-1/0/1 or 0/1/2 or 
#'              -10/0/10). It is usually used in regression methods 
#' @examples   
#' genodata = rbind(c("AA","AB","--","CC")
#'                 ,c("AA","--","BB","TT")
#'                 ,c("AB","AB","AB","--"))
#' monoAlle = geno.biAL2monoAL(genodata)
#' monoAlle.AB = geno.AGCT2AB(monoAlle)
#' geno.AB = geno.monoAL2biAL(monoAlle.AB)
#' (mat.AB2Num = geno.AB2Num(geno.AB,"--"))
#' # C1 C2 C3 C4
#' # R1  1  0 NA  1
#' # R2  1 NA  1 -1
#' # R3  0  0  0 NA
#' geno.impu.ave(mat.AB2Num)
#' #    C1 C2  C3 C4
#' # R1  1  0 0.5  1
#' # R2  1  0 1.0 -1
#' # R3  0  0 0.0  0
# -----------------------------------------------------------------------------

  geno.impu.ave = function(genodata){
      pop_mean = matrix(colMeans(genodata,na.rm=T),nrow = 1)
      pop_mean.mat = matrix(1,nrow(genodata),1) %*% pop_mean
      mat_mean = (is.na(genodata)) * pop_mean.mat
      genodata[is.na(genodata)] = 0
      return(genodata + mat_mean)
  }

# -----------------------------------------------------------------------------
# UPDATED Aug 25, 2015 5:28 PM
# FUNCTION:     geno.mono2PED(monodata,IndID,lite,FamilyID,PID,MID,Sex,Pheno)
#' @title       Create PED files
#' @param       monodata A monoAL object
#' @param       IndID    A vector of Individual ID
#' @param       lite     If using a lite version of PED file 
#' @param       FamilyID A vector of family ID
#' @param       PID      A vector of paternal ID
#' @param       MID      A vector of maternal ID
#' @param       Sex      A vector of individual sex
#' @param       Pheno    A vector of phenotypic value
#' @return      A data frame that can be written as a PED file
# -----------------------------------------------------------------------------
#' @export      
#' @note        This function does some simple tests and combine available data
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Jul 10, 2015 19:37
# By default, each line of the MAP file describes a single marker and must 
# contain exactly 4 columns: 
#     chromosome (1-22, X, Y or 0 if unplaced) 
#     rs# or snp identifier 
#     Genetic distance (morgans) 
#     Base-pair position (bp units)
# Genetic distance can be specified in centimorgans with the --cm flag. 
# Alternatively, you can use a MAP file with the genetic distance 
# excluded by adding the flag --map3
# -----------------------------------------------------------------------------

  geno.mono2PED = function(monodata,
                           IndID,
                           lite     = T,
                           FamilyID = NULL,
                           PID      = NULL,
                           MID      = NULL,
                           Sex      = NULL,
                           Pheno    = NULL
                           ){
    # check the animal id
    if (any(IndID != monodata$names.ind)){ 
      stop("wrong ID!")
    }
    
    # Geno data
    Geno = matrix("",monodata$dim[1],monodata$dim[2]*2)
    Geno[,c(TRUE,FALSE)] = monodata$Allele_1
    Geno[,c(FALSE,TRUE)] = monodata$Allele_2
    
    if (lite){
      # result
      res.PED = data.frame(IndID    = IndID,
                           Sex      = Sex,
                           Pheno    = Pheno,
                           Geno     = Geno,
                           stringsAsFactors = FALSE)
    } else {
      # FAKE info
      nind = length(IndID)
      if (is.null(FamilyID)) FamilyID = paste0("Family_",seq(nind))
      if (is.null(PID))      PID      = rep(0,nind)
      if (is.null(MID))      MID      = rep(0,nind)
      if (is.null(Sex))      Sex      = rep(1,nind)
      if (is.null(Pheno))    Pheno    = rep(1,nind)
      # result
      res.PED = data.frame(FamilyID = FamilyID,
                           IndID    = IndID,
                           PID      = PID,
                           MID      = MID,
                           Sex      = Sex,
                           Pheno    = Pheno,
                           Geno     = Geno,
                           stringsAsFactors = FALSE
                           )
      }
    return(res.PED)
  }

# -----------------------------------------------------------------------------
# Apr 20, 2015 19:54        rrBLUP
# GENO-DATA (with colnames c("marker","chrom","pos",ID.Ind))
#   Col 1 Marker ID
#   Col 2 Chromosome
#   Col 3 position
#   other Cols: Genotypes with code [1/0/-1/NA]
# PHENO-DATA
#   Col 1 ID.Ind
#   Including all phenotypes and fixed effects (how about covariate?)
# -----------------------------------------------------------------------------
