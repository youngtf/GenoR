#------------------------------------------------------------------------------
# R file  --Tianfu Yang
# Type:              Functions
# Subtype/Project:   GenoCleaning
# Descripetions:     functions used to clean the genotype data
# Last Update:       2014-05-20
# Contents:
# =============================================================================
#
#    ______     
#   | "TG" |  <-------------------------
#   |______|  ------------------------  | 
#      | ^                         <5>| | <5>
#  <1> | | <2>                        | |
#      v |                            | v
#    ______        ______           ______        __         __________ 
#   |"T""G"| ---> |"0""1"| ------> |-1/0/1| ---> |QC| ----> |Imputation| 
#   |______| <3>  |______|   <4>   |______| <6>  |__|  <7>  |__________|
#      |              |                                          |
#      |<8>           |<8>                                       |
#      v              v                                          v
#    ______        ______                                    __________
#   | PED  |      | PED  |                                  |  Gensel  |
#   |______|      |______|                                  |__________|
#
# =============================================================================
# 1                                                       Updated (MAR10,2015)
#   Function:    geno.biAL2monoAL(genodata,na.geno,na.output)
#   Description: divide bi-allelic genotype data into a mono-allelic one
# 2                                                       Updated (MAR10,2015)
#   Function:    geno.monoAL2biAL(monoAlle,na.geno,na.output)
#   Description: combine the allele data into a Bi-allele data
# 3                                                       Updated (MAR10,2015)
#   Function:    geno.AGCT2AB(monoAlle,na.output)
#   Description: Transform the coding of monoAL data to 1/2, which is the 
#                index in the allele list.
# 4                                                       Updated (MAR10,2015)
#   Function:    geno.monoAL2Num(monoAlle_AB,na.geno,na.output)
#   Description: transform the coding from mono-allele to -1/0/1/-5
# 5
#   Function:    geno.AB2Num(genodata)
#   Description: transform the coding from AA/AB/BB tp -1/0/1
# 6
#   Function:    geno.QC(genodata)
#   Description: Quality control for genotype data
# 7 
#   Function:    geno.impu.ave(genodata)
#   Description: Impute the missing cells with average
# 8
#   Function:    geno.mono2PED(genodata)
#   Description: Create PED files
# 9 
#   Function:    geno.biAL2rrb(genodata)
#   Description: Create geno and pheno files for R/rrBLUP
#------------------------------------------------------------------------------

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
#------------------------------------------------------------------------------
#- 1 --------------------------------------------------------------------------
# Updated Feb 18, 2015 16:12
# Function:    geno.biAL2monoAL(genodata,na.geno,na.output)
# Description: divide bi-allelic genotype data into mono-allelic
# input:   genodata: n(nind) x m(nmar)
# ouput:   A LIST (see unit test)
#------------------------------------------------------------------------------
### Unit test -----------------------------------------------------------------
# genodata = rbind(c("AA","AB","--")
#                 ,c("AA","--","BB")
#                 ,c("AB","AB","AB"))
# monoAlle = geno.biAL2monoAL(genodata)
# monoAlle
#
# $dim
# Row Col 
# 3   3 
# 
# $names.ind
# [1] "R1" "R2" "R3"        # auto naming
# 
# $names.mar
# [1] "C1" "C2" "C3"        # auto naming
# 
# $names.base
#   $names.base$C1
#   [1] "A" "B"
# 
#   $names.base$C2
#   [1] "A" "B"
# 
#   $names.base$C3
#   [1] "B" "A"
# 
# 
# $Allele_1
# C1  C2  C3 
# R1 "A" "A" "-"
# R2 "A" "-" "B"
# R3 "A" "A" "A"
# 
# $Allele_2
# C1  C2  C3 
# R1 "A" "B" "-"
# R2 "A" "-" "B"
# R3 "B" "B" "B"
#
#==============================================================================

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
    warning("more than 2 possible alleles: ",idx.innormal,)
  }
  ## result: allele list
    res.data$names.base = allele.list
  ## Attributes
    attributes(res.data)$class = "monoAL"
  ## Return  
    res.data
}

#- 2 --------------------------------------------------------------------------
# Updated Mar 10, 2015 14:55
# Function:    geno.monoAL2biAL(monoAlle,na.geno,na.output)
# Description: combine the allele data into a Bi-allele data
# input:   monoAL data
# ouput:   A matrix (see unit test)
#------------------------------------------------------------------------------
### Unit test -----------------------------------------------------------------
# genodata = rbind(c("AA","AB","--")
#                 ,c("AA","--","BB")
#                 ,c("AB","AB","AB"))
# monoAlle = geno.biAL2monoAL(genodata)
# geno.monoAL2biAL(monoAlle)
#
#     C1   C2   C3  
# R1 "AA" "AB" "--"
# R2 "AA" "--" "BB"
# R3 "AB" "AB" "AB"
### ---------------------------------------------------------------------------

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

#- 3 --------------------------------------------------------------------------
# Updated Mar 10, 2015 15:57
# Function:     geno.AGCT2AB(monoAlle,na.output)
# Description:  Transform the coding of monoAL data to 1/2/0, which is the 
#               index in the allele list / missing.
#               the A/B code is only determined by the allele list in the 
#               raw monoAL data (which is also kept in the result)
# input:   monoAL data
# ouput:   monoAL / monoAL_AB data
#------------------------------------------------------------------------------
### Unit test -----------------------------------------------------------------
# genodata = rbind(c("AA","AB","--")
#                 ,c("AA","--","BB")
#                 ,c("AB","AB","AB"))
# monoAlle = geno.biAL2monoAL(genodata)
# geno.AGCT2AB(monoAlle)
# 
# $dim
# Row Col 
# 3   3 
# 
# $names.ind
# [1] "R1" "R2" "R3"
# 
# $names.mar
# [1] "C1" "C2" "C3"
# 
# $names.base
# $names.base$C1
# [1] "A" "B"
# 
# $names.base$C2
# [1] "A" "B"
# 
# $names.base$C3
# [1] "B" "A"
# 
# 
# $Allele_1
# [,1] [,2] [,3]
# [1,]    1    1    0
# [2,]    1    0    1
# [3,]    1    1    2
# 
# $Allele_2
# [,1] [,2] [,3]
# [1,]    1    2    0
# [2,]    1    0    1
# [3,]    2    2    1
# 
# attr(,"class")
# [1] "monoAL"    "monoAL_AB"
### ---------------------------------------------------------------------------
geno.AGCT2AB = function(monoAlle,na.output = 0 ){
  
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
                             ,function(x) match(monoAlle$Allele_1[,x]
                                               ,allele_List[[x]])
                             )
    res.list$Allele_2 = apply(t(seq(nmar)), 2
                             ,function(x) match(monoAlle$Allele_2[,x]
                                               ,allele_List[[x]])
                             )
  ## All NA is missing
    res.list$Allele_1[is.na(res.list$Allele_1)] = na.output
    res.list$Allele_2[is.na(res.list$Allele_2)] = na.output
  ## Attributes
    attributes(res.list)$class = c("monoAL"
                                   ,"monoAL_AB" # implies its coding format
                                  )
  ## return the result
    res.list
}

#- 4 --------------------------------------------------------------------------
# Updated Mar 10, 2015 16:27
# Function:    geno.monoAL2Num(monoAlle_AB,na.geno,na.output)
# Description: transform the coding from mono-allele to -1/0/1/-5
# input:   monoAL_AB data
# ouput:   A matrix (see unit test)
#------------------------------------------------------------------------------
### Unit test -----------------------------------------------------------------
# genodata = rbind(c("AA","AB","--")
#                 ,c("AA","--","BB")
#                 ,c("AB","AB","AB"))
# monoAlle    = geno.biAL2monoAL(genodata)
# monoAlle_AB = geno.AGCT2AB(monoAlle)
# geno.monoAL2Num(monoAlle_AB)
# 
#    C1 C2 C3
# R1  1  0 -5
# R2  1 -5  1
# R3  0  0  0
### ---------------------------------------------------------------------------
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

#- 5 --------------------------------------------------------------------------
# Function:    AB2Num(genodata)
# Description: transform the coding from AA/AB/BB tp -1/0/1
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Mar 13, 2015 14:32
# 1 Don't use it until it can handle missing data - check!
#------------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Unit test 
# genodata = rbind(c("AA","AB","--")
#                 ,c("AA","--","BB")
#                 ,c("AB","AB","AB"))
# geno.AB2Num(genodata)
# #      [,1] [,2] [,3]
# # [1,]    1    0    0
# # [2,]    1    0   -1
# # [3,]    0    0    0
# geno.AB2Num(genodata,"--")
# #      [,1] [,2] [,3]
# # [1,]    1    0   NA
# # [2,]    1   NA   -1
# # [3,]    0    0    0
# geno.Num2AB(geno.AB2Num(genodata,"--"))
# #      [,1] [,2] [,3]
# # [1,] "AA" "AB" NA  
# # [2,] "AA" NA   "BB"
# # [3,] "AB" "AB" "AB"
# a = geno.AB2Num(genodata)
# a[2,2] = -9
# a
# #      [,1] [,2] [,3]
# # [1,]    1    0    0
# # [2,]    1   -9   -1
# # [3,]    0    0    0
# geno.Num2AB(a,missing = -9)
# #      [,1] [,2] [,3]
# # [1,] "AA" "AB" "AB"
# # [2,] "AA" NA   "BB"
# # [3,] "AB" "AB" "AB"
# -----------------------------------------------------------------------------

geno.AB2Num = function(genodata,missing = NULL){
  if ((!is.null(missing)) & (any(genodata == missing))){
    genodata[(genodata == missing)] = NA    
  }
  return((genodata == "AA") - (genodata == "BB"))
}
geno.Num2AB = function(genodata,
                       code=c(1,0,-1),
                       missing = NULL){
  if (missing %in% code){
    stop("Invalid missing code (ambiguous code)")
  }
  if ((!is.null(missing)) & (any(genodata == missing))){
    genodata[(genodata == missing)] = NA
  }
  if (any(genodata == code[3])){
    genodata[(genodata == code[3])] = "BB"
  }
  if (any(genodata == code[2])){
    genodata[(genodata == code[2])] = "AB"
  }
  if (any(genodata == code[1])){
    genodata[(genodata == code[1])] = "AA"
  }
  genodata
}
#- 6 --------------------------------------------------------------------------
# Function:    geno.QC(genodata)
# Description: Quality control for genotype data
#------------------------------------------------------------------------------
geno.QC = function(genodata){

  # make a matrix, 1 for 1/0/-1 ,and 0 for NA or any other value
  mat.geno = (genodata == 1) + (genodata == 0) + (genodata == -1)
  # number of valid genotypes
  ngeno = colSums(mat.geno,na.rm=T) 

  EnAA_AB_BB = cbind((AAF^2) * nValid 
                      ,(AAF * BAF * 2) * nValid
                      ,(BAF^2) * nValid
                    )
  res = colSums(((nAA_AB_BB - EnAA_AB_BB)^2) / EnAA_AB_BB)
}


geno.QC.Missing = function(genodata){
    nind = nrow(genodata)
    rate.missing = colSums(is.na(genodata)) / nind
    rate.missing
}

geno.QC.freq = function(genodata){                 # Major genotype frequency
  # make a matrix, 1 for 1/0/-1 ,and 0 for NA or any other value
  mat.geno = (genodata == 1) + (genodata == 0) + (genodata == -1)
  # number of valid genotypes
  ngeno = colSums(mat.geno,na.rm=T) 
  #------------------------------------------------------------------------------
  # Mar 17, 2015 11:02
  # 1 when ngeno == 0 the MAF may be not valid (which may be removed in "missing")
  #------------------------------------------------------------------------------
  mat.alle.A   = (genodata == 1) * 2 + (genodata == 0)            # n.allele A
  freq.alle.A  = colSums(mat.alle.A,na.rm=T) / (2 * ngeno)
  freq.alle.B  = 1 - freq.alle.A
  MAF = pmin(freq.alle.A,freq.alle.B)

  freq.alle.AA  = colSums(genodata == 1, na.rm=T)
  freq.alle.AB  = colSums(genodata == 0, na.rm=T)
  freq.alle.BB  = colSums(genodata == -1,na.rm=T)
  MGF = pmax(freq.alle.AA,freq.alle.AB,freq.alle.BB)
  
  res.freq = data.frame(ngeno  = ngeno,
                        AAF    = freq.alle.A,
                        BAF    = freq.alle.B,
                        MAF    = MAF,
                        FreqAA = freq.alle.AA,
                        FreqAB = freq.alle.AB,
                        FreqBB = freq.alle.BB,
                        MGF    = MGF,
                        row.names        = colnames(genodata),
                        stringsAsFactors = FALSE
                        )
}

#- 7 --------------------------------------------------------------------------
# Function:    geno.impu.ave(genodata)
# Description: Impute the missing cells with average
#------------------------------------------------------------------------------
  geno.impu.ave = function(genodata){
      pop_mean = colMeans(genodata,na.rm=T)
      mat_mean = (is.na(genodata)) * repmat(pop_mean,nrow(genodata),1)
      genodata[is.na(genodata)] = 0
      genodata + mat_mean
  }
#- 8 --------------------------------------------------------------------------
# Updated Feb 18, 2015 16:12
# Function:    geno.mono2PED     (monodata,FamilyID,IndID,PID,MID,Sex,Pheno)
# Function:    geno.mono2PED.lite(monodata,         IndID,        Sex,Pheno)
# Description: Create PED files
#------------------------------------------------------------------------------
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
                      FamilyID = NULL,
                      IndID    = NULL,
                      PID      = NULL,
                      MID      = NULL,
                      Sex      = NULL,
                      Pheno    = NULL
                      ){
    # check the animal id
    if (any(IndID != monodata$names.ind)) stop("wrong ID!")
    Geno = matrix("",monodata$dim[1],monodata$dim[2]*2)
    Geno[,c(TRUE,FALSE)] = monodata$Allele_1
    Geno[,c(FALSE,TRUE)] = monodata$Allele_2
    
    if (is.null(FamilyID)) FamilyID = paste0("Family_",seq(nind))
    if (is.null(PID))      PID      = rep(0,nind)
    if (is.null(MID))      MID      = rep(0,nind)
    if (is.null(Sex))      Sex      = rep(1,nind)
    if (is.null(Pheno))    Pheno    = rep(1,nind)
    
    if (is.null(IndID)){
      IndID.Raw = rownames(monodata)[1:nind]
      len.names = nchar(IndID.Raw)   
      IndID = substr(IndID.Raw,1,(len.names-2))
    }
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

  geno.mono2PED.lite = function(monodata,
                                IndID,
                                Sex    = 1,
                                Pheno  = 1
                                ){
    # check the animal id
    if (any(IndID != monodata$names.ind)) stop("wrong ID!")
    Geno = matrix("",monodata$dim[1],monodata$dim[2]*2)
    Geno[,c(TRUE,FALSE)] = monodata$Allele_1
    Geno[,c(FALSE,TRUE)] = monodata$Allele_2
    res.PED = data.frame(IndID    = IndID,
                         Sex      = Sex,
                         Pheno    = Pheno,
                         Geno     = Geno,
                         stringsAsFactors = FALSE
    )
  }

#- 9 --------------------------------------------------------------------------
# Updated Apr 20, 2015 19:53
# Function:    geno.biAL2rrb(monodata,map)
# Description: Create PED files
#------------------------------------------------------------------------------
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
  