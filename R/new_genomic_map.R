# -----------------------------------------------------------------------------
# R FILE             --Tianfu Yang
# Type:              Functions/Analysis/Visualization
# Subtype/Project:   
# Descriptions:    
# -----------------------------------------------------------------------------
# Contents:
# -----------------------------------------------------------------------------
# To-do
# create_genomic_map: lost attributes
# geno_map_check_dup: print too long
# -----------------------------------------------------------------------------
# Pre-load
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FUNCTION:     create_genomic_map
#' DESCRIPTION
#' @rdname      create_genomic_map
#' @title       Create a geno.map object
#' @param       SNPID a character vector of snp ID
#' @param       chromosome a character vector of chromosome ID
#' @param       position a integer vector of physical position
#' @param       allosome a character vector of allosomes, default: excluded
#' @param       autosome a character vector of autosomes
#' @param       mitochondria a character vector of mitochondria, default: excluded
#' @return      a geno.map object
# -----------------------------------------------------------------------------
#' @export      
#' @note        A "factory unction" for geno.map object
# @examples    
# -----------------------------------------------------------------------------
create_genomic_map = function(
  SNPID, 
  chromosome, 
  position, 
  autosome,
  allosome     = NULL, 
  mitochondria = NULL
){
  res_map = as.tbl(
    data.frame(
      SNPID            = SNPID, 
      CHR              = chromosome, 
      POS              = position, 
      stringsAsFactors = F)
  )

  ## add attributes
  attributes(res_map)$chrom_id     = autosome
  attributes(res_map)$allosome     = autosome
  if (!is.null(allosome)){
    attributes(res_map)$allosome   = allosome
    attributes(res_map)$chrom_id   = c(attributes(res_map)$chrom_id, 
                                       allosome)
  }
  if (!is.null(mitochondria)){
    attributes(res_map)$MT         = mitochondria
    attributes(res_map)$chrom_id   = c(attributes(res_map)$chrom_id, 
                                       mitochondria)
  }
  attributes(res_map)$n_chromosome = length(attributes(res_map)$chrom_id)
  
  ## remove unmapped
  res_map = filter(res_map, CHR %in% attributes(res_map)$chrom_id)
  attributes(res_map)$chr_len      = tapply(res_map$POS, res_map$CHR, max)
  
  res_map$CHR = naturalfactor(res_map$CHR)
  res_map_sorted = res_map[order(res_map$CHR, res_map$POS),]
  
  return(res_map_sorted)
}

# -----------------------------------------------------------------------------
# FUNCTION:     geno_map_check_dup
#' DESCRIPTION
#' @rdname      geno_map_check_dup
#' @title       check duplicated positions
#' @param       geno_map a geno_map object
#' @return      a geno_map object
# -----------------------------------------------------------------------------
#' @export      
#' @note        check and randomly remove duplicated positions. 
# @examples    
# -----------------------------------------------------------------------------

geno_map_check_dup = function(geno_map){
  if_dup = duplicated(geno_map[,c("CHR","POS")])
  print("The following SNPs were removed due to duplicated positions:")
  for (idx in which(if_dup)){
    print(geno_map[idx,])
  }
  new_geno_map = geno_map[!if_dup,]
  return(new_geno_map)
}

