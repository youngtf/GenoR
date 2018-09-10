# -----------------------------------------------------------------------------
# R FILE             --Tianfu Yang
# Type:              Functions/Analysis/Visualization
# Subtype/Project:   
# Descriptions:    
# -----------------------------------------------------------------------------
# Contents:
# -----------------------------------------------------------------------------
# To-do
# -----------------------------------------------------------------------------
# Pre-load
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FUNCTION:     create_genomic_map
#' DESCRIPTION
#' @rdname      # the name of the rd file
#' @title       # the title of the function
#' @param       \
#' @return      \item{name a}{description} \
# -----------------------------------------------------------------------------
#' @export      \
#' @note        \
#' @examples    \
# -----------------------------------------------------------------------------
create_genomic_map = function(
  SNPID, 
  chromosome, 
  position, 
  allosome, 
  autosome     = NULL,
  mitochondria = NULL
  ){
  chr_factor = naturalfactor(chromosome)
  res_map = as.tbl(
    data.frame(
      SNPID            = SNPID, 
      CHR              = chr_factor, 
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
  
  res_map_sorted = res_map[order(res_map$CHR, res_map$POS),]
  
  return(res_map_sorted)
}
