# -----------------------------------------------------------------------------
# R file             -- Tianfu Yang
# Type:              Functions
# Subtype/Project:   Survival analysis
# Descriptions:      
# Last Update:       Mar 24, 2015
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Updated Mar 24, 2015 12:56
# Function:     SUR.RCWeight
# Description:  Calculating the weight in Right-Censored data analysis
# input:        
# ouput:        
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Unit test 
# df.pheno = data.frame(AID   = as.character(1:5),
#                       pheno = c(9,7,5,3,1),
#                       IF    = c(1,1,0,0,1))
# SUR.RCWeight(df.pheno)
# -----------------------------------------------------------------------------
SUR.RCWeight = function(df.pheno  # df 3 cols: AID, pheno and IF.NONcensor 
                    ){
  # Check data:
  if (any(names(df.pheno) != c("AID","pheno","IF"))) {
    stop("Please name the factors in the data frame as 'AID','pheno' and 'IF'.")
  }
  
  # basic info
  n = nrow(df.pheno)
  
  # reorder the observation
  idx.pheno.order = order(df.pheno$pheno)
  df.pheno.order = df.pheno[idx.pheno.order,]
  
  # initial weights
  w = vector("numeric",n)
  w[1] = df.pheno.order$IF[1] / n
  
  # calculate other weights
  for (i in 2:n) {
    p = df.pheno.order$IF[i] / (n - i + 1)
    for (j in 1:(i-1)) {
      p = p * ((n - j)/(n - j + 1)) ^ df.pheno.order$IF[j]
    }
    w[i] = p
  }
  
  # set result
  res = df.pheno.order
  res$weight = w

  res
}
