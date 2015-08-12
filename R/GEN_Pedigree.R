# -----------------------------------------------------------------------------
# R file             --Tianfu Yang
# Type:              Functions
# Subtype/Project:   Pedigree-related functions
# Descriptions:      These are some functions for preparing pedigree data, such 
#                    as sorting, checking, and building additive relationship 
#                    matrix.
#
# |----------|   |--------------|   |----------|   |------|   |-----------|
# | Ped info | > | Pedigree Obj | > | CHECK ID | > | Sort | > | CHECK SEX |
# |----------|   |--------------|   |----------|   |------|   |-----------|
#
# Last Update:       Jul 29, 2015 5:14 PM
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Updated Jul 28, 2015 5:20 PM
# Function:     Gen_Pedigree(AID, DID, SID, SEX, na.string = NULL)
# Description:  The basic function to generate a Pedigree object
# input:        AID                # ID of Offspring
#               DID                # ID of Dams
#               SID                # ID of Sires
#               SEX                # Sex: F/M/U for female/male/unknown
#               na.string = NULL   # string for unknown ID, NA for default
# ouput:        Pedigree data
# -----------------------------------------------------------------------------
Gen_Pedigree = function(AID, DID, SID, 
                        SEX = c("F","M","U"), 
                        na_string = NULL,
                        check_ID  = TRUE,
                        check_sex = TRUE,
                        sort_data = TRUE){
  # check the data
  if (any(!(SEX %in% c("F","M","U")))){
    stop("Please code the SEX as F/M/U for female/male/unknown.")
  }
  if (!(length(AID) == length(DID) & 
        length(DID) == length(SID) &
        length(SID) == length(SEX))){
    stop(paste("AID, DID and SID should have the same length,",
                "but the input vectors have length:",
                length(AID), length(DID), length(SID),length(SEX)))
  }
  if (any(is.na(AID))){
    stop("AID should all be valid string and no NA is allowed")
  }
  if (!is.null(na_string) & any(AID == na_string)){
    stop("AID should all be valid string and no NA is allowed")
  }
  # build Pedigree
  res_ped = data.frame(AID = as.character(AID),
                       DID = as.character(DID),
                       SID = as.character(SID),
                       SEX = as.character(SEX),
                       stringsAsFactors = FALSE) 
  if (!is.null(na_string)){
    message("All ID that code as",na_string,"were set as missing.")
    res_ped[res_ped == na_string] = NA
  }
  attributes(res_ped)$class = c("Pedigree","data.frame") 
  # output
  return(res_ped)
}
# -----------------------------------------------------------------------------
# Updated Jul 28, 2015 5:34 PM
# Function:     Gen_Ped_check_ID(pedigree_data)
# Description:  Check the ID of animals:
#               1. similar ID;
#               2. all blanks
#               3. end with blanks.
# input:        pedigree_data, 
# ouput:        pedigree_data (checked)
# -----------------------------------------------------------------------------
Gen_Ped_check_ID = function(pedigree_data,
                            similar_ID = TRUE,
                            with_blanks = TRUE,
                            end_with_blanks = FALSE){
  message("Checking the ID of animals...")
  # check the data 
  if (!inherits(pedigree_data,"Pedigree")){
    stop(paste("The input data should be a Pedigree object,", 
               "which is likely obtained from function Gen_Pedigree() "))
  }
  res_check_id = vector("list",0)
  
  # all ID
  ID_all = unique(c(pedigree_data$DID,pedigree_data$SID,pedigree_data$AID))
  ID_all = ID_all[!is.na(ID_all)]
  n_ID = length(ID_all)
  nchar_ID = nchar(ID_all)
  
  # with blanks
  if (with_blanks){
    idx.with_blank = grep(" ", ID_all)
    cat("  Number of ID(s) with blanks:",length(idx.with_blank),"\n")
    res_check_id$with_blank = ID_all[idx.with_blank]
  }
  # end with blanks
  if(end_with_blanks){
    last_letter = substr(ID_all,nchar_ID,nchar_ID)
    idx.end_with_blank = which(last_letter == " ")
    cat("  Number of ID(s) ending with blanks:",length(idx.end_with_blank),"\n")
    res_check_id$end_with_blank = ID_all[idx.end_with_blank]
  }
  # similar_ID
  if(similar_ID){
    res_similar_ID = data.frame(ID_A = character(0),ID_B = character(0))
    for (i in 1:(n_ID - 1)){
      ID_A = ID_all[i]
      ID_similar_temp = agrep(ID_A, ID_all[-i], 
                              max.distance = list(cost = 0,
                                                  insertions = 0, 
                                                  deletions = 2,
                                                  substitutions = 0),
                              ignore.case = TRUE,value = TRUE)
      if (length(ID_similar_temp) > 0){
        ID_similar_temp_df = data.frame(ID_A,ID_similar_temp,
                                        stringsAsFactors = FALSE)
        res_similar_ID = rbind(res_similar_ID,ID_similar_temp_df)
      }
    }
    cat("Number of pairs of similar IDs:",nrow(res_similar_ID),"\n")
    res_check_id$similar_ID = res_similar_ID
  }
  message("Checking the ID of animals...Done!")
  return(res_check_id)
}

# -----------------------------------------------------------------------------
# Updated Jul 28, 2015 5:34 PM
# Function:     Gen_Ped_sort(pedigree_data)
# Description:  Sort Pedigree object (parents first).
# input:        pedigree_data
# ouput:        pedigree_data (sorted)
# -----------------------------------------------------------------------------
Gen_Ped_sort = function(pedigree_data){
  message("Sorting the entries of relationship...")
  # check the data 
  if (!inherits(pedigree_data,"Pedigree")){
    stop(paste("The input data should be a Pedigree object,", 
               "which is likely obtained from function Gen_Pedigree() "))
  }
  # get all aninmals
  ID_all = unique(c(pedigree_data$DID,pedigree_data$SID,pedigree_data$AID))
  
  message("Done!")
  return(res_ped)
}
# -----------------------------------------------------------------------------
# Updated Jul 28, 2015 5:34 PM
# Function:     Gen_Ped_check_sex(pedigree_data, auto_fill = TRUE)
# Description:  Check the Sex of animals
# input:        pedigree_data (should be sorted first)
#               auto_fill = TRUE
# ouput:        pedigree_data (checked)
# -----------------------------------------------------------------------------
Gen_Ped_check_sex = function(pedigree_data, auto_fill = TRUE){
  message("Checking the sex of animals...")
  # check the data 
  if (!inherits(pedigree_data,"Pedigree")){
    stop(paste("The input data should be a Pedigree object,", 
               "which is likely obtained from function Gen_Pedigree() "))
  }
  Sire_ID   = as.character(na.omit(pedigree_data$SID))
  Dam_ID    = as.character(na.omit(pedigree_data$DID))
  male_ID   = as.character(na.omit(pedigree_data$AID[pedigree_data$SEX == "M"]))
  female_ID = as.character(na.omit(pedigree_data$AID[pedigree_data$SEX == "F"]))
  u_sex_ID  = as.character(na.omit(pedigree_data$AID[pedigree_data$SEX == "U"]))
  
  res_check_sex = vector("list",0)
  
  # check the column and sex
  Sires_are_male = (Sire_ID %in% male_ID)
  if (!(all(Sires_are_male,na.rm = T))){
    Sires_are_not_male = Sire_ID[!Sires_are_male]
    cat("  Sire(s) with contradictory sex:", Sires_are_not_male)
    res_check_sex$Sires_are_not_male = Sires_are_not_male
  }
  Dams_are_female = Dam_ID %in% female_ID
  if (!(all(Dams_are_female,na.rm = T))){
    Dams_are_not_female = Dam_ID[!Dams_are_female]
    cat("  Dam(s) with contradictory sex:", Dams_are_not_female)
    res_check_sex$Dams_are_not_female = Dams_are_not_female
  }
  # check the overlap of two column
  Both_Sire_Dam = intersect(Sire_ID,Dam_ID)
  if (length(Both_Sire_Dam > 0)){
    cat("  Animals in both Sire column and Dam column:", Both_Sire_Dam)
    res_check_sex$Both_Sire_Dam = Both_Sire_Dam
  }
  if (auto_fill & length(u_sex_ID) > 0){
    res_check_sex$possible_male = u_sex_ID[u_sex_ID %in% Sire_ID]
    res_check_sex$possible_female = u_sex_ID[u_sex_ID %in% Dam_ID]
  }
  message("Checking the sex of animals...Done!")
  return(res_check_sex)
}

