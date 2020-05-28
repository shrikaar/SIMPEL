#' Function to process XCMS data and start getting isotopic information
#' from the formula an other information
#'
#' This function allows you to express your love of cats.
#' @param XCMS_data Do you love cats? Defaults to TRUE.
#' @param compounds_data Do you love cats? Defaults to TRUE.
#' @keywords metabolomics
#' @export
#' @examples
#' get_table_objects()

get_table_objects <- function(XCMS_data, compounds_data){



  XCMS_data <- XCMS_data %>%
    mutate(comp_result = NA)

  XCMS_data <- XCMS_data %>%
    mutate(carbon = NA)

  XCMS_data <- XCMS_data %>%
    mutate(nitrogen = NA)

  XCMS_data <- XCMS_data %>%
    mutate(isotope_numbers = NA)

  XCMS_data <- XCMS_data %>%
    mutate(Bin = NA)


  XCMS_data <- XCMS_data %>%
    mutate(ID = NA)



  #looping through the list of isotopologs we're working with
  for(i in 1:nrow(compounds_data))
  {
    #for each isotopolog, we're going to see if it's in the data
    #the XCMS data, that is
    compounds_data[i, ] %>%
      unlist() %>%
      as.vector() %>%
      print()




    ##for each compound in the metadata (id'd by rt and mass)
    #determine its formula and mass
    comp_name <- as.character(compounds_data[i, "formula"])


    r_time    <- as.numeric(compounds_data[i, "rt"])

    #send that mass to look up all of the C and N isotopes
    #and get the upper and lower bounds on each one
    comp_lookup_table <- get_comp_mz_lookup(compounds_data, comp_name, r_time, ppm)

    #information for the bin as found in the prefix
    myBin = as.character(compounds_data[i, "prefix"])


    myID = as.character(compounds_data[i, "formula"])

    #looping through the xcms data

    #for each feature in the mass spec data, id'd by RT and mz
    #determine if it is one of the C or N isotopes by whether or not
    #it falls within the range of the lower and upper bounds which
    #are stored in comp_lookup_table
    for(j in 1:nrow(XCMS_data)){
      if(is.na(XCMS_data[j, "comp_result"]))
      {
        #print(j)
        x = as.numeric(XCMS_data[j, "mz"])
        y = as.numeric(XCMS_data[j, "rt"])
        #determine if this feature, id's by mz and RT is one of the C and N isotopes
        #whose lower and upper bounds are stored in comp_lookup_table
        #print("in the second loop")

        #match up the xcms an compound_data tables

        #add the number of N's and the number of C's
        val <- get_comp_stage(x, y, comp_lookup_table, r_time, rt_tolerance)
        #print(val)
        #print("is val")

        #val <- ifelse(is.null(val), NA, val)

        if(is.null(val) == FALSE)
        {
          XCMS_data[j, "comp_result"] <- val$compound
          XCMS_data[j, "carbon"] <- val$carbon
          XCMS_data[j, "nitrogen"] <- val$nitrogen

          isotopeNumbers = val$carbon$carbon + val$nitrogen$nitrogen

          XCMS_data[j, "isotope_numbers"] <- isotopeNumbers
          XCMS_data[j, "Bin"] <- myBin
          XCMS_data[j, "ID"] <- myID
          print("we have isotopes now")
          print(i)
          print("is i")
          print(j)
          print("is j")
          print(val)
          print("is val")
        }
        if(is.null(val) == TRUE)
        {
          XCMS_data[j, "comp_result"] <- NA
          XCMS_data[j, "carbon"] <- 0
          XCMS_data[j, "nitrogen"] <- 0
          XCMS_data[j, "isotope_numbers"] <- 0
          XCMS_data[j, "Bin"] <- myBin
          XCMS_data[j, "ID"] <- myID
        }



        #print(XCMS_data[j,])
      }
    }
  }






  #subset by just the column that shrikaar wants, now
  #only labeled data
  XCMS_data = XCMS_data[is.na(XCMS_data$comp_result) == FALSE,]


  #get tge xcms data as input
  XCMS_dataMToBin = XCMS_data

  #get the bins
  allBins = XCMS_dataMToBin$Bin
  allBins = unique(allBins)

  #have a vector to keep track of all
  #the compounds without an MID
  vecOfNoMIDs = vector()

  #collapse all the MIDs into one
  XCMS_dataAvgs = XCMS_data


  #we're going to have a column by which we bin
  #and we don't need a lot of the information in here
  XCMS_dataAvgs = subset(XCMS_dataAvgs,  select = -c(carbon, nitrogen, isotope_numbers, comp_result))

  #only one average per bin
  XCMS_dataAvgs = XCMS_dataAvgs[!duplicated(XCMS_dataAvgs$Bin),]


  #subset for each bin
  for(i in 1:length(allBins))
  {




    myBin = allBins[i]
    print(myBin)
    print("is my bin")
    #subset




    #focus on just the xcms datas associated with a bin
    AllMIDSubsBeforeCategories = subset(XCMS_dataMToBin,  Bin == myBin)

    #determine the formula of whatever is in the bin
    myFormula = unique(AllMIDSubsBeforeCategories$ID)

    #get the carbons an nitrogens in the formula
    #carbonsFormula = get_element_count( myFormula )[['C']]
    #nitrogenFormula = get_element_count( myFormula )[['N']]



    #get all of the isotopologues
    isotopologueList = AllMIDSubsBeforeCategories$comp_result

    #for the denominator of the averages, we're going to want to have
    #all of the possible C's and N's which could be labeled
    #CandNnumbers = carbonsFormula + nitrogenFormula


    #get the mid0 for the bin
    justM0 = subset(AllMIDSubsBeforeCategories,  carbon == 0 & nitrogen == 0, select = -c(Bin,ID, mz, rt, carbon, nitrogen, isotope_numbers, comp_result))

    #get the total number of nitrogens
    nitrogens = AllMIDSubsBeforeCategories$nitrogen
    #get the total number of carbons
    carbons = AllMIDSubsBeforeCategories$carbon

    #now we've got everything that we need
    AllMIDSubsBeforeCategories = subset(XCMS_dataMToBin,  Bin == myBin, select = -c(Bin,ID, mz, rt, carbon, nitrogen, isotope_numbers, comp_result) )




    #justToSum = subset(AllMIDSubs, select = -c(Bin,ID, mz, rt, carbon, nitrogen, isotope_numbers, comp_result))
    #justToSum = subset(AllMIDSubs, select = -c(Bin,ID, mz, rt, carbon, nitrogen, isotope_numbers, comp_result)

    #justToSum = colSums(justToSum)

    #there's no M0, so we're gonna have to remove this
    if(nrow(justM0) == 0)
    {
      vecOfNoMIDs = c(vecOfNoMIDs, myBin)

    }

    #if we have mutiple MIDs, we're gonna start binning them
    if(nrow(justM0) > 0)
    {



      #get the replicates

      #see if we have replicates for multiple categories
      reps = unique(data_cleanEnd(colnames(AllMIDSubsBeforeCategories)))

      #get these replicates
      reps = unique(reps)

      #and get the categories as well, which may be treatment coniditions
      categories = unique(data_cleanII(colnames(AllMIDSubsBeforeCategories)))
      categories = unique(categories)


      #subset each set of mids row by its replicate
      for(rep in 1:length(reps))
      {
        #in each category
        for(category in 1:length(categories))
        {


          justReps = colnames(AllMIDSubsBeforeCategories)[colnames(AllMIDSubsBeforeCategories) %like% paste0("_",reps[rep], sep = NULL)]

          repsAndCategories = justReps[justReps %like% categories[category]]
          print(repsAndCategories)


          whichRep = reps[rep]
          whichCategory = categories[category]
          AllMIDSubs = subset(AllMIDSubsBeforeCategories, select = repsAndCategories)

          #this will take us across each timepoint for the given
          #category and replicate - how does this sample, ID'd by
          #category and replicate change across time as evidenced by the
          #behavior of its MIDs
          for(j in 1:ncol(AllMIDSubs))
          {
            myToGetPercent = AllMIDSubs[,j]
            justToSum = colSums(AllMIDSubs)[j]

            print(myToGetPercent)

            justPercent =  (myToGetPercent / justToSum) * 100
            #print(justPercent)

            print("starting")

            #effectively "refill" the table
            #and also calculate the avgs information

            #store all the avgs information

            #will store all of the information at a single timepoint for all the MIDs
            vectorOfAvgsInfo = vector()
            allMIDsNumber = sum(AllMIDSubsBeforeCategories$isotope_numbers)

            #this will take us across all the isotoplogs for a timepoint
            for(k in 1:length(justPercent))
            {
              nameOfColumn = colnames(AllMIDSubs)[j]
              toReplace = justPercent[k]

              myIsoBin = isotopologueList[k]

              #fill in all the elements in the jth row
              #by each of their columns
              XCMS_dataMToBin[XCMS_dataMToBin$comp_result == myIsoBin, colnames(XCMS_dataMToBin) == nameOfColumn] = toReplace


              #multiple by M[k]
              AvgsInfo = toReplace * carbons[k] + nitrogens[k]

              #sum up all of the isotopologues times their percentages
              vectorOfAvgsInfo = c(vectorOfAvgsInfo, AvgsInfo)

              print(myBin)
              print(nameOfColumn)
              print(toReplace)


              #if(toReplace == 100)
              #{
              #  print(toReplace)
              #  print("is toReplace")
              #}
            }

            #multiple each isotopologue by the number of carbons/nitrogens but all
            #divide everything by the total number of isotopologues

            print(vectorOfAvgsInfo)
            print("is vectorOfAvgsInfo")

            CandNnumbers = sum(carbons) + sum(nitrogens)

            #collapse all the isotopolg info into one datapoint
            vectorOfAvgs = sum(vectorOfAvgsInfo) / CandNnumbers

            if(sum(CandNnumbers) == 0)
            {
              vectorOfAvgs = 0

            }


            #repopulate the XCMS_dataAvgs table with all of the updated data
            XCMS_dataAvgs[XCMS_dataAvgs$Bin == myBin, colnames(XCMS_dataAvgs) == nameOfColumn] = vectorOfAvgs



            print(vectorOfAvgs)
            print("is vectorOfAvgsInfo")
          }
        }
      }
    }
  }


  #remove everything matching the bins to drop, now
  XCMS_dataMToBin = XCMS_dataMToBin[!XCMS_dataMToBin$Bin %in% vecOfNoMIDs,]


  #initialize with one of the other columns
  XCMS_dataMToBin$Isotopologue = XCMS_dataMToBin$isotope_numbers


  for(i in 1:nrow(XCMS_dataMToBin))
  {

    carbonNum = XCMS_dataMToBin[i,]$carbon
    nitrogenNum = XCMS_dataMToBin[i,]$nitrogen


    nameToAdd = paste0("C", carbonNum, "N", nitrogenNum)
    #justCarbonAndNitro = subset(AllMIDSubs, select = c(carbon, nitrogen))

    XCMS_dataMToBin[i, colnames(XCMS_dataMToBin) == "Isotopologue"] = nameToAdd


  }

  #colnames(AllMIDSubs)[colnames(AllMIDSubs) %like% data_clean(colnames(AllMIDSubs))]


  #return the filed
  listOfFiles = list()

  #put out the MIDs
  write.table(XCMS_dataMToBin, file = "row_summed_MIDs_softcodedIII.txt", sep = "\t")



  XCMS_dataAvgs = XCMS_dataAvgs[!XCMS_dataAvgs$Bin %in% vecOfNoMIDs,]

  #put out the avgs
  write.table(XCMS_dataAvgs, file = "xcms_avgs_no_nans.txt", sep = "\t")

  listOfFiles[[1]] = XCMS_dataMToBin
  listOfFiles[[2]] = XCMS_dataAvgs


  return(listOfFiles)
}
