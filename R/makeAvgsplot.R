##########do the plotting for all of the avgs
makeAvgsplot <- function(mydata, whichTissue, limitPlot)
{


  #read in the data table
  mydata = mydata
  #parameter about to which percent we want to plot
  #for example, if limitPlot= 65 we will only plot up
  #until 65 percent
  limitPlot = limitPlot

  #combine the multiple columns into a single column with unique names
  names = make.names(paste(mydata$Bin,mydata$rt,mydata$m.z,mydata$ID,sep = "_"), unique = TRUE)
  rownames(mydata) = names

  #what is the identifier for this group (i.e. is it a condition or tissue group)
  whichTissue = whichTissue

  #pull out the condition/species specific information
  #by only those columns which match the whichTissue
  mydata =mydata[, !(colnames(mydata) %in% colnames(mydata)[(colnames(mydata) %like% whichTissue)])]

  #read all the data in
  name = paste(whichTissue, "allAvgsPlot.pdf", sep = "_")
  pdf(name)

  par(mfrow=c(3,3))
  for(i in 1:nrow(mydata))
  {

    #reset the plotting frame every 9 plots
    #use modular division to determine if we're at the 9th
    if(i %% 9 == 0)
    {
      par(mfrow=c(3,3))
    }
    theVecOfInfo = mydata[i,]


    #if using id's for title
    #set title before we make it numeric
    title = paste(theVecOfInfo[colnames(theVecOfInfo) == "ID"][1,], theVecOfInfo[colnames(theVecOfInfo) == "Bin"][1,])

    ##calculate each time point and each std dev

    #get the set of times from the column names
    myTimes = data_clean(as.character(colnames(mydata)))
    #pull out any non-numeric elements of the timepoints
    myTimes =  as.numeric(gsub("X","",as.character(myTimes)))
    #remove any non-valid timepoints
    myTimes = myTimes[is.na(myTimes) == FALSE]
    #get a unique vector of the timepoints
    myTimesUnique = unique(myTimes)
    #isotopologueName = as.character(myVecOfAllInfo$Isotopologue)

    #hold all the Avgs and standard deviations
    allMyInfoAvg = vector()
    allMyInfoSD = vector()

    ###get the average and the sd for each of the metabolites

    #for each time
    for(k in 1:length(myTimesUnique))
    {

      #get the time
      myTime = myTimesUnique[k]



      #pull out first field
      myTimepointMatch = gsub("(*.*)_.*_.*", "\\1", names(theVecOfInfo))

      #pull out non-numeric stuff
      myTimepointMatch = gsub("[^0-9.-]", "", myTimepointMatch)

      #rename the colnames with just the time information
      #I think that this is redundant, actually
      names(theVecOfInfo) = gsub("[^0-9.-]", "", names(theVecOfInfo))

      #make sure that it's numeric now
      myTimepointMatch = as.numeric(myTimepointMatch)

      #set the names to just the timepoint now, it's all we need
      names(theVecOfInfo) = as.numeric(myTimepointMatch)

      #calculate the mean for this timepoint
      myTimeMean = mean(as.numeric(theVecOfInfo[,as.numeric(colnames(theVecOfInfo)) %in% c(myTime)]))

      #calculate the SD for this timepoint
      myTimeSD = sd(as.numeric(theVecOfInfo[,as.numeric(colnames(theVecOfInfo)) %in% c(myTime)]))

      allMyInfoAvg = c(allMyInfoAvg, myTimeMean)
      allMyInfoSD = c(allMyInfoSD, myTimeSD)
      #df = rbind(df, data.frame(Isotopologue=isotopologueName, Time=myTime, Mean=myTimeMean, stdDev=myTimeSD))
    }

    #declare the times which will be the x-axis of the plots
    x = myTimesUnique



    #make the plots (time on the x-axis)
    #the averages on the y-axis
    plot(x, allMyInfoAvg,
         ylim=range(c(0, limitPlot)),
         pch=19, xlab="Time", ylab="% 13C enrichment",
         main=title
    )

    #add in the error bars as specified by the clever
    #stack overflow post
    # hack: we draw arrows but with very special "arrowheads"
    arrows(x, allMyInfoAvg-allMyInfoSD, x, allMyInfoAvg+allMyInfoSD, length=0.05, angle=90, code=3)

  }
  dev.off()

}
