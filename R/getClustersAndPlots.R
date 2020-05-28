
#we'll need the list of outputs
#the name of the category
#the number of clusters
#the kmeans output
getClustersAndPlots <- function(overallToAddPlants,tissueSubset, tissueName, AvgClusNum, cluster_infoAveragesTissue)
{


  #bring the kmeans and the PCA in here





  #have one object to store all of the plots
  plotsList <- list()


  #get the tissue name, pare down these variables, note well
  tissueName = tissueName
  tissue = tissueName

  #how many k-means clusters
  AvgClusNum = AvgClusNum



  #this data frame is going to have the global information for a complete k-mean cluster
  #across times
  dfAllClusters = createDataframe(c("Time","mean","sd","Cluster"))
  dfAllClusters$Time = as.numeric(dfAllClusters$Time)
  dfAllClusters$mean = as.numeric(dfAllClusters$mean)
  dfAllClusters$sd = as.numeric(dfAllClusters$sd)
  dfAllClusters$Cluster = as.numeric(dfAllClusters$Cluster)
  cluster_infoAveragesTissue=cluster_infoAveragesTissue

  #for each of the clusters, we will want to know all of the data
  #to plot
  for(i in 1:AvgClusNum)
  {

    #cluster specific information
    theCluster = i

    #we reuse i before making it up here, so have a placeholder for the cluster number
    #since we will need it after the resetting
    theClusterAvgNum = i

    #plot out everything associated with each cluster
    #this is the kmeans results on the Avgs table
    toSubsetTable = names(cluster_infoAveragesTissue$cluster[which(cluster_infoAveragesTissue$cluster == i)])

    #for each cluster, pull out just those metabolites which fall in that cluster
    for_plottingFurther = subset(tissueSubset, rownames(tissueSubset) %in% toSubsetTable)

    #get their MIDs
    #make sure that the naming is consistent, i.e. if the avgs are named by bin do the same for MIDs
    #and likewise if the avgs are named by compound

    ##beware of the ambiguity regarding whether or not we want to use ID or Bin
    #here as well
    #may want to softcode this, note well


    #for_plottingMIDs <- metabsTableMIDsMerged[ which(metabsTableMIDsMerged$Bin %in% rownames(for_plottingFurther)), ]


    #we have to match up the avgs and the MIDs
    for_plottingMIDs <- metabsTableMIDsMerged[ which(metabsTableMIDsMerged$ID %in% rownames(for_plottingFurther)), ]


    #get the bin info
    justBinInfo = for_plottingMIDs$Bin

    #make sure that it's unique
    justBinInfo = make.names(justBinInfo, unique = TRUE)

    for_plottingMIDs$Bin = NULL

    #for_plottingMIDs <- metabsTableMIDsMerged[ which(metabsTableMIDsMerged$ID %in% rownames(for_plottingFurther)), ]
    #for_plottingMIDs$ID = NULL

    #for_plotting = metabsTableMIDsMerged[metabsTableMIDsMerged$ID %in% toSubsetTable, ]

    #softcode all of this directly from the table
    vecOfExpTimes = unique(data_clean(colnames(for_plottingFurther)))

    ####get all the specific information about each Avg cluster
    for(timepoint in 1:length(vecOfExpTimes))
    {
      whichTimepoint = vecOfExpTimes[timepoint]

      #this is going to be the mean and avg for all of the metabolies associated with a cluster
      #in total, we will want the avg and std at each time point for all of the compounds
      #in that cluster
      myMeanToAdd = mean(rowMeans(for_plottingFurther[,colnames(for_plottingFurther) %like% whichTimepoint]))
      mySdToAdd = sd(rowMeans(for_plottingFurther[,colnames(for_plottingFurther) %like% whichTimepoint]))
      dfAllClusters = rbind(dfAllClusters, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd,Cluster=theClusterAvgNum))
    }


    #get the number of clusters for the MIDs kmeans
    numberOfClusters = as.integer(sqrt(nrow(for_plottingMIDs) / 2))

    #numberOfClusters = calcTheCutOff(for_plottingMIDs)
    #run the function
    #generateInfo(for_plottingCot, "Cotyledon", as.integer(numberOfClusters))

    #aha, make sure that our MIDs are just our tissue of interest
    for_plottingMIDs = for_plottingMIDs[,colnames(for_plottingMIDs) %like% tissue]

    #now, do kmeans on the many MIDs associated
    #with the Avg concentrations that cluster together
    cluster_infoAveragesFurtherClustering = kmeans(for_plottingMIDs, numberOfClusters,nstart = howmanyN)


    allAverageAndMIDSCluster = list()

    if(nrow(for_plottingMIDs) > ncol(for_plottingMIDs))
    {

      clusplotName = paste("clusplot for MIDs k-means on averages cluster ", i)

    }


    #frame to store all the MID-specific information

    #store each MID cluster's information

    #for each of the clusters of MIDs
    #pool together all of the information
    #of the MIDs in that cluster
    dfAllClustersMIDs <- data.frame(
      Time=numeric(),
      mean=numeric(),
      sd=numeric(),
      Cluster=factor()
    )

    #within each cluster of MIDs
    #keep track of all of the abundance data
    dfMIDs <- data.frame(
      Time=numeric(),
      mean=numeric(),
      sd=numeric(),
      Compound=factor(),
      cluster=factor()
    )

    #store out here -- we're going to build up everyhting for all the MIDs clusters
    MidsPlotList = list()


    ##for each MIDs cluster
    #we're going to go through and collect the information to fill in the data frames
    for(j in 1:length(unique(cluster_infoAveragesFurtherClustering$cluster)))
    {

      print(j)
      toSubsetTableII = names(cluster_infoAveragesFurtherClustering$cluster[which(cluster_infoAveragesFurtherClustering$cluster == j)])
      for_plottingSub = subset(metabsTableMIDsMerged, rownames(metabsTableMIDsMerged) %in% toSubsetTableII)
      for_plottingSub= for_plottingSub[,colnames(for_plottingSub) %like% tissue]

      #theMIDs cluster is here
      theCluster = j

      ####get all the specific information about each Avg cluster
      for(timepoint in 1:length(vecOfExpTimes))
      {
        whichTimepoint = vecOfExpTimes[timepoint]
        myMeanToAdd =  mean(rowMeans(for_plottingSub[,colnames(for_plottingSub) %like% whichTimepoint]))
        mySdToAdd =  sd(as.vector(as.matrix(for_plottingSub[,colnames(for_plottingSub) %like% whichTimepoint])))

        #get all the MIDs associated with that
        dfAllClustersMIDs = rbind(dfAllClustersMIDs, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd,Cluster=theCluster))
      }


      #Get the compound specific information in each MIDs cluster
      for(k in 1:length(rownames(for_plottingSub)))
      {
        theCompound = rownames(for_plottingSub)[k]
        for(timepoint in 1:length(vecOfExpTimes))
        {
          whichTimepoint = vecOfExpTimes[timepoint]

          myMeanToAdd =  mean(as.numeric(for_plottingSub[rownames(for_plottingSub) == theCompound,colnames(for_plottingSub) %like% whichTimepoint]))
          mySdToAdd =  sd(as.numeric(for_plottingSub[rownames(for_plottingSub) == theCompound,colnames(for_plottingSub) %like% whichTimepoint]))
          #we're looking at all of the MIDs in a single cluster
          dfMIDs = rbind(dfMIDs, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd, Compound=theCompound, cluster = j))
        }

      }

      theTotalName = paste(tissue, "Avg Cluster", theClusterAvgNum, "MIDs cluster", j)
      whichMIDsCluster = j
      toSubsetPlotMIDs = dfMIDs[dfMIDs$cluster == j,]
      toSubsetPlotMIDs$Time = as.numeric(gsub("X","",as.character(toSubsetPlotMIDs$Time)))

      #we've re-ordered after sorting by the slope

      ###do the slope and sorting now
      #obv we'll want this to be function since it is repeated
      #at the mids and averages level
      #note that we've hardcoded this, for now
      for(l in 1:length(unique(toSubsetPlotMIDs$Compound)))
      {


        #extract all of the MIDs
        myCompound = unique(toSubsetPlotMIDs$Compound)[l]

        #get the slope through the midpoint time
        midTime = unique(toSubsetPlotMIDs)$Time[length(unique(toSubsetPlotMIDs$Time)) / 2]
        startTime = unique(toSubsetPlotMIDs$Time)[1]
        p1 = subset(toSubsetPlotMIDs, subset=(Time=="0" & Compound== myCompound))$mean
        p2 = subset(toSubsetPlotMIDs, subset=(Time==midTime & Compound== myCompound))$mean


        theSlope = p2 - p1
        print(theSlope)

        #we'll want to change this to ClusterSlope, actually
        toSubsetPlotMIDs[toSubsetPlotMIDs$Compound ==  myCompound, "ClusterSlope"]  = rep(theSlope,length(toSubsetPlotMIDs$Compound[toSubsetPlotMIDs$Compound == myCompound]))

      }

      ##sort by the slope
      toSubsetPlotMIDs = toSubsetPlotMIDs[order(toSubsetPlotMIDs$ClusterSlope,decreasing=T),]

      labeledSortVec = vector()
      for(i in 1:length(unique(toSubsetPlotMIDs$Compound)))
      {
        #this will ultimately need to be softcoded, but should
        #be okay for now
        labeledSortVec=  c(labeledSortVec,rep(i,length(unique(vecOfExpTimes))))
      }

      toSubsetPlotMIDs$slopeSorted = labeledSortVec

      #set the factor levels based off of the reordered by slope, now
      toSubsetPlotMIDs$SortNames = toSubsetPlotMIDs$Compound %>% factor(levels = unique(toSubsetPlotMIDs$Compound))

      #plot all of the MIDs from this cluster of MIDs associated with the Avgs Cluster
      MidsPlotList[[j]] =  ggplot(toSubsetPlotMIDs, aes(x=Time,y=mean,colour=SortNames,group=SortNames)) + geom_line() + ggtitle(theTotalName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(SortNames)), show.legend = F, alpha=.1)

      #ggplot(toSubsetPlotMIDs, aes(x=Time,y=mean,colour=SortNames,group=SortNames)) + geom_line() + ggtitle(theTotalName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(Compound)), show.legend = F, alpha=.1)

      #leave out +  ylim(0,100)

      #ggplot(toSubsetPlotMIDs, aes(x=Time,y=mean,colour=Compound,group=Compound)) + geom_line() + ylim(0,100) + ggtitle(theTotalName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(Compound)), show.legend = F, alpha=.1)

    }


    dfAllClustersMIDs$Cluster = as.factor(dfAllClustersMIDs$Cluster)
    MIDsAllClustersPlotName = paste(" Plot for Across all MIDs clusters At Each Time in Avg Cluster", theClusterAvgNum)

    #plot everything for the MIDs cluster, now
    dfAllClustersMIDs$Time = as.numeric(gsub("X","",as.character(dfAllClustersMIDs$Time)))


    ##do the ranking here, as well
    #get the slope of the total cluster.  Although we will want to condense all of this
    #ultimately, this is the beginning
    #note that we've hardcoded this, for now

    #note that we've hardcoded this, for now
    #we are going to want to softcode this
    #because it will be important that we can do this no matter how many different timepoints
    #there are


    #for each cluster of MIDs, calculate the collective slope
    for(j in 1:length(unique(dfAllClustersMIDs$Cluster)))
    {

      #lm(y ~ x)$coeff[[2]]
      p1 = subset(dfAllClustersMIDs, subset=(Time==unique(dfAllClustersMIDs$Time)[1] & Cluster==j))$mean
      p2 = subset(dfAllClustersMIDs, subset=(Time==unique(dfAllClustersMIDs$Time)[(length(unique(dfAllClustersMIDs$Time)) / 2)] & Cluster==j))$mean

      theSlope = p2 - p1
      print(theSlope)

      dfAllClustersMIDs[dfAllClustersMIDs$Cluster == j, "ClusterSlope"]  = rep(theSlope,length(dfAllClustersMIDs$Cluster[dfAllClustersMIDs$Cluster == j]))

      #get the slop of the first to numerbs
      #lm(dfAllClusters[dfAllClusters$Cluster == j,]$mean ~ dfAllClusters[dfAllClusters$Cluster == j,]$Time)$coeff[[2]]


    }


    dfAllClustersMIDs = dfAllClustersMIDs[order(dfAllClustersMIDs$ClusterSlope,decreasing=T),]


    labeledSortVec = vector()

    #it will be each cluster which contains the information
    #of many compounds
    for(i in 1:length(unique(dfAllClustersMIDs$Cluster)))
    {
      labeledSortVec=  c(labeledSortVec,rep(i,  length(unique(vecOfExpTimes))))
    }

    dfAllClustersMIDs$slopeSorted = labeledSortVec

    #make sure that we have the legends sorted here, now
    dfAllClustersMIDs$SortNames = dfAllClustersMIDs$Cluster %>% factor(levels = unique(dfAllClustersMIDs$Cluster))

    listToAdd1 = list()
    listToAdd1[[1]] = ggplot(dfAllClustersMIDs, aes(x=Time,y=mean,colour=SortNames,group=SortNames)) + geom_line() + ggtitle(MIDsAllClustersPlotName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(SortNames)), show.legend = F, alpha=.1)

    #ggplot(dfAllClusters, aes(x=Time,y=mean,colour=Cluster,group=as.factor(Cluster))) + geom_line() + ggtitle(forAllAvgsCluster) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(Cluster)), show.legend = F, alpha=.1) + geom_dl(aes(label = slopeSorted), method = list(dl.combine("top.points"), cex = 1))
    #store the average information before the MIDs for each cluster
    #kinetic information on each MIDs cluster
    allAverageAndMIDSCluster = prepend(MidsPlotList,listToAdd1)
    plotsList = prepend(allAverageAndMIDSCluster,plotsList)
  }


  #we will have to do the labeling here, note well

  #determine the cluster by slope
  #working with the averages
  dfAllClusters$ClusterSlope = dfAllClusters$Cluster
  dfAllClusters$Cluster = as.factor(dfAllClusters$Cluster)
  #determine which cluster has the best slope

  #subset by each cluster
  #note that we've hardcoded this, for now
  for(j in 1:length(unique(dfAllClusters$Cluster)))
  {

    #lm(y ~ x)$coeff[[2]]
    p1 = 0
    p2 = 0

    p1 = subset(dfAllClusters, subset=(Time==unique(dfAllClusters$Time)[1] & Cluster==j))$mean
    p2 = subset(dfAllClusters, subset=(Time==unique(dfAllClusters$Time)[(length(unique(dfAllClusters$Time)) / 2)] & Cluster==j))$mean



    theSlope = 0
    theSlope = as.numeric(p2) - as.numeric(p1)
    print(theSlope)

    dfAllClusters[dfAllClusters$Cluster == j, "ClusterSlope"]  =  theSlope

    #This was just a place holder
    #rep(theSlope,length(dfAllClusters$Cluster[dfAllClusters$Cluster == j]))

    #get the slop of the first to numerbs
    #lm(dfAllClusters[dfAllClusters$Cluster == j,]$mean ~ dfAllClusters[dfAllClusters$Cluster == j,]$Time)$coeff[[2]]


  }

  forAllAvgsCluster = paste(tissueName, ": K-means Clustering")


  #dfAllClusters = dfAllClusters[order(dfAllClusters$ClusterSlope,decreasing=T),]


  ##for now, we want the ordering to be consistent with the kmeans biplot
  dfAllClusters = dfAllClusters[order(dfAllClusters$Cluster,decreasing=F),]

  #we've re-ordered after sorting by the slope

  #this is going to need to be softcoded as well.  instead of being '6' we will change
  #this to the length of timepoints

  ###rep(mytestrep,  length(unique(vecOfExpTimes)) )
  #this is going to be our approach instead
  allTimesVector = vector()
  for(time in 1:AvgClusNum)
  {
    allTimesVector  = c(allTimesVector,rep(time,  length(unique(vecOfExpTimes))))
  }
  dfAllClusters$slopeSorted = allTimesVector


  #allColorsVector = vector()
  #for(time in 1:AvgClusNum)
  #{
  #  color = myColors[time]
  #  allColorsVector  = c(allColorsVector,rep(color,  length(unique(vecOfExpTimes))))
  #}


  #dfAllClusters$slopeSorted = c(rep(1,  length(unique(vecOfExpTimes)) ),rep(2,  length(unique(vecOfExpTimes)) ), rep(3,  length(unique(vecOfExpTimes))), rep(4,  length(unique(vecOfExpTimes))), rep(5,  length(unique(vecOfExpTimes))))
  dfAllClusters$slopeSorted = as.factor(dfAllClusters$slopeSorted)
  dfAllClusters$Cluster = as.factor(dfAllClusters$Cluster)
  dfAllClusters$Time = as.numeric(gsub("X","",as.character(dfAllClusters$Time)))


  #note that all we need to do is change this from slopeSorted
  allClusterForBeginning = list()


  #plotInfoColors
  allClusterForBeginning[[1]] = ggplot(dfAllClusters, aes(x=Time,y=mean,colour=Cluster,group=Cluster)) + geom_line() + ggtitle(forAllAvgsCluster) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(Cluster)), show.legend = F, alpha=.1)

  #don't sort now  ggplot(dfAllClusters, aes(x=Time,y=mean,colour=slopeSorted,group=slopeSorted)) + geom_line() + ggtitle(forAllAvgsCluster) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(slopeSorted)), show.legend = F, alpha=.1)


  #detemine its relative rank

  #dfAllClusters[order(dfAllClusters$ClusterSlope),]

  ##prepend all averages information this to the front
  allTotClusters = prepend(allClusterForBeginning,plotsList)


  ##add the overall kmeans plot of the averages
  allTotClusters = prepend(allTotClusters,allClusterForBeginning)

  ##Now, loop through and plot !!
  #pdf("test_loop_to_print_final_feb18.pdf")
  #clusplot(tissueSubset, cluster_infoAveragesTissue$cluster, main='2D representation of the Cluster solution',
  #        color=TRUE, shade=TRUE,
  #         labels=3, lines=0)


  heatmapName = paste("Heatmap for :", tissue, sep = " ")
  #print(heatmap(as.matrix(dat), Colv = NA, scale = c("none"), col = cm.colors(300), main = heatmapName))


  #add the dendrogram here


  #for(i in allTotClusters){print(i)}
  #dev.off()

  return(allTotClusters)
}
