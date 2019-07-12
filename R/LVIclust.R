#'@title Local variables importance (LVI) clustering from GWRFC outputs
#'@description CHECK...
#'@param input_shapefile string. Input shapefile with dependent and independent variables. It can be polygons or points based.
#'@param input_GWRFC string. Input shapefile of GWRFC outputs.
#'@param method_hc string. A method to uso for hierarchical clustering with hclust. It can be: "ward.D", "ward.D2","single","complete","average","mcquitty","median","centroid"
#'@param num_clusters numeric. Number of clusters for summarize local variables importance (LVI). If it is defined as 'auto' (default), it is calculated via Calinski-Harabasz Index and the second peak from 2,3,4...20 cluster computations. If numeric, it cannot be less than 2.
#'@param plots logical. If true, plots and summary reports for each cluster are created.
#'@param output_folder string. Output folder where GWRFC outputs will be stored.
#'@export

LVIclust <- function(
  input_shapefile,
  input_GWRFC,
  method_hc="ward.D2",
  num_clusters="auto",
  plots=T,
  output_folder
){

  ##### LIBRARIES #####

  get.libraries <- function(libraries){
    for (i in libraries){
      if(!require(i,character.only=T)){
        install.packages(i)
        library(i,character.only=T)}
    }
  }

  get.libraries(c("raster","stringr","zoo","ggplot2","rgeos","scales",
                  "plyr","reshape","fpc","pracma","rgdal","gtools","foreign"))

  ##### DEBUGGING #####

  debug <- F
  if(debug){
    input_shapefile = "C:/DATA/poli/GWRFC/shp/Puntos_2001.shp"
    input_GWRFC = "C:/DATA/poli/GWRFC/corrida/3_test/GWRFC_ADP_100_exponential.shp"
    num_clusters = "auto"
    plots=T
    method_hc="ward.D2"
    output_folder = "C:/DATA/poli/GWRFC/corrida/3_test"
  }

  ##### PREPARE DATA #####

  print("Reading data...")

  #random + test
  set.seed(666)
  if(num_clusters==1){
    stop("LVI Clustering can´t be less than 2")
  }
  #folder
  dir.create(output_folder,showWarnings = F, recursive = T)
  #read GWRFC + extract LVI data
  gwrfc.shp <- shapefile(input_GWRFC)
  gwrfc.dep <- gwrfc.shp@data[,which(names(gwrfc.shp@data) %in% c("BEST","DEP","PRED","PROB","FAIL","KAPPA","BW"))]
  gwrfc.shp@data <- gwrfc.shp@data[,which(!names(gwrfc.shp@data) %in% c("BEST","DEP","PRED","PROB","FAIL","KAPPA","BW"))]
  gwrfc.na <- complete.cases(gwrfc.shp@data)
  #read Raw data + extract variables
  raw.data <- read.dbf(gsub(".shp$",".dbf",input_shapefile))
  raw.data <- raw.data[,which(names(raw.data) %in% names(gwrfc.shp@data))]

  #### FUNCTIONS ####

  radar.lvi <- function(x){
    #fix polar coord
    cp <- coord_polar(theta = "y")
    cp$is_free <- function() TRUE
    #plot
    p <- ggplot(x, aes(x=value,y=variable,
                       group=CLUSTER,
                       fill=CLUSTER)) +
      theme_bw() +
      theme(plot.title = element_text(size = 15),
            #panel
            panel.grid.major = element_line(size = 0.5, linetype = 'dotted',colour = "black"),
            #legend
            legend.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.key = element_rect(colour = "black"),
            legend.position="bottom",
            #axis x
            axis.title.x=element_text(size=15),
            axis.text.x=element_text(hjust = 1,size=15),
            #axis y
            axis.title.y=element_text(size=15),
            axis.text.y=element_text(hjust = 1,size=15),
            strip.text.y=element_text(size=15),
            #facet titles off
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            #margins
            plot.margin=unit(c(0,0,0,0),"cm")) +
      #legend columns
      guides(col = guide_legend(nrow=2,byrow=T)) +
      guides(fill = "none",linetype="none") +
      #shapes
      geom_polygon(size=0.5) +
      #other
      labs(y="Variables",x="Average importance [%]") +
      cp +
      facet_wrap(~CLUSTER)
    return(p)
  }

  boxplot.raw <- function(x,y_range=c(-5,5)){
    p <- ggplot(x,aes(x=reorder(variable,value), y=value,
                      color=CLUSTER,
                      group=interaction(variable,CLUSTER))) +
      theme_bw() +
      theme(plot.title = element_text(size = 15),
            legend.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.key = element_rect(colour = "black"),
            legend.position="bottom",
            axis.title.x=element_text(size=15),
            axis.text.x=element_text(hjust = 1,size=15,angle=45),
            strip.text.x = element_blank(),
            axis.title.y=element_text(size=15),
            axis.text.y=element_text(hjust = 1,size=15),
            strip.text.y=element_text(size=15),
            strip.background = element_blank(),
            panel.grid.minor = element_blank()) +
      #boxplot geom
      geom_boxplot(position=position_dodge(width = 0.75),
                   show.legend = T) +
      #others
      geom_hline(yintercept=0,linetype="dotted",color="black",size=0.75) +
      coord_cartesian(ylim=y_range) +
      labs(x="Variables",y="Z-scores")
    return(p)
  }

  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  x.summary <- function(x,fun){
    if(fun=="Mode"){
      x.agr <- aggregate(x[,1:(ncol(x)-1),drop=F],
                         by = list(x$CLUSTER),
                         FUN = Mode)
    }else{
      x.agr <- aggregate(x[,1:(ncol(x)-1),drop=F],
                         by = list(x$CLUSTER),
                         FUN = fun)
    }
    names(x.agr)[1] <- "CLUSTER"
    x.agr <- melt(x.agr,id.vars="CLUSTER")
    names(x.agr)[3] <- fun
    return(x.agr)
  }

  get.elbow <- function(x, y, threshold) {
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative
    indices <- which(abs(d2) > threshold)
    return(indices)
  }

  #### CLUSTERING ####

  print("Start clustering...")

  #apply clustering
  gwrfc.clus <- hclust(dist(gwrfc.shp@data[gwrfc.na,]), method = method_hc)
  #get recommended clusters
  if(num_clusters=="auto"){
    cal.vals <- list()
    for(i in 2:20){
      cal.vals[[i]] <- calinhara(gwrfc.shp@data[gwrfc.na,],cutree(gwrfc.clus, k = i))
    }
    num_clusters <- get.elbow(unlist(cal.vals),1:19,2)[1]
  }
  #add data to LVI
  gwrfc.shp@data$CLUSTER <- NA
  gwrfc.shp@data[gwrfc.na,]$CLUSTER <- cutree(gwrfc.clus, k = num_clusters)
  lvi.data <- gwrfc.shp@data
  lvi.data <- lvi.data[gwrfc.na,]
  lvi.data$CLUSTER <- factor(lvi.data$CLUSTER)
  gwrfc.shp@data <- data.frame(CLUSTER=gwrfc.shp@data$CLUSTER)
  #output name
  output.name <- paste0(output_folder,"/LVI_",num_clusters,"clus.shp")
  shapefile(gwrfc.shp,output.name,overwrite=T)

  #### PLOTS ####

  if(plots){

    print("Making plots and reports...")

    #add data to RAW
    raw.data$DEP <- factor(gwrfc.dep$DEP)
    raw.data$CLUSTER <- NA
    raw.data[gwrfc.na,]$CLUSTER <- cutree(gwrfc.clus, k = num_clusters)
    raw.data <- raw.data[gwrfc.na,]
    raw.data$CLUSTER <- factor(raw.data$CLUSTER)

    #### LVI DATA REPORT ####

    #lvi values by cluster
    lvi.report <- x.summary(lvi.data,"mean")
    names(lvi.report)[3] <- "mean"
    lvi.report$sd <- x.summary(lvi.data,"sd")[,3]
    lvi.report$min <- x.summary(lvi.data,"min")[,3]
    lvi.report$max <- x.summary(lvi.data,"max")[,3]
    #radar plot
    lvi.plot <- melt(lvi.data,id.vars="CLUSTER")
    lvi.plot <- radar.lvi(lvi.plot)
    #save
    output.name <- paste0(output_folder,"/LVI_",num_clusters,"clus_radarPlot.jpg")
    ggsave(filename=output.name,lvi.plot,dpi = 300, width=10*num_clusters,height=10*num_clusters,units="cm")
    rm(lvi.plot,lvi.data)

    #### VARIABLES REPORT ####

    #raw data summary
    raw.data <- raw.data[gwrfc.na,]
    raw.type <- sapply(raw.data,class)
    #get summaries of quantitative variables
    quanti.df <- raw.data[,raw.type %in% "numeric",drop=F]
    quanti.df$CLUSTER <- raw.data$CLUSTER
    quanti.report <- x.summary(quanti.df,"mean")
    quanti.report$sd <- x.summary(quanti.df,"sd")[,3]
    quanti.report$min <- x.summary(quanti.df,"min")[,3]
    quanti.report$max <- x.summary(quanti.df,"max")[,3]
    #boxplot
    quanti.df[,1:(ncol(quanti.df)-1)] <- scale(quanti.df[,1:(ncol(quanti.df)-1)])
    quanti.plot <- melt(quanti.df,id.vars="CLUSTER")
    quanti.plot <- boxplot.raw(quanti.plot)
    output.name <- paste0(output_folder,"/LVI_",num_clusters,"clus_varsBoxplot.jpg")
    ggsave(filename=output.name,quanti.plot,dpi = 300, width=10*num_clusters,height=20,units="cm")
    #get summaries of qualitative variables
    quali.df <- raw.data[,raw.type %in% "factor",drop=F]
    quali.report <- x.summary(quali.df,"Mode")

    #### CLUSTER REPORT ####

    clus.data <- gwrfc.dep[gwrfc.na,]
    clus.data$CLUSTER <- NA
    clus.data[gwrfc.na,]$CLUSTER <- raw.data$CLUSTER
    clus.report <- x.summary(clus.data[,c(4,6,7,8)],"mean")
    clus.report$sd <- x.summary(clus.data[,c(4,6,7,8)],"sd")[,3]
    clus.report$min <- x.summary(clus.data[,c(4,6,7,8)],"min")[,3]
    clus.report$max <- x.summary(clus.data[,c(4,6,7,8)],"max")[,3]

    #### DEPENDENT REPORT ####

    dep.report <- table(clus.data$CLUSTER,clus.data$DEP)

    #### SAVE REPORTS ####

    #save reports
    report.data <- list(LVI=lvi.report,
                        VARS_quanti=quanti.report,
                        VARS_quali=quali.report,
                        CLUSTERS=clus.report,
                        DEPENDENT=dep.report)
    output.name <- paste0(output_folder,"/LVI_",num_clusters,"clus_report.rds")
    saveRDS(report.data,output.name)
    return(report.data)
  }

  print("****LVIclust end sucessfully*****")

}




