#'@title Local variables importance (LVI) clustering from GWRFC outputs
#'@description CHECK...
#'@param input_shapefile string or Spatial-class. Input shapefile with dependent and independent variables. It can be the filename of the shapefile or an object of class SpatialPolygonsDataFrame or SpatialPointsDataFrame.
#'@param input_GWRFC string or Spatial-class. Input shapefile of GWRFC outputs. It can be the filename of the shapefile or an object of class SpatialPolygonsDataFrame or SpatialPointsDataFrame.
#'@param method_hc string. A method to use for hierarchical clustering with hclust. It can be:"ward.D","ward.D2","single","complete","average","mcquitty","median" or "centroid".
#'@param clus_data string. Data which should be used in clustering. It can be "LVI" to refer to all LVI variables or any specific column at input_GWRFC. In this case, hierarchical clustering is not applied and target variable is reclassified into quantiles if it is numerical.
#'@param clus_num numeric or string. If applies hierarchical clustering, then it is the number of clusters; otherwise it is the number of quantiles to use. If it is defined as 'auto' (default), then a number is calculated via Calinski-Harabasz Index in hierarchical clustering or defined in 5 for quantiles.
#'@param plots logical. If true, plots are added to clusters report.
#'@param output_folder string. Output folder where LVIclust outputs will be stored.
#'@export

LVIclust <- function(
  input_shapefile,
  input_GWRFC,
  method_hc="ward.D2",
  clus_data="LVI",
  clus_num="auto",
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
    input_shapefile = shapefile("C:/DATA/demo/educacion/cant_bach.shp")
    input_GWRFC = shapefile("C:/DATA/demo/educacion/resultados/GWRFC_ADP_50_exponential.shp")
    clus_data="P_APROB"
    clus_num="auto"
    plots=T
    method_hc="ward.D2"
    output_folder = "C:/DATA/demo/educacion/resultados"
  }

  ##### PREPARE DATA #####

  print("Reading data...")

  #random + test
  set.seed(666)
  if(clus_num==1){
    stop("LVI Clustering can´t be less than 2")
  }
  #folder
  dir.create(output_folder,showWarnings = F, recursive = T)
  #read GWRFC
  if(class(input_GWRFC)=="SpatialPolygonsDataFrame"|class(input_GWRFC)=="SpatialPointsDataFrame"){
    gwrfc.shp <- input_GWRFC
    na.gwrfc <- complete.cases(gwrfc.shp@data)
    gwrfc.shp <- gwrfc.shp[na.gwrfc,]
  }else{
    gwrfc.shp <- shapefile(input_GWRFC)
    na.gwrfc <- complete.cases(gwrfc.shp@data)
    gwrfc.shp <- gwrfc.shp[na.gwrfc,]
  }
  #read RAW + filter rows
  if(class(input_shapefile)=="SpatialPolygonsDataFrame"|class(input_shapefile)=="SpatialPointsDataFrame"){
    raw.data <- input_shapefile@data
  }else{
    raw.data <- read.dbf(gsub(".shp$",".dbf",input_shapefile))
  }
  rownames(raw.data) <- 1:nrow(raw.data)
  raw.data <- raw.data[rownames(raw.data) %in% gwrfc.shp@data$ID_row,]
  gwrfc.shp@data$ID_row <- NULL
  #get LVI & model results + filter columns in RAW
  target.names <- names(gwrfc.shp@data)[grep("P_",names(gwrfc.shp@data))]
  target.names  <- c("BEST","DEP","PRED",target.names,"FAIL","KAPPA","BW")
  gwrfc.dep <- gwrfc.shp@data[,which(names(gwrfc.shp@data) %in% target.names)]
  gwrfc.shp@data <- gwrfc.shp@data[,which(!names(gwrfc.shp@data) %in% target.names)]
  raw.data <- raw.data[,which(names(raw.data) %in% names(gwrfc.shp@data))]

  #### FUNCTIONS ####

  bar.lvi <- function(x){
    #plot
    p <- ggplot(x, aes(x=reorder(variable, value),y=value,
                       group=CLUSTER,
                       fill=CLUSTER)) +
      theme_bw() +
      theme(plot.title = element_text(size = 20),
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
      guides(col = guide_legend(nrow=1,byrow=T)) +
      #shapes
      geom_col(position = "dodge") +
      #other
      coord_flip() +
      facet_wrap(~CLUSTER) +
      labs(title="Average importance of independent variables, grouped by clusters",
           x="Variables",
           y="Average importance (%)")

    return(p)
  }

  boxplot.quanti <- function(x,y_range=c(-3,3)){
    p <- ggplot(x,aes(x=reorder(variable,value), y=value,
                      fill=CLUSTER,
                      group=interaction(variable,CLUSTER))) +
      theme_bw() +
      theme(plot.title = element_text(size = 20),
            #legend
            legend.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.key = element_rect(colour = "black"),
            legend.position="bottom",
            #x axis
            axis.title.x=element_text(size=15),
            axis.text.x=element_text(hjust = 1,size=15,angle=45),
            #y axis
            axis.title.y=element_text(size=15),
            axis.text.y=element_text(hjust = 1,size=15),
            strip.text.y=element_text(size=15),
            #facet titles off
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            #margins
            plot.margin=unit(c(0,0,0,0),"cm")) +
      #boxplot geom
      geom_boxplot(position=position_dodge(width = 0.75),outlier.alpha=0.33,
                   show.legend = T) +
      #others
      geom_hline(yintercept=0,linetype="dotted",color="black",size=0.75) +
      coord_cartesian(ylim=y_range) +
      labs(title="Standarized values for quantitative variables, grouped by clusters",
           x="Variables (Mean ± SD)",
           y="Z-scores")
    return(p)
  }

  hist.quali <- function(x){
    p <- ggplot(x,aes(x=value,
                      fill=CLUSTER)) +
      theme_bw() +
      theme(plot.title = element_text(size = 20),
            #legend
            legend.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.key = element_rect(colour = "black"),
            legend.position="bottom",
            #x axis
            axis.title.x=element_text(size=15),
            axis.text.x=element_text(hjust = 1,size=15,angle=45),
            strip.text.x=element_text(size=15),
            #y axis
            axis.title.y=element_text(size=15),
            axis.text.y=element_text(hjust = 1,size=15),
            strip.text.y=element_text(size=15),
            #margins
            plot.margin=unit(c(0,0,0,0),"cm")) +
      #bar geom
      geom_bar(position = "dodge") +
      #other
      coord_flip() +
      facet_wrap(~variable,scales="free",ncol=4) +
      labs(title="Classes counts for qualitative variables, grouped by clusters",
           x="Count",
           y="Classes")
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

  if(clus_data=="LVI"){
    #get recommended clusters
    gwrfc.clus <- hclust(dist(gwrfc.shp@data), method = method_hc)
    if(clus_num=="auto"){
      cal.vals <- list()
      for(i in 2:20){
        cal.vals[[i]] <- calinhara(gwrfc.shp@data,cutree(gwrfc.clus, k = i))
      }
      clus_num <- get.elbow(unlist(cal.vals),1:19,0.015)[1]
      if(is.na(clus_num)){
        warning("elbow not found. Assumed 2 clusters...")
        clus_num <- 2
        plot(1:19,unlist(cal.vals),xlab="Number of clusters",ylab="Calinhara index")
        abline(v=clus_num,col="red",lty=3)
      }else{
        plot(1:19,unlist(cal.vals),xlab="Number of clusters",ylab="Calinhara index")
        abline(v=clus_num,col="red",lty=3)
      }
    }
    #add data to LVI
    gwrfc.shp@data$CLUSTER <- cutree(gwrfc.clus, k = clus_num)
  }else{
    clus_num <- ifelse(clus_num=="auto",5,clus_num)
    if(length(grep(clus_data,names(gwrfc.dep)))==0){
      stop("clus_data not found")
    }
    gwrfc.shp@data$CLUSTER <- gwrfc.dep[,grep(clus_data,names(gwrfc.dep))]
    if(class(gwrfc.shp@data$CLUSTER)!="factor"|class(gwrfc.shp@data$CLUSTER)!="character"){
      gwrfc.shp@data$CLUSTER <- base::cut(gwrfc.shp@data$CLUSTER,
                                          breaks=quantile(gwrfc.shp@data$CLUSTER,probs=seq(0,1,length.out=clus_num)),
                                          labels=paste0("Q",floor(seq(0,1,length.out=clus_num)*100))[2:clus_num],
                                          include.lowest=T)
    }else{
      gwrfc.shp@data$CLUSTER <- factor(gwrfc.shp@data$CLUSTER)
    }
    clus_num <- nlevels(gwrfc.shp@data$CLUSTER)
  }
  #output file
  cluster.shp <- gwrfc.shp
  cluster.shp@data <- data.frame(CLUSTER=gwrfc.shp@data$CLUSTER)
  output.name <- paste0(output_folder,"/",clus_data,"_",clus_num,"clus.shp")
  shapefile(cluster.shp,output.name,overwrite=T)

  #### LVI DATA REPORT ####

  print("Making report...")

  #get LVI data
  lvi.data <- gwrfc.shp@data
  lvi.data$CLUSTER <- factor(lvi.data$CLUSTER)
  #lvi values by cluster
  lvi.report <- x.summary(lvi.data,"mean")
  names(lvi.report)[3] <- "mean"
  lvi.report$sd <- x.summary(lvi.data,"sd")[,3]
  lvi.report$min <- x.summary(lvi.data,"min")[,3]
  lvi.report$max <- x.summary(lvi.data,"max")[,3]
  #bar plot
  if(plots){
    lvi.plot <- melt(lvi.data,id.vars="CLUSTER")
    lvi.plot <- bar.lvi(lvi.plot)
    output.name <- paste0(output_folder,"/",clus_data,"_",clus_num,"clus_LVIPlot.jpg")
    ggsave(filename=output.name,lvi.plot,dpi = 300, width=10*clus_num,height=10*clus_num,units="cm")
  }

  #### RAW DATA REPORT ####

  #get RAW data
  raw.data$DEP <- factor(gwrfc.dep$DEP)
  raw.data$CLUSTER <- factor(gwrfc.shp@data$CLUSTER)
  raw.type <- sapply(raw.data,class)
  #QUANTI variables by cluster
  quanti.df <- raw.data[,raw.type %in% "numeric",drop=F]
  quanti.df$CLUSTER <- raw.data$CLUSTER
  quanti.report <- x.summary(quanti.df,"mean")
  quanti.report$sd <- x.summary(quanti.df,"sd")[,3]
  quanti.report$min <- x.summary(quanti.df,"min")[,3]
  quanti.report$max <- x.summary(quanti.df,"max")[,3]
  #QUALI variables by cluster
  quali.df <- raw.data[,raw.type %in% c("factor","character"),drop=F]
  quanti.df$CLUSTER <- raw.data$CLUSTER
  quali.report <- x.summary(quali.df,"Mode")
  #plots
  if(plots){
    #for quantitative
    scale.vals <- scale(quanti.df[,1:(ncol(quanti.df)-1)])
    quanti.df[,1:(ncol(quanti.df)-1)] <- scale.vals
    scale.vals <- paste0(names(quanti.df)[1:(ncol(quanti.df)-1)],"\n","(",
                         round(attr(scale.vals,"scaled:center"),2)," ± ",
                         round(attr(scale.vals,"scaled:scale"),2),")")
    names(quanti.df)[1:(ncol(quanti.df)-1)] <- scale.vals
    quanti.plot <- melt(quanti.df,id.vars="CLUSTER")
    quanti.plot <- boxplot.quanti(quanti.plot)
    output.name <- paste0(output_folder,"/",clus_data,"_",clus_num,"clus_quantiPlot.jpg")
    ggsave(filename=output.name,quanti.plot,dpi = 300, width=10*clus_num,height=30,units="cm")
    #for qualitative
    quali.plot <- melt(quali.df,id.vars = "CLUSTER")
    quali.plot <- hist.quali(quali.plot)
    output.name <- paste0(output_folder,"/",clus_data,"_",clus_num,"clus_qualiPlot.jpg")
    set.size <- sum(c(floor((ncol(quali.df)-1)/4),(ncol(quali.df)-1) %% 4))
    ggsave(filename=output.name,quali.plot,dpi = 300, width=40,height=10*set.size,units="cm",limitsize=F)
  }

  #### ACCURACY REPORT ####

  #get accuracy data
  acc.data <- gwrfc.dep
  acc.data$CLUSTER <- raw.data$CLUSTER
  acc.data <- acc.data[,names(acc.data) %in% c("KAPPA","BW","CLUSTER")]
  acc.report <- x.summary(acc.data,"mean")
  acc.report$sd <- x.summary(acc.data,"sd")[,3]
  acc.report$min <- x.summary(acc.data,"min")[,3]
  acc.report$max <- x.summary(acc.data,"max")[,3]

  #### PROBABILITIES REPORT ####

  #get probabilities data
  prob.data <- gwrfc.dep
  prob.data$CLUSTER <- raw.data$CLUSTER
  prob.data <- prob.data[,names(prob.data) %in% c(grep("^P_",names(prob.data),value=T),"CLUSTER")]
  prob.report <- x.summary(prob.data,"mean")
  prob.report$sd <- x.summary(prob.data,"sd")[,3]
  prob.report$min <- x.summary(prob.data,"min")[,3]
  prob.report$max <- x.summary(prob.data,"max")[,3]

  #### SAVE REPORTS ####

  #merge
  report.data <- list(ACCURACY=acc.report,
                      PROBABILITIES=prob.report,
                      QUANTI_VARS=quanti.report,
                      QUALI_VARS=quali.report,
                      LVI=lvi.report)
  #save reports
  output.name <- paste0(output_folder,"/",clus_data,"_",clus_num,"clus_report.rds")
  saveRDS(report.data,output.name)
  return(report.data)
  print("****LVIclust end sucessfully*****")
}




