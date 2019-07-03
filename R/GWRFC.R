#'@title Geographically weighted Random Forest Classification (GWRFC)
#'@description CHECK...
#'@param input_shapefile spatialpolygonsdataframe. Input shapefile with dependent and independent variables
#'@param remove_columns string. Remove specific variables from 'input_shapefile'. Variables are identified by column name. NA ignores column remove.
#'@param dependent_varName string. Dependent variable name. Must exists at 'input_shapefile'.
#'@param kernel_type string. Kernel type to apply in GWRFC. It can be: 'gaussian', 'exponential', 'bisquare' or 'tricube'.
#'@param kernel_adaptative logical. Is the kernel adaptative? otherwise it is considered as fixed.
#'@param kernel_bandwidth numeric. Defines kernel bandwidth. If 'kernel_adaptative' is TRUE, define the number of local observations in the kernel, otherwise define its distance.
#'@param clusters_LVI numeric. Number of clusters for summarize local variables importance (LVI). If it is defined as 'auto' (default), it is calculated via Nbclust package.
#'@param number_cores numeric. Number of cores for parallel processing. Cores are register and operated via doParallel, foreach and parallel packages.
#'@param output_folder string. Output folder where GWRFC results are stored.
#'@export

GWRFC <- function(
  input_shapefile,
  remove_columns = NA,
  dependent_varName,
  kernel_type = "exponential",
  kernel_adaptative = T,
  kernel_bandwidth,
  clusters_LVI = 'auto',
  number_cores = 1,
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

  get.libraries(c("raster","GWmodel","caret","stringr","ranger","zoo","ggplot2",
                  "rgeos","scales","doParallel","NbClust","spgwr","NbClust",
                  "parallel","plyr","spdep","reshape","rgdal","mclust","gtools"))

  ##### DEBUGGING #####

  debug <- F
  if(debug){
    input_shapefile = "C:/DATA/poli/GWRFC/shp/Puntos_2001.shp"
    remove_columns = NA
    dependent_varName = "Class"
    kernel_type = "exponential"
    kernel_adaptative = T
    kernel_bandwidth = 50
    clusters_LVI = "auto"
    number_cores = 3
    output_folder = "C:/DATA/poli/GWRFC/corrida/3_test"
  }

  ##### PREPARE DATA #####

  print("Reading data...")

  #random
  set.seed(666)
  #folder
  dir.create(output_folder,showWarnings = F, recursive = T)
  if(file.exists(paste0(output_folder,"/progress.txt"))){
    unlink(paste0(output_folder,"/progress.txt"),recursive = T, force = T)
  }
  #read shp
  model.shp <- shapefile(input_shapefile)
  #remove columns?
  if(!is.na(remove_columns)[1]){
    if(length(grep(paste(remove_columns,collapse="|"),names(model.shp))) != 0){
      model.shp <- model.shp[,!names(model.shp) %in% remove_columns]
    }else{
      stop("'remove_columns' not found in 'input_shapefile'")
    }
  }
  #get dependent/independent columns
  model.dep <- grep(dependent_varName,names(model.shp))
  if(length(dependent_varName)==0){
    stop("'dependent_varName' not found")
  }
  model.ind <- names(model.shp)[!grepl(dependent_varName,names(model.shp))]
  model.ind <- grep(paste(model.ind,collapse="|"),names(model.shp))
  #prepare data + distance matrix
  model.shp <- model.shp[,c(model.dep,model.ind)]
  model.shp@data[,1] <- factor(model.shp@data[,1])
  dmat <- gw.dist(dp.locat=coordinates(model.shp),rp.locat=coordinates(model.shp))

  ##### FUNCTIONS ####

  zero.variance <- function(x){
    nzv <- nearZeroVar(x[-1], saveMetrics= TRUE)
    nzv.0 <- rownames(nzv)[nzv$zeroVar]
    if(length(nzv.0)!=0){
      x <- x[,!names(x) %in% nzv.0]
    }
    if(length(nzv)!=0){
      x <- x[,!names(x) %in% nzv]
    }
    return(x)
  }
  corrupted.cases <- function(){
    cell.names <- c(cell.vars,"BEST","DEP","PRED","ACC","KAPPA")
    cell.out <- as.data.frame(matrix(nrow=1,ncol=length(cell.vars)+5))
    names(cell.out) <- cell.names
    cell.out$DEP <- as.character(cell.i[,1])
    return(cell.out)
  }
  apply.adaptative <- function(){
    #first run
    cell.data <- model.shp@data
    cell.data$dist <- dmat[i,]
    cell.data <- cell.data[order(cell.data$dist)[1:kernel_bandwidth],]
    num.classes <- length(unique(cell.data[,1]))
    #expand kernel
    kernel.expand <- kernel_bandwidth
    while(num.classes==1){
      cell.data <- model.shp@data
      cell.data$dist <- dmat[i,]
      kernel.expand <- kernel.expand + (kernel_bandwidth/10)
      cell.data <- cell.data[order(cell.data$dist)[1:kernel.expand],]
      num.classes <- length(unique(cell.data[,1]))
    }
    return(cell.data)
  }
  apply.fixed <- function(){
    #first run
    cell.data <- model.shp@data
    cell.data$dist <- dmat[i,]
    cell.data <- cell.data[cell.data$dist < kernel_bandwidth,]
    num.classes <- length(unique(cell.data[,1]))
    #expand kernel
    kernel.expand <- kernel_bandwidth
    while(num.classes==1){
      cell.data <- model.shp@data
      cell.data$dist <- dmat[i,]
      kernel.expand <- kernel.expand + (kernel_bandwidth/10)
      cell.data <- cell.data[cell.data$dist < kernel.expand,]
      num.classes <- length(unique(cell.data[,1]))
    }
    return(cell.data)
  }
  radar.plot <- function(){
    #fix polar coord
    cp <- coord_polar(theta = "y")
    cp$is_free <- function() TRUE
    #plot
    p <- ggplot(gwc.report, aes(x=value,y=variable,
                                linetype=clusters,
                                group=clusters,
                                color=clusters)) +
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
      geom_polygon(fill=NA,size=0.8) +
      #other
      #scale_linetype_manual("",values=c("solid","dashed")) +
      labs(y="Variables",x="Average importance [%]") +
      cp
    return(p)
  }

  #### GW RANDOM FOREST ####

  print("Start processing...")

  #process
  cl <- makeCluster(number_cores)
  registerDoParallel(cl)
  gwc.extract <- foreach(i=1:nrow(model.shp@data),.packages=c("ranger","scales","caret","GWmodel"),.errorhandling="pass") %dopar% {
    #subset by kernel_bandwidth
    if(kernel_adaptative){
      cell.data <- apply.adaptative()
    }else{
      cell.data <- apply.fixed()
    }
    #get 'i' observation + drop unused levels
    cell.i <- cell.data[1,]
    cell.data[,1] <- droplevels(cell.data[,1])
    #CORRUPTED CASE 1: only one observation
    if(nrow(cell.data)==1){
      cell.out <- corrupted.cases()
    }else{
      #balance + set levels in classification
      cell.data <- upSample(cell.data[,2:ncol(cell.data)],cell.data[,1],yname=names(cell.data)[1])
      cell.data <- cell.data[,c(ncol(cell.data),1:ncol(cell.data)-1)]
      #distance weights
      cell.weights <- gw.weight(cell.data$dist,
                                bw=kernel_bandwidth,
                                kernel=kernel_type,
                                adaptive=kernel_adaptative)
      #assign dependent + remove zero variance
      cell.data <- cell.data[,1:(ncol(cell.data)-1)]
      cell.vars <- sort(names(cell.data[,2:ncol(cell.data)]))
      cell.data <- cbind(cell.data[,1],zero.variance(cell.data[,-1]))
      names(cell.data)[1] <- dependent_varName
      #apply ranger: FIRST TIME
      cell.formula <- formula(paste0(dependent_varName," ~ ."))
      cell.rf <- ranger(formula=cell.formula,
                        data=cell.data,
                        replace=F,
                        min.node.size=1,
                        scale.permutation.importance=T,
                        case.weights=cell.weights,
                        importance="permutation")
      #extract important variables
      cell.best <- names(cell.rf$variable.importance)[cell.rf$variable.importance > 0]
      #CORRUPTED CASE 2: null calculation after run
      if(length(cell.best)==0){
        cell.out <- corrupted.cases()
      }else{
        #apply ranger: SECOND TIME
        if(length(cell.best)!=length(cell.vars)){
          cell.formula <- formula(paste0(dependent_varName," ~ ",paste(cell.best,collapse=" + ")))
          cell.rf <- ranger(formula=cell.formula,
                            data=cell.data,
                            replace=F,
                            min.node.size=1,
                            scale.permutation.importance=T,
                            case.weights=cell.weights,
                            importance="permutation")
        }
        #extract importance + set 0 to negative
        cell.out <- data.frame(t(matrix(cell.rf$variable.importance)))
        names(cell.out) <- names(cell.rf$variable.importance)
        cell.out[,cell.out < 0] <- 0
        #add eliminated columns as 0 + sort
        cell.rem <- cell.vars[!cell.vars %in% names(cell.out)]
        if(length(cell.rem!=0)){
          cell.val <- as.data.frame(t(matrix(rep(0,length(cell.rem)))))
          names(cell.val) <- cell.rem
          cell.out <- cbind(cell.out,cell.val)
        }
        cell.out <- cell.out[sort(names(cell.out))]
        #get three best variables
        cell.out$BEST <- paste(head(names(cell.out[order(cell.out,decreasing=T)]),3),collapse="+")
        #get prediction
        cell.out$DEP <- as.character(cell.i[,1])
        cell.out$PRED <- as.character(predict(cell.rf,cell.i[,2:(ncol(cell.i)-1)])$predictions)
        #get confusion matrix metrics
        if(any(is.na(colnames(cell.rf$confusion.matrix)))){
          cell.rf$confusion.matrix <- cell.rf$confusion.matrix[,!is.na(colnames(cell.rf$confusion.matrix))]
        }
        cell.out$ACC <- caret::confusionMatrix(cell.rf$confusion.matrix)$overall[1]
        cell.out$KAPPA <- caret::confusionMatrix(cell.rf$confusion.matrix)$overall[2]
      }
    }
    #status text
    cat(paste0(as.character(i)," of ",nrow(model.shp@data)," \n"),
        file=paste0(output_folder,"/progress.txt"), append=TRUE)
    #end
    return(cell.out)
  }
  stopCluster(cl)
  #extract data results
  gwc.data <- do.call("rbind.data.frame",gwc.extract)

  #### CLUSTERING ####

  print("Start clustering...")

  #get variables & NA not for use in clustering
  index.acc <- (length(gwc.data) - 4):length(gwc.data)
  index.na <- complete.cases(gwc.data)
  #get recommended clusters
  if(clusters_LVI=="auto"){
    gwc.clus <- gwc.data[index.na,1:(index.acc[1]-1)]
    clusters_LVI <- NbClust(gwc.clus,
                            distance = "euclidean",
                            method="kmeans",
                            index="gap")$Best.nc[1]
    gwc.clus <- Mclust(gwc.clus, G=clusters_LVI)
    gwc.data$CLUSTER <- NA
    gwc.data[index.na,]$CLUSTER <- gwc.clus$classification
  }else{
    #based in specific number of clusters
    gwc.clus <- gwc.data[index.na,1:(index.acc[1]-1)]
    gwc.clus <- Mclust(gwc.clus, G=clusters_LVI)
    gwc.data$CLUSTER <- NA
    gwc.data[index.na,]$CLUSTER <- gwc.clus$classification
  }

  #### CREATE REPORT ####

  print("Start report...")

  #from cluster means
  gwc.report <- as.data.frame(t(gwc.clus$parameters$mean))
  gwc.report$clusters <- paste0("CLUS_",1:clusters_LVI)
  gwc.report <- melt(gwc.report,id.vars = "clusters")
  output.name <- paste0(output_folder,"/clusters_VarImp.jpg")
  ggsave(filename=output.name,radar.plot(),dpi = 300, width=32,height=32,units="cm")

  #### SAVE SHAPEFILE ####

  print("Start saving...")

  #save shapefile
  model.shp@data <- gwc.data
  output.name <- paste0(output_folder,"/GWRFC_",
                        ifelse(kernel_adaptative,"ADP_","FIX_"),
                        kernel_bandwidth,"_",
                        kernel_type,".shp")
  shapefile(model.shp,output.name,overwrite=T)
  #warning
  if(length(which(is.na(model.shp@data$PRED)))!=0){
    warning(paste0("kernel_bandwidth to short for some iterations and ",
                   length(which(is.na(model.shp@data$PRED))),
                   " observations were not possible to evaluate."))
  }
  #end
  print(paste0("check file: ",basename(output.name)))
  print("****GWRFC end sucessfully*****")
}
