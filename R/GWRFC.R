#'@title Geographically weighted Random Forest Classification (GWRFC)
#'@description CHECK...
#'@param input_shapefile string or Spatial-class. Input shapefile with dependent and independent variables.  It can be the filename of the shapefile or an object of class SpatialPolygonsDataFrame or SpatialPointsDataFrame.
#'@param remove_columns string. Remove specific variables from input_shapefile. Variables are identified by column name. NA ignores column remove.
#'@param dependent_varName string. Dependent variable name. Must exists at input_shapefile and should be categorical (with not more than 20 classes).
#'@param kernel_type string. Kernel type to apply in GWRFC. It can be: 'gaussian', 'exponential', 'bisquare' or 'tricube'.
#'@param kernel_adaptative logical. Is the kernel adaptative? otherwise it is considered as fixed.
#'@param kernel_bandwidth numeric. Defines kernel bandwidth. If kernel_adaptative is TRUE, then you should define the number of local observations in the kernel, otherwise define you should define a distance to specify kernel bandwidth.
#'@param number_cores numeric. Number of cores for parallel processing. Cores are register and operated via doParallel, foreach and parallel packages. Be careful with increasing numbers of cores, as RAM memory may be not enough.
#'@param output_folder string. Output folder where GWRFC outputs will be stored.
#'@export

GWRFC <- function(
  input_shapefile,
  remove_columns = NA,
  dependent_varName,
  kernel_type = "exponential",
  kernel_adaptative = T,
  kernel_bandwidth,
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
                  "rgeos","scales","doParallel","spgwr","parallel","plyr",
                  "spdep","reshape","rgdal","gtools"))

  ##### DEBUGGING #####

  debug <- F
  if(debug){
    input_shapefile = "C:/DATA/demo/educacion/cant_bach.shp"
    remove_columns = c("ids", "DPA_V", "DPA_A", "DPA_C", "DPA_DESC", "DPA_P", "DPA_DESP",
                       "CODIG", "INSTI", "REGIO", "REGIM", "ZONA_ DISTR", "CIRCU", "PROVI", "CANTO",
                       "PARRO","ZONA_","DISTR")
    dependent_varName = "CRITE"
    kernel_type = "exponential"
    kernel_adaptative = T
    kernel_bandwidth = 50
    number_cores = 3
    output_folder = "C:/DATA/demo/educacion/resultados"
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
  if(class(input_shapefile)=="SpatialPolygonsDataFrame"|class(input_shapefile)=="SpatialPointsDataFrame"){
    model.shp <- input_shapefile
  }else{
    model.shp <- shapefile(input_shapefile)
  }
  #assign rownames
  rownames(model.shp@data) <- 1:length(model.shp)
  #remove columns?
  if(!is.na(remove_columns)[1]){
    if(length(grep(paste(remove_columns,collapse="|"),names(model.shp))) != 0){
      model.shp <- model.shp[,!names(model.shp) %in% remove_columns]
    }else{
      stop("remove_columns not found in input_shapefile")
    }
  }
  #test + get dependent column
  model.dep <- grep(dependent_varName,names(model.shp))
  if(length(model.dep)==0){
    stop("dependent_varName not found")
  }else if(length(model.dep)>=2){
    stop("Found two or more columnsnames in input_shapefile for the specified dependent_varName")
  }else if(length(unique(model.shp@data[,model.dep]))>=21){
    stop(paste0("dependent_varName has ",length(unique(model.shp@data[,model.dep])),
                " classes. Consider to reduce it into 20 classes"))
  }
  #get independent columns
  model.ind <- names(model.shp)[!grepl(dependent_varName,names(model.shp))]
  model.ind <- grep(paste(model.ind,collapse="|"),names(model.shp))
  #prepare data + put as factor dep
  model.shp <- model.shp[,c(model.dep,model.ind)]
  model.shp@data[,1] <- factor(model.shp@data[,1])
  #complete cases
  pos.NA <- which(!complete.cases(model.shp@data))
  if(length(pos.NA)!=0){
    model.shp <- model.shp[which(complete.cases(model.shp@data)),]
    warning(paste0("input_shapefile has ",length(pos.NA)," incomplete case(s). Removing it/them..."))
  }
  #distance matrix
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
    cell.names <- sort(names(model.shp@data[,2:ncol(model.shp@data)]))
    probs.class <- paste0("P_",levels(model.shp@data[,1]))
    cell.names <- c(cell.names,"BEST","DEP","PRED",probs.class,"FAIL","KAPPA","BW")
    cell.out <- as.data.frame(matrix(nrow=1,ncol=length(cell.names)))
    names(cell.out) <- cell.names
    return(cell.out)
  }
  apply.adaptative <- function(){
    #first run
    cell.data <- model.shp@data
    cell.data$dist <- dmat[i,]
    cell.data <- cell.data[order(cell.data$dist)[1:kernel_bandwidth],]
    num.classes <- table(cell.data[,1])
    #expand kernel
    kernel.expand <- kernel_bandwidth
    while(any(num.classes<=4)){
      cell.data <- model.shp@data
      cell.data$dist <- dmat[i,]
      kernel.expand <- kernel.expand + (kernel_bandwidth/10)
      cell.data <- cell.data[order(cell.data$dist)[1:kernel.expand],]
      num.classes <- table(cell.data[,1])
    }
    return(list(cell.data,kernel.expand))
  }
  apply.fixed <- function(){
    #first run
    cell.data <- model.shp@data
    cell.data$dist <- dmat[i,]
    cell.data <- cell.data[cell.data$dist < kernel_bandwidth,]
    num.classes <- table(cell.data[,1])
    #expand kernel
    kernel.expand <- kernel_bandwidth
    while(any(num.classes<=4)){
      cell.data <- model.shp@data
      cell.data$dist <- dmat[i,]
      kernel.expand <- kernel.expand + (kernel_bandwidth/10)
      cell.data <- cell.data[cell.data$dist < kernel.expand,]
      num.classes <- table(cell.data[,1])
    }
    return(list(cell.data,kernel.expand))
  }
  get.probs <- function(){
    x <- as.data.frame(predict(cell.rf,cell.i[,2:(ncol(cell.i)-1)])$predictions)
    x.pred <- names(x[,which.max(x),drop=F])
    x.prob <- x[,grep(cell.out$DEP,names(x))]
    return(list(x.pred,x.prob))
  }
  get.probsClass <- function(){
    x <- as.data.frame(predict(cell.rf,cell.i[,2:(ncol(cell.i)-1)])$predictions)
    x.pred <- names(x[,which.max(x),drop=F])
    return(list(x.pred,x))
  }
  get.kappa <- function(){
    predictions <- factor(apply(cell.rf$predictions,1,function(x){
      x <- names(x)[which.max(x)]
    }),levels=levels(cell.data[,1]))
    conf.m <- caret::confusionMatrix(table(cell.data[,1],predictions))
    return(conf.m$overall[2])
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
      ker.bw <- cell.data[[2]]
      cell.data <- cell.data[[1]]
    }else{
      cell.data <- apply.fixed()
      ker.bw <- cell.data[[2]]
      cell.data <- cell.data[[1]]
    }
    #get 'i' observation
    cell.i <- cell.data[1,]
    #CORRUPTED CASE 1: only one observation
    if(nrow(cell.data)==1){
      cell.out <- corrupted.cases()
    }else{
      #drop unused levels
      cell.data[,1] <- droplevels(cell.data[,1])
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
                        probability=T,
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
                            probability=T,
                            importance="permutation")
        }
        #extract importance
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
        cell.out$BEST <- head(names(cell.out[order(cell.out,decreasing=T)]),3)[1]
        #get dependent
        cell.out$DEP <- as.character(cell.i[,1])
        #get prediction
        cell.out$PRED <- get.probsClass()[[1]]
        #get probabilities
        cell.out.prob <- get.probsClass()[[2]]
        names(cell.out.prob) <- paste0("P_", names(cell.out.prob))
        cell.out <- cbind(cell.out,as.data.frame(cell.out.prob))
        #get accuracies metrics
        cell.out$FAIL <- ifelse(cell.out$DEP!=cell.out$PRED,"yes","no")
        cell.out$KAPPA <- get.kappa()
        #get bandwidth applied
        cell.out$BW <- ker.bw
      }
    }
    #status text
    cat(paste0(as.character(i)," of ",nrow(model.shp@data)," \n"),
        file=paste0(output_folder,"/progress.txt"), append=TRUE)
    #end
    return(cell.out)
  }
  stopCluster(cl)
  #extract data results + reoder
  gwc.data <- do.call("rbind.data.frame",gwc.extract)
  #add original rownames
  gwc.data$ID_row <- rownames(model.shp@data)

  #### SAVE SHAPEFILE ####

  print("Start saving...")

  #save shapefile
  model.shp@data <- gwc.data
  output.name <- paste0(output_folder,"/GWRFC_",
                        ifelse(kernel_adaptative,"ADP_","FIX_"),
                        kernel_bandwidth,"_",
                        kernel_type,".shp")
  shapefile(model.shp,output.name,overwrite=T)
  #error warning
  if(length(which(is.na(model.shp@data$PRED)))!=0){
    warning(paste0(length(which(is.na(model.shp@data$PRED)))," observations were not possible to evaluate during Random Forest execution."))
  }
  #end
  print(paste0("check file: ",basename(output.name)))
  print("****GWRFC end sucessfully*****")
}
