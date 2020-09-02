#'@title Geographically weighted Random Forest Classification
#'@description GWRFC is a software for analyze and explore spatial data. It constructs geographically weighted models (GW; Fotheringham et al. 1998) to train random forest (RF; Breiman 2001) and report local models with partial depende plots (PDP, Greenwell, 2019). Prediction results and accurancy metrics (ACC) are also representated accondingly.
#'@param input_shapefile string or Spatial-class. Input shapefile with dependent and independent variables.  It can be the filename of the shapefile or an object of class {SpatialPolygonsDataFrame} or {SpatialPointsDataFrame}.
#'@param remove_columns string. Remove specific variables from \strong{input_shapefile}. Variables are identified by column name. NA ignores column remove.
#'@param dependent_varName string. Dependent variable name. Must exists at \strong{input_shapefile} and should be categorical (with not more than 20 classes).
#'@param kernel_function string. Kernel type to apply in GWRFC. It can be: 'gaussian', 'exponential', 'bisquare' or 'tricube'.
#'@param kernel_adaptative logical. Is the kernel adaptative? otherwise it is considered as fixed (larger processing time).
#'@param kernel_bandwidth numeric. Defines kernel bandwidth. If \strong{kernel_adaptative} is TRUE, then you should define the number of local observations in the kernel, otherwise you should define a distance to specify kernel bandwidth.
#'@param upsampling logical. If TRUE, upsampling is applied before random forest training, otherwise it is downsampled. Consider that upsampling is a bit more computing demanding but accuracy is improved.
#'@param save_models logical. If TRUE, random forest models are stored at \strong{output_folder} as a RDS file. Beware it can be large, therefore storage requires hard drive memory and can slow down algorithm exit.
#'@param enable_pdp logical. --EXPERIMENTAL-- If TRUE, partial dependence plots YHAT maximun, together with its correspondent independent variable value (PDP) are calculated.
#'@param number_cores numeric. Number of cores for parallel processing. Cores are register and operated via doParallel, foreach and parallel packages. Be careful with increasing numbers of cores, as RAM memory may be not enough.
#'@param output_folder string. Output folder where GWRFC outputs will be stored.
#'@return As a result, four shapefiles are created whose prefixes refer to: \enumerate{
#'                            \item LVI: Local variables importance. Calculated via permutation for each variable.
#'                            \item PDP: Independent variables local maxima (class or value). Identified when YHAT reach its maximum during RF model marginalization. Calculated for each variable.
#'                            \item YHAT: Prediction result for \strong{dependent_varName} when PDP local maxima is applied. Calculated for each variable.
#'                            \item ACC: Prediction and accuracies: predicted class, kappa from Out-of-Bag, classes probabilities and prediction failures.
#'                           }
#'        In all shapefiles cases, a column called 'ID_row' refers to rownames of \strong{input_shapefile}. In addition, processing evolution can be monitored at \strong{output_folder} as: data_progress.txt
#'@examples
#'
#'#view deforestation data
#'
#'data("deforestation")
#'tmap_mode("view")
#'tm_basemap("OpenStreetMap") +
#'  tm_shape(deforestation) +
#'  tm_polygons(col="fao",style="cat",title="Annual deforestation rate  2000-2010 (FAO) - categorical (quantiles)",palette="YlOrRd")
#'
#'#run GWRFC
#'
#'GWRFC(input_shapefile = deforestation, #can be a spatial dataframe (points or polygons) or the complete filename of the shapefile to analyze.
#'      remove_columns = c("ID_grid","L_oth"), #for remove variables if they are not informative. Put NA to avoid removal.
#'      dependent_varName = "fao", #the depedent variable to evaluate. It should be of factor or character data type.
#'      kernel_function = "exponential", #the weightening function. See help for other available functions.
#'      kernel_adaptative = T, #use TRUE for adaptative kernel distance or FALSE for a fixed kernel distance.
#'      kernel_bandwidth = 400, #as the kernel is adaptative, 400 refers to the minimun number of observations to use in modelling.
#'      upsampling = T, #improves accuracy (recommended) but is a bit more computing costly.
#'      save_models = T, #save RF models. Beware of hard disk space and extra processing time.
#'      enable_pdp = F, #experimental, use with caution as is sensible to noise.
#'      number_cores = 3, #defines the number of CPU cores to use
#'      output_folder = "E:/demo/deforestation") #check this folder for GWRFC outputs.
#'
#'@export

GWRFC <- function(
  input_shapefile,
  remove_columns = NA,
  dependent_varName,
  kernel_function = "exponential",
  kernel_adaptative = T,
  kernel_bandwidth,
  upsampling = T,
  save_models = F,
  enable_pdp = F,
  number_cores = 1,
  output_folder
){

  ##### PREPARE DATA #####

  print("Reading data...")

  #random
  set.seed(666)
  #folder
  dir.create(output_folder,showWarnings = F, recursive = T)
  if(file.exists(paste0(output_folder,"/data_progress.txt"))){
    unlink(paste0(output_folder,"/data_progress.txt"),recursive = T, force = T)
  }
  #read shp
  if(class(input_shapefile)=="SpatialPolygonsDataFrame"|class(input_shapefile)=="SpatialPointsDataFrame"){
    model.shp <- input_shapefile
  }else{
    model.shp <- shapefile(input_shapefile)
  }
  #check logical structure
  if(length(model.shp)<= kernel_bandwidth){
    stop("kernel_bandwidth too large for input_shapefile features number")
  }
  #assign rownames
  rownames(model.shp@data) <- 1:length(model.shp)
  #remove columns?
  if(!is.na(remove_columns)[1]){
    if(length(grep(paste(remove_columns,collapse="|"),names(model.shp))) != 0){
      model.shp <- model.shp[,!names(model.shp) %in% remove_columns]
    }else{
      stop("remove_columns not found at input_shapefile. Verify its names.")
    }
  }
  #test + get dependent column
  model.dep <- grep(paste0("^",dependent_varName,"$"),names(model.shp))
  if(length(model.dep)==0){
    stop("dependent_varName not found")
  }else if(length(model.dep)>=2){
    stop("Found two or more column names at input_shapefile for the specified dependent_varName. Rename it.")
  }else if(length(unique(model.shp@data[,model.dep]))>=21){
    warning(paste0("dependent_varName has ",length(unique(model.shp@data[,model.dep])),
                   " classes. Procede with caution or reduce them into around 10 interpretable classes"))
  }
  #get independent columns
  model.ind <- names(model.shp)[!grepl(dependent_varName,names(model.shp))]
  model.ind <- grep(paste(model.ind,collapse="|"),names(model.shp))
  #put as factor
  model.shp <- model.shp[,c(model.dep,model.ind)]
  fac.vars <- sapply(model.shp@data,class)=="character"
  if(any(fac.vars)){
    model.shp@data[fac.vars] <- as.data.frame(lapply(model.shp@data[fac.vars],as.factor))
  }
  #complete cases
  pos.NA <- which(!complete.cases(model.shp@data))
  if(length(pos.NA)!=0){
    model.shp <- model.shp[which(complete.cases(model.shp@data)),]
    warning(paste0("input_shapefile has ",length(pos.NA)," incomplete case(s). Removing it/them..."))
  }
  #number and names of classses
  dep.len <- nlevels(model.shp@data[,1])
  dep.nam <- levels(model.shp@data[,1])

  ##### FUNCTIONS ####

  save.shp <- function(x,outN){
    #merge
    gwc.data <- lapply(gwc.extract,"[[",x)
    gwc.data <- do.call("rbind.data.frame",gwc.data)
    gwc.data$ID_row <- rownames(model.shp@data)
    gwc.data <- gwc.data[,c(ncol(gwc.data),1:(ncol(gwc.data)-1))]
    #save
    output.shp <- model.shp
    output.shp@data <- gwc.data
    output.name <- paste0(output_folder,"/GWRFC_",
                          ifelse(kernel_adaptative,"ADP_","FIX_"),
                          kernel_bandwidth,"_",
                          conv.name(kernel_function),
                          paste0("_",outN,".shp"))
    shapefile(output.shp,output.name,overwrite=T)
    print(paste0(output.name," * stored sucessfully!"))
  }
  conv.name <- function(x){
    if(x=='gaussian'){
      x <- "GA"
    }else if(x=='exponential'){
      x <- "EX"
    }else if(x=='bisquare'){
      x <- "BI"
    }else if(x=='bisquare'){
      x <- "TR"
    }
    return(x)
  }
  zero.variance <- function(x){
    nzv <- caret::nearZeroVar(x[-1], saveMetrics= TRUE)
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
    cell.vars <- sort(names(model.shp@data[,2:ncol(model.shp@data)]))
    #lvi
    cell.lvi <- as.data.frame(matrix(nrow=1,ncol=length(cell.vars)))
    names(cell.lvi) <- cell.vars
    #yhat
    cell.yhat <- as.data.frame(matrix(nrow=1,ncol=length(cell.vars)))
    names(cell.yhat) <- cell.vars
    #pdp
    cell.pdp <- as.data.frame(matrix(nrow=1,ncol=length(cell.vars)))
    names(cell.pdp) <- cell.vars
    #acc
    cell.acc <- as.data.frame(matrix(nrow=1,ncol=length(levels(model.shp@data[,1]))+4))
    names(cell.acc) <- c("DEP","PRED",paste0("P_",levels(model.shp@data[,1])),"FAIL","KAPPA")
    #rf model
    cell.rf <- NA
    #out
    cell.out <- list(cell.lvi,cell.yhat,cell.pdp,cell.acc,cell.rf)
    return(cell.out)
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
  #start cluster
  cl <- makeCluster(number_cores)
  registerDoParallel(cl)
  gwc.extract <- foreach(i=1:length(model.shp),
                         .packages=c("ranger","scales","caret","GWmodel","pdp","raster","sp"),
                         .errorhandling="pass") %dopar% {
    #timer
    init.time <- proc.time()
    #order data by distance
    cell.data <- model.shp@data
    cell.data$dist <- GWmodel::gw.dist(dp.locat=coordinates(model.shp),rp.locat=coordinates(model.shp[i,]))[,1]
    if(kernel_adaptative){
      cell.data <- cell.data[order(cell.data$dist)[1:kernel_bandwidth],]
    }else{
      cell.data <- cell.data[cell.data$dist < kernel_bandwidth,]
    }
    #get 'i' observation & rows index
    cell.i <- cell.data[1,]
    cell.pos <- as.numeric(rownames(cell.data))
    #CORRUPTED CASE 1: only one observation
    if(nrow(cell.data)==1){
      cell.out <- corrupted.cases()
    }else{
      #drop unused levels
      cell.data[,1] <- droplevels(cell.data[,1])
      #balance + set levels in classification
      if(upsampling){
        cell.data <- caret::upSample(x=cell.data[,2:ncol(cell.data)],
                                     y=cell.data[,1],
                                     yname=names(cell.data)[1])
      }else{
        cell.data <- caret::downSample(x=cell.data[,2:ncol(cell.data)],
                                       y=cell.data[,1],
                                       yname=names(cell.data)[1])
      }
      cell.data <- cell.data[,c(ncol(cell.data),1:(ncol(cell.data)-1))]
      #distance weights
      cell.weights <- GWmodel::gw.weight(cell.data$dist,
                                         bw=kernel_bandwidth,
                                         kernel=kernel_function,
                                         adaptive=kernel_adaptative)
      #assign dependent + remove zero variance
      cell.data <- cell.data[,1:(ncol(cell.data)-1)]
      cell.vars <- sort(names(cell.data[,2:ncol(cell.data)]))
      cell.data <- cbind(cell.data[,1],zero.variance(cell.data[,-1]))
      names(cell.data)[1] <- dependent_varName
      #apply ranger: FIRST TIME
      cell.formula <- formula(paste0(dependent_varName," ~ ."))
      cell.rf <- ranger::ranger(formula=cell.formula,
                                data=cell.data,
                                replace=T,
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
          cell.rf <- ranger::ranger(formula=cell.formula,
                                    data=cell.data,
                                    replace=T,
                                    min.node.size=1,
                                    scale.permutation.importance=T,
                                    case.weights=cell.weights,
                                    probability=T,
                                    importance="permutation")
        }
        #extract LVI
        cell.lvi <- data.frame(t(matrix(cell.rf$variable.importance)))
        names(cell.lvi) <- names(cell.rf$variable.importance)
        cell.lvi[,cell.lvi < 0] <- 0
        cell.rem <- cell.vars[!cell.vars %in% names(cell.lvi)]
        if(length(cell.rem!=0)){
          cell.val <- as.data.frame(t(matrix(rep(NA,length(cell.rem)))))
          names(cell.val) <- cell.rem
          cell.lvi <- cbind(cell.lvi,cell.val)
        }
        cell.lvi <- cell.lvi[sort(names(cell.lvi))]
        #apply PDP
        if(enable_pdp){
          cell.pdp <- lapply(names(cell.rf$variable.importance),function(x){
            x <- pdp::partial(cell.rf,
                              pred.var = x,
                              grid.resolution = 5,
                              train = cell.data,
                              type = "classification",
                              parallel = F)
            x <- x[which.max(x[,2]),]
            return(x)
          })
          #extract yhat
          cell.yhat <- as.data.frame(lapply(cell.pdp,"[[",2),stringsAsFactors=F)
          names(cell.yhat) <- names(cell.rf$variable.importance)
          if(length(cell.rem!=0)){
            cell.val <- as.data.frame(t(matrix(rep(NA,length(cell.rem)))))
            names(cell.val) <- cell.rem
            cell.yhat <- cbind(cell.yhat,cell.val)
          }
          cell.yhat <- cell.yhat[sort(names(cell.yhat))]
          #extract max values
          cell.pdp <- as.data.frame(lapply(cell.pdp,"[[",1),stringsAsFactors=F)
          names(cell.pdp) <- names(cell.rf$variable.importance)
          if(length(cell.rem!=0)){
            cell.val <- as.data.frame(t(matrix(rep(NA,length(cell.rem)))))
            names(cell.val) <- cell.rem
            cell.pdp <- cbind(cell.pdp,cell.val)
          }
          cell.pdp <- cell.pdp[sort(names(cell.pdp))]
        }else{
          cell.yhat <- as.data.frame(matrix(nrow=1,ncol=length(cell.vars)))
          names(cell.yhat) <- cell.vars
          cell.pdp <- as.data.frame(matrix(nrow=1,ncol=length(cell.vars)))
          names(cell.pdp) <- cell.vars
        }
        #extract accuracies
        cell.acc <- data.frame(DEP=as.character(cell.i[,1]))
        cell.acc$PRED <- get.probsClass()[[1]]
        cell.acc.prob <- get.probsClass()[[2]]
        if(length(cell.acc.prob)!=dep.len){
          miss.class <- as.data.frame(matrix(rep(0,dep.len),nrow=1,ncol=dep.len))
          names(miss.class) <- dep.nam
          miss.class[which(dep.nam %in% names(cell.acc.prob))] <- cell.acc.prob
          cell.acc.prob <- miss.class
        }
        names(cell.acc.prob) <- paste0("P_", names(cell.acc.prob))
        cell.acc <- cbind(cell.acc,as.data.frame(cell.acc.prob))
        cell.acc$FAIL <- ifelse(cell.acc$DEP!=cell.acc$PRED,"yes","no")
        cell.acc$KAPPA <- get.kappa()
      }
      #close
      cell.out <- list(cell.lvi,cell.yhat,cell.pdp,cell.acc,cell.rf)
    }
    #status text
    init.time <- round(proc.time()-init.time,2)
    cat(paste0(as.character(i)," of ",length(model.shp)," - time elapsed: ",init.time[3],"\n"),
        file=paste0(output_folder,"/data_progress.txt"), append=TRUE)
    #end
    return(cell.out)
  }
  stopCluster(cl)

  #### SAVE SHAPEFILES ####

  print("Start saving...")

  save.shp(1,"LVI")
  if(enable_pdp){
    save.shp(2,"YHAT")
    save.shp(3,"PDP")
  }
  save.shp(4,"ACC")

  #### SAVE RF MODELS ####

  if(save_models){
    rf.models <- lapply(gwc.extract,"[[",5)
    rf.models <- rf.models[!is.na(rf.models)]
    output.name <- paste0(output_folder,"/GWRFC_",
                          ifelse(kernel_adaptative,"ADP_","FIX_"),
                          kernel_bandwidth,"_",
                          conv.name(kernel_function),
                          paste0("_MODELS.rds"))
    saveRDS(rf.models,output.name)
    print(paste0(output.name," * stored sucessfully!"))
  }

  #error warning
  if(length(which(is.na(model.shp@data$PRED)))!=0){
    warning(paste0(length(which(is.na(model.shp@data$PRED)))," observations were not possible to evaluate during Random Forest execution."))
  }
  #end
  print("****GWRFC end sucessfully*****")
}
