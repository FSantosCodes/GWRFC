#'@title Geographically weighted Random Forest Classification (GWRFC)
#'@description GWRFC is function that replaces the linear regression model of the Geographically Weighted Regression (GWR; Fotheringham, Charlton, and Brunsdon 1998) with the random forest algorithm (RF; Breiman 2001). For this, it applies case weights according to the weightening scheme of GWR in the bagging step of RF. As a result, GWRFC  produces spatial representations of variables importance, classification probabilities and accuracy of RF models at local level.
#'@param input_shapefile string or Spatial-class. Input shapefile with dependent and independent variables.  It can be the filename of the shapefile or an object of class {SpatialPolygonsDataFrame} or {SpatialPointsDataFrame}.
#'@param remove_columns string. Remove specific variables from \strong{input_shapefile}. Variables are identified by column name. NA ignores column remove.
#'@param dependent_varName string. Dependent variable name. Must exists at \strong{input_shapefile} and should be categorical (with not more than 20 classes).
#'@param kernel_type string. Kernel type to apply in GWRFC. It can be: 'gaussian', 'exponential', 'bisquare' or 'tricube'.
#'@param kernel_adaptative logical. Is the kernel adaptative? otherwise it is considered as fixed (larger processing time).
#'@param kernel_bandwidth numeric. Defines kernel bandwidth. If \strong{kernel_adaptative} is TRUE, then you should define the number of local observations in the kernel, otherwise you should define a distance to specify kernel bandwidth.
#'@param upsampling logical. If TRUE, upsampling is applied before random forest training, otherwise it is downsampled. Consider that upsampling is a bit more computing demanding but accuracy is improved.
#'@param number_cores numeric. Number of cores for parallel processing. Cores are register and operated via doParallel, foreach and parallel packages. Be careful with increasing numbers of cores, as RAM memory may be not enough.
#'@param output_folder string. Output folder where GWRFC outputs will be stored.
#'@return As a result, a shapefile is created whose attribute table contains: \enumerate{
#'                            \item Local variables importance: calculated via permutation for each variable and derived from each RF local models.
#'                            \item BEST: most important variable in the local RF model.
#'                            \item DEP: original value from dependent variable.
#'                            \item PRED: predicted class result.
#'                            \item P_: classification probabilities for each dependent variable classes.
#'                            \item FAIL: Is prediction result correct (as is compared with DEP)?
#'                            \item KAPPA: accuracy according kappa index obtained from RF local models.
#'                            \item ID_row: an identifier to link rows with the original dataset if incomplete cases are found.
#'                          }
#'        Aditionally, you can check processing evolution for the parallel computing at \strong{output_folder} as: data_progress.txt
#'@examples
#'#with deforestation dataset
#'data(deforestation)
#'
#'#as the dependent variable is a continuous value, here it is splited in quantiles
#'deforestation@data$fao <- factor(cut(deforestation@data$fao,breaks=quantile(deforestation@data$fao,probs=seq(0,1,length.out=5)),labels=c("Q1","Q2","Q3","Q4"),include.lowest=T))
#'
#'#and then GWRFC is run (see coments for each parameter)
#'GWRFC(input_shapefile = deforestation, #can be also a complete filename of .shp extension.
#'       remove_columns = c("ID_grid","L_oth"), #these variables are ignored in the analysis as are not informative.
#'       dependent_varName = "fao", #the depedent variable to evaluate, should be of class factor or character.
#'       kernel_type = "exponential", #the weightening function. See above for other functions.
#'       kernel_adaptative = T, #TRUE for adaptative or FALSE for a fixed kernel distance.
#'       kernel_bandwidth = 400, #as the kerner is adaptative, 400 refers to the minimun number of observations.
#'       upsampling = T, #improves accuracy but is a bit more computing costly.
#'       number_cores = 3, #3 cores used in an AMD A6/16 GB RAM computer.
#'       output_folder = "C:/DATA/demo/deforestation") #check this folder for outputs
#'@export

GWRFC <- function(
  input_shapefile,
  remove_columns = NA,
  dependent_varName,
  kernel_type = "exponential",
  kernel_adaptative = T,
  kernel_bandwidth,
  upsampling = T,
  number_cores = 1,
  output_folder
){

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
    cell.names <- sort(names(model.shp@data[,2:ncol(model.shp@data)]))
    probs.class <- paste0("P_",levels(model.shp@data[,1]))
    cell.names <- c(cell.names,"BEST","DEP","PRED",probs.class,"FAIL","KAPPA")
    cell.out <- as.data.frame(matrix(nrow=1,ncol=length(cell.names)))
    names(cell.out) <- cell.names
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
  gwc.extract <- foreach(i=icount(length(model.shp)),.packages=c("ranger","scales","caret","GWmodel"),.errorhandling="pass") %dopar% {
    init.time <- proc.time()
    #order data by distance
    cell.data <- model.shp@data
    cell.data$dist <- GWmodel::gw.dist(dp.locat=coordinates(model.shp),rp.locat=coordinates(model.shp[i,]))[,1]
    if(kernel_adaptative){
      cell.data <- cell.data[order(cell.data$dist)[1:kernel_bandwidth],]
    }else{
      cell.data <- cell.data[cell.data$dist < kernel_bandwidth,]
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
                                         kernel=kernel_type,
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
        #get best variable
        cell.out$BEST <- head(names(cell.out[order(cell.out,decreasing=T)]),3)[1]
        #get dependent
        cell.out$DEP <- as.character(cell.i[,1])
        #get prediction
        cell.out$PRED <- get.probsClass()[[1]]
        #get probabilities
        cell.out.prob <- get.probsClass()[[2]]
        if(length(cell.out.prob)!=dep.len){
          miss.class <- as.data.frame(matrix(rep(0,dep.len),nrow=1,ncol=dep.len))
          names(miss.class) <- dep.nam
          miss.class[which(dep.nam %in% names(cell.out.prob))] <- cell.out.prob
          cell.out.prob <- miss.class
        }
        names(cell.out.prob) <- paste0("P_", names(cell.out.prob))
        cell.out <- cbind(cell.out,as.data.frame(cell.out.prob))
        #get accuracies metrics
        cell.out$FAIL <- ifelse(cell.out$DEP!=cell.out$PRED,"yes","no")
        cell.out$KAPPA <- get.kappa()
      }
    }
    #status text
    init.time <- round(proc.time()-init.time,2)
    cat(paste0(as.character(i)," of ",length(model.shp)," - time elapsed: ",init.time[3],"\n"),
        file=paste0(output_folder,"/data_progress.txt"), append=TRUE)
    #end
    return(cell.out)
  }
  stopCluster(cl)
  #consolidate results
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
