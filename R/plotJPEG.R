#'@title Plot JPEG file using base graphics
#'@description Internal function to plot a JPEG file in R
#'@param input_jpeg string. Input filename of JPG picture.

plotJPEG <- function(input_jpeg,add=FALSE){
  jpg = readJPEG(input_jpeg, native=T) # read the file
  res = dim(jpg)[2:1] # get the resolution, [x, y]
  if (!add) # initialize an empty plot area if add==FALSE
    plot(1,1,xlim=c(1,res[1]),ylim=c(1,res[2]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(jpg,1,1,res[1],res[2])
}
