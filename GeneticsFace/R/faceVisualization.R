#
#	faceVisualization.R
#	(c) 2015 Brunila Balliu

dstPointMatxPoints=function(p1,feature){
  # compute the distance of a point p1 from each row of a matrix of points
  # features is a matrix containing coordinates pairs   
  sqrt(rowSums((p1-feature)^2))
}

centerLine=function(structure,meanGraph) {
  # compute the center of a distance
  # structure is the integer of the two points comprising the distance
  # meanGraph contains the coordinate pairs (per row) of average individual
  colSums(meanGraph[structure,])/2
}

centroidTriangle=function(structure,meanGraph) {
  # compute the centroid of a triangle
  # structure is the integer of the three points comprising the triangle
  # meanGraph contains the coordinate pairs (per row) of average individual
  colSums(meanGraph[structure,])/3
}

colorPoint=function(regressCoefficients,distances) {
  # compute color coefficients for one point in a grid
  # regressCoefficients: vector of regression coefficients from glmnet for a feature
  # distances: vector of distances as computed from dstPointMatxPoints
  sum(abs(regressCoefficients)/(1+distances))            
}

gridColor=function(meanGraph, model, modelDesc, pars){
  # Compute color coefficients for each point in a grid based on each feature
  # Need to write description of arguments meanGraph=gr$graphs ; model= rClass$model; modelDesc=dataFeature$desc  
  
  colorCoefficients=list()
  #### Construct a grid of Npoints x Npoints points
  x.coo=seq(meanGraph[16,1], min(meanGraph[,1]),length.out=pars$Npoints/2)
  y.coo=seq(min(meanGraph[,2]),max(meanGraph[,2]),length.out=pars$Npoints)
  grid=as.matrix(expand.grid(x.coo,y.coo)[(length(x.coo)*length(y.coo)):1,])
  
  #### Extract coefficients 
  cfs=lapply(modelDesc$features, function(feature) extractFeatureCoefficients(model, feature , type = 'feature', modelDesc));
  names(cfs)=modelDesc$features
  
  #### Compute color coefficients 
  coordinate.center=meanGraph[rep(cfs$coordinate$structure,each = 2),] 
  dst_p_coordinate.center=t(apply(grid,1,dstPointMatxPoints,feature=coordinate.center))
  colorCoefficients[["coordinate"]]=apply(dst_p_coordinate.center,1,colorPoint,regressCoefficients=cfs$coordinate$coefficients)    
  
  distance.center=t(apply(cfs$distance$structure,1,centerLine,meanGraph))    
  dst_p_distance.center=t(apply(grid,1,dstPointMatxPoints,feature=distance.center))
  colorCoefficients[["distance"]]=apply(dst_p_distance.center,1,colorPoint,regressCoefficients=cfs$distance$coefficients)    
  
  area.centroid=t(apply(cfs$area$structure,1,centroidTriangle,meanGraph))
  dst_p_area.centroid=t(apply(grid,1,dstPointMatxPoints,feature=area.centroid))
  colorCoefficients[["area"]]=apply(dst_p_area.centroid,1,colorPoint,regressCoefficients=cfs$area$coefficients)    
  
  angle.vertex= meanGraph[c(t(cfs$angle$structure)),]
  dst_p_angle.vertex=t(apply(grid,1,dstPointMatxPoints,feature=angle.vertex))
  colorCoefficients[["angle"]]=apply(dst_p_angle.vertex,1,colorPoint,regressCoefficients=cfs$angle$coefficients)      
  
  colorCoefficients[['all']]=Reduce('+', colorCoefficients)
  
  colorCoefficients[["x.coo"]]=x.coo  
  colorCoefficients[["y.coo"]]=y.coo  
  
  return(colorCoefficients)
}  


coloredPlots=function(feature, meanGraph, modelDesc, colorCoefficients, pars){
  # Construct colored plots
  colorFeature=colorCoefficients[[feature]]  
  color=hsv(h = seq(0.7,0,length.out=pars$n.col), s = pars$s, v = pars$v, alpha=pars$alpha)
  
  if(pars$STAND&(var(colorFeature)!=0)) colorFeature=scale(x = colorFeature)
  if(pars$NORM) {colorFeature=(exp(colorFeature*pars$SCALAR))/(1+exp(colorFeature*pars$SCALAR))}
  colorFeature=matrix(colorFeature,ncol=pars$Npoints/2,byrow=T)
  if(pars$MIRRORED) {
    colorFeature=cbind(colorFeature[,1:ncol(colorFeature)],colorFeature[,ncol(colorFeature):1])
    colorFeature=colorFeature[nrow(colorFeature):1,]
  }
  x=seq(min(meanGraph[,1]),max(meanGraph[,1]),length.out=pars$Npoints)
  y=seq(min(meanGraph[,2]),max(meanGraph[,2]),length.out=pars$Npoints)
  if(pars$TRIANGULATION){
    par(mar=rep(0, 4), xpd = NA) 
    image(x,y,z = t(colorFeature), col=color,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n') #, asp=800/800
    points(meanGraph[,1], meanGraph[,2],col="red",xlab="",ylab="", main="")
    textxy(meanGraph[,1], meanGraph[,2],labs=1:nrow(meanGraph), cex = 1, col = "red")
    for(i in 1:nrow(modelDesc$structure$area)) {
      ena=modelDesc$structure$area[i,]
      ena[4]=modelDesc$structure$area[i,1]
      lines(meanGraph[ena,1], meanGraph[ena,2],lwd=2)
    }
  }else{
    par(mar=rep(0, 4), xpd = NA) 
    image(x,y,z = t(colorFeature), col=color,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n') #, asp=800/800
  }
}


collapseColorAverg=function(feature,Average,pars){
  MIXcol=1-pars$MIXave
  Color=resize(readImage(paste(feature,'png',sep = ".")), dim(Average)[1], dim(Average)[2])
  Color2=resize(pars$MIXave*(Average)^pars$POWERave+MIXcol*(Color[,,1:3])^pars$POWERcol,dim(Average)[1], dim(Average)[2])
  writeImage(Color2, paste(feature, "comb.png",sep=".") , quality=100)
}

visualizeClassfifier = function(meanGraph, model, modelDesc, pars) {
  colorCoefficients=gridColor(meanGraph, model, modelDesc, pars)
  lapply(c(modelDesc$features, 'all'), function(feature) { 
    png(paste(feature,'png',sep = ".")) 
    coloredPlots(feature, meanGraph, modelDesc, colorCoefficients,pars) 
    dev.off()
    collapseColorAverg(feature,Average,pars)
  })
}
