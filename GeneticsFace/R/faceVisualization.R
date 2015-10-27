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

coloredPlots = function(feature, meanGraph, modelDesc, colorCoefficients, pars){
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

	tmpImg = tempfile(fileext = '.png');
	png(tmpImg);
	par(mar=rep(0, 4), xpd = NA) 
	image(x,y,z = t(colorFeature), col=color,bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n') #, asp=800/800
	if (pars$TRIANGULATION) {
		points(meanGraph[,1], meanGraph[,2],col="red",xlab="",ylab="", main="")
		#textxy(meanGraph[,1], meanGraph[,2],labs=1:nrow(meanGraph), cex = 1, col = "red")
		for(i in 1:nrow(modelDesc$structure$area)) {
			ena=modelDesc$structure$area[i,]
			ena[4]=modelDesc$structure$area[i,1]
			lines(meanGraph[ena,1], meanGraph[ena,2],lwd=2)
		}
	}
	dev.off();
	readImage(tmpImg);
}

collapseColorAverg=function(input, meanGraph, average, pars) with(pars, {
	gdim = graphDimensions(meanGraph);
	odim = dim(average)[1:2];
	MIXcol = 1 - MIXave;
	importance = resize(readImage(input), gdim['extend', 1], gdim['extend', 2]);
	# Image lives in IVth quadrant
	importanceT = translate(importance, -c(gdim['mn', 1], odim[2] - gdim['mx', 2]), output.dim = odim);
	importancePlot = (MIXave * average^POWERave + MIXcol * importanceT[,,1:3]^POWERcol);
	importancePlot
})

parsDefault = list(
	Npoints = 100,			# number of points used to construct the (Npoints x Npoints) grid 
	halfFace = FALSE,		# should grid only be constructed for half of the face
	midPoint=16,			# if halfFace=T, what is the point in the center of the X axis?
	n.col = 100,			# number of different hues (see example below) for function hsv()
	s=1, v=1,				# numeric values in the range [0, 1] for saturation and hue value 
							# to be combined to form a vector of colors for function hsv()
	alpha=.5,				# numeric vector of values in the range [0, 1] for alpha transparency 
							# channel (0 means transparent and 1 means opaque). Needed for choosing 
							# transparency of the colored graph in order to see the average black 
							# and white image in the background
	STAND = TRUE,			# should the color coefficients be standardized. need for plots to be comparable
	NORM = TRUE,			# should the color coefficients be normalized need for plots to be comparable
	SCALAR = 1,				# numeric value in the range of [0, Inf] used for normalization of color coefficient 
							# exp(color*SCALAR)/(1/exp(color*SCALAR))
	MIRRORED = TRUE,		# should the colors from half of the face be mirrored to the other half? 
	TRIANGULATION = T,	# should the plot of the average individual appear on top of the colored image? 
					        # Useful to see if plots are aligned properly
	MIXave=.5,				# numeric values in the range [0, 1] for mixing parameter for the colored 
							# and background photo. 
	POWERave = 2,			# numeric value in the range [0, Inf] for intensity of colors
							#	for the background photo.
	POWERcol = 2			# numeric value in the range [0, Inf] for intensity of colors
							#	for the importance image photo.
);

importancePlot = function(meanGraph, model, modelDesc, pars = list(), output, average = NULL) {
	pars = merge.lists(parsDefault, pars);
	colorCoefficients = gridColor(meanGraph, model, modelDesc, pars)
	r = lapply(c(modelDesc$features, 'all'), function(feature) {
		pathImportance = Sprintf('%{output}s-%{feature}s.png');
		writeImage(
			coloredPlots(feature, meanGraph, modelDesc, colorCoefficients, pars)
		, pathImportance);
		if (!is.null(average)) writeImage(
			collapseColorAverg(pathImportance, meanGraph, average, pars)
		, Sprintf('%{output}s-%{feature}s-background.png'))
	});
	r
}

importancePlots = function(coords, groups, models, modelDesc, pars = list(), output,
	averageInput = NULL, globalExtend = 512, outputImportance = '06_importance', plotTriangulation = TRUE) {
	outputDir = Sprintf('%{output}s/%{outputImportance}s');
	Dir.create(outputDir);
	ilapply(levels(groups), function(group, i) {
		grs = symmetrizedAverageGraph(coords[,, groups == group], flip = TRUE, extend = globalExtend);
		average = channel(readImage(files = Sprintf('%{averageInput}s/%{group}s.tif')), 'rgb');  
		importancePlot(grs, models[,i , drop = F], modelDesc, average = average,
			output = Sprintf('%{outputDir}s/importance-%{group}s'),
			pars = list(TRIANGULATION = plotTriangulation));
	});
}

