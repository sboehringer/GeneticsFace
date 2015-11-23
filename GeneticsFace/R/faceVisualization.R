#
#	faceVisualization.R
#	(c) 2015 Brunila Balliu

dstPointMatxPoints = function(p1, feature){
	# compute the distance of a point p1 from each row of a matrix of points
	# features is a matrix containing coordinates pairs
	sqrt(colSums((p1 - t(feature))^2))
}

centerLine = function(structure, meanGraph) {
  # compute the center of a distance
  # structure is the integer of the two points comprising the distance
  # meanGraph contains the coordinate pairs (per row) of average individual
  colMeans(meanGraph[structure, ])
}

centroidTriangle = function(structure, meanGraph) {
  # compute the centroid of a triangle
  # structure is the integer of the three points comprising the triangle
  # meanGraph contains the coordinate pairs (per row) of average individual
  colMeans(meanGraph[structure, ])
}

colorPoint = function(distances, regressCoefficients, distanceFudge = 1) {
  # compute color coefficients for one point in a grid
  # regressCoefficients: vector of regression coefficients from glmnet for a feature
  # distances: vector of distances as computed from dstPointMatxPoints
  sum(abs(regressCoefficients)/(distanceFudge + distances))            
}

gridCoords = function(meanGraph, Npoints, flip = FALSE, flop = FALSE, byIndex = FALSE) {
	x.coo = if (byIndex) 1:Npoints else seq(min(meanGraph[,1]), max(meanGraph[,1]), length.out = Npoints);
	if (flip) x.coo = rev(x.coo);
	y.coo = if (byIndex) 1:Npoints else seq(min(meanGraph[,2]), max(meanGraph[,2]), length.out = Npoints);
	if (flop) y.coo = rev(y.coo);
	grid = as.matrix(expand.grid(x.coo,y.coo)[1:(length(x.coo) * length(y.coo)), ]);
	grid
}

gridCoordsSymmetry = function(Npoints) {
	coordsO = gridCoords(Npoints = Npoints, byIndex = TRUE);
	coordsO = data.frame(coordsO, iO = 1:dim(coordsO)[1]);
	coordsF = gridCoords(Npoints = Npoints, byIndex = TRUE, flip = TRUE);
	coordsF = data.frame(coordsF, iF = 1:dim(coordsF)[1]);
	coordsM = merge(coordsO, coordsF);
	iSymm = coordsM[order(coordsM$iF), ]$iO;
	iSymm
}

# identity
featureCenterCoordinate = function(meanGraph, structure)meanGraph[rep(structure[, 1], each = 2), ];
# midpoint line segment
featureCenterDistance = function(meanGraph, structure)t(apply(structure, 1, centerLine, meanGraph));
# area centroid
featureCenterArea = function(meanGraph, structure)t(apply(structure, 1, centroidTriangle, meanGraph))
# angle vertex
featureCenterAngle = function(meanGraph, structure)meanGraph[c(t(structure)), ]


# Compute color coefficients for each point in a grid based on each feature
# Need to write description of arguments meanGraph=gr$graphs ; model= rClass$model; modelDesc=dataFeature$desc  
gridColor=function(meanGraph, model, modelDesc, pars){
	#### Construct a grid of Npoints x Npoints points
	grid = gridCoords(meanGraph, pars$Npoints);
	iSymm = gridCoordsSymmetry(pars$Npoints);
  
	#### Compute color coefficients
	colorCoefficients = nlapply(modelDesc$features, function(feature) {
		cfs = extractFeatureCoefficients(model, feature , type = 'feature', modelDesc);

		featureCenterFct = get(Sprintf('featureCenter%{feature}u'));
		featureCenter = featureCenterFct(meanGraph, cfs$structure);

		distFeature = t(apply(grid, 1, dstPointMatxPoints, feature = featureCenter));
		Logs("Feature %{feature}s; #:%{countFeatures}d; #coef: %{countCoeff}d",
			countFeatures = dim(featureCenter)[1], countCoeff = length(cfs$coefficients), logLevel = 5);

		cls = apply(distFeature, 1, colorPoint, regressCoefficients = cfs$coefficients);
		if (pars$Symmetrize) cls = cls + cls[iSymm];
		cls
	});
	colorCoefficients[['all']] = Reduce('+', colorCoefficients)
	return(colorCoefficients)
}  

coloredPlots = function(feature, meanGraph, modelDesc, colorCoefficients, pars){
	# Construct colored plots
	colorFeature = colorCoefficients[[feature]];
	color = hsv(h = seq(0.7,0,length.out=pars$n.col), s = pars$s, v = pars$v, alpha=pars$alpha);

	if (pars$STAND && (var(colorFeature)!=0)) colorFeature = scale(x = colorFeature)
	if (pars$NORM) {colorFeature=(exp(colorFeature*pars$SCALAR))/(1+exp(colorFeature*pars$SCALAR))}

	colorFeature = matrix(colorFeature, ncol=pars$Npoints, byrow = T)
	x = seq(min(meanGraph[,1]),max(meanGraph[,1]),length.out=pars$Npoints)
	y = seq(min(meanGraph[,2]),max(meanGraph[,2]),length.out=pars$Npoints)

	tmpImg = tempfile(fileext = '.png');
	png(tmpImg);
	par(mar=rep(0, 4), xpd = NA);
	image(x, y, z = t(colorFeature),
		col = color, bty ="n", axes=F, frame.plot=F, xaxt='n', ann=FALSE, yaxt='n') #, asp=800/800
	if (pars$TRIANGULATION) {
		points(meanGraph[, 1], meanGraph[, 2], col = "red", xlab="", ylab="", main="")
		textxy(meanGraph[, 1], meanGraph[, 2], labs = 1:nrow(meanGraph), cex = 1, col = "red")
		for (i in 1:nrow(modelDesc$structure$area)) {
			ena = c(modelDesc$structure$area[i, ], modelDesc$structure$area[i, 1])
			lines(meanGraph[ena, 1], meanGraph[ena, 2], lwd = 2)
		}
	}
	dev.off();
	readImage(tmpImg);
}

collapseColorAverg=function(input, meanGraph, average, pars) with(pars, {
	gdim = graphDimensions(meanGraph);
	odim = dim(average)[1:2];
	MIXcol = 1 - MIXave;
	imgInput = readImage(input);
	importance = resize(imgInput, gdim['extend', 1], gdim['extend', 2]);
	# Image lives in IVth quadrant
	importanceT = translate(importance, -c(gdim['mn', 1], odim[2] - gdim['mx', 2]), output.dim = odim);
	importancePlot = (MIXave * average^POWERave + MIXcol * importanceT[,, 1:3]^POWERcol);
	importancePlot
})

parsDefault = list(
	Npoints = 100,			# number of points used to construct the (Npoints x Npoints) grid 
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
	Symmetrize = TRUE,		# copy importance to other half of face
	TRIANGULATION = T,		# should the plot of the average individual appear on top of the colored image? 
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
		, Sprintf('%{output}s-background-%{feature}s.png'))
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

