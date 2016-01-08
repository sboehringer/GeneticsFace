#
#	faceAnalysis.R
#Fri Jan  8 17:24:31 CET 2016

if (0) {
	fts = extractFeaturesArray(r$coords, c('coordinate', 'distance', 'area', 'angle'),
		structure = struct, symmetries = symms);
	dataFeature = dataComponents(fts, 'feature');
	dataFull = dataComponents(fts, c('feature', 'asymm'));
}

if (0) {
	group = as.integer(FetchRegexpr('mut_neg', r$images, ret.all = T) == '');
	rClass = classifyFaceFeatures(group, dataFeature$data, Nrepeat = 1);
	rClassT = cbind(group, unlist(rClass$prediction));
	print(rClassT);
	print(sum(matrix.same(rClassT))/nrow(rClassT));
}

if (0) {
	cfsL = lapply(c('coordinate', 'distance', 'area', 'angle'), function(feature) {
		print(feature);
		extractFeatureCoefficients(rClass$model, feature, type = 'feature', dataFeature$desc)
	});
}
if (0) {
	save(gr, rClass, dataFeature, file = "~/Dropbox/Consultation/Bruna/classifier.RData");
}
if (0) {
# 
# 	# Example with number of colors
# 	n.col=pars$n.col
# 	color=hsv(h = seq(0.7,0,length.out=n.col), s = 1, v = 1, alpha=0.5)
# 	plot(c(0,1,1,0,0),c(0,0,1,1,0),type='l',ylim=c(0,n.col), xaxt='n', yaxt='n', ann=FALSE )
# 	for (i in 1:n.col) polygon(c(0,1,1,0,0),c(0,0,1,1,0)+i-1, col =color[i], lty = 1, lwd = 2, border = 1)
# 	axis(2,at=c(.5,n.col-1),labels=c('Least Important','Most Important'))
	averDir = r$outputAverages;
	outputDir = Sprintf('%{output}s/06_importance', output = r$output);
	Dir.create(outputDir);
	sapply(levels(groups), function(group) {
		grs = symmetrizedAverageGraph(r$coords[,,groups == group], flip = TRUE, extend = globalExtend);
		average = channel(readImage(files = Sprintf('%{averDir}s/%{group}s.tif')), 'rgb');  
		importancePlot(grs, rClass$model, dataFeature$desc, average = average,
			output = Sprintf('%{outputDir}s/importance-%{group}s'),
			pars = list(TRIANGULATION = T));
	});
}
if (1) {
	#model = cbind(rClass$model, rClass$model);
	m1 = m = glmnetModel(rClass$model);
	#m1[T] = 0;
	#m1[coefficientIndeces('angle', 'feature', dataFeature$desc), ] = 1;
	#m1[coefficientIndeces('distance', 'feature', dataFeature$desc), ] = 1;
	importancePlots(
		r$coords, r$group, m1, dataFeature$desc,
		output = r$output, averageInput = r$outputAverages
	);
}
if (0) {
	save(r, model, dataFeature, file = "~/Dropbox/Consultation/Bruna/classifier-data.RData");
}

if (0) {
	grAS = symmetrizedAverageGraph(r$coords, flip = TRUE, extend = globalExtend);
	nodes = pointsGrob(grAS[, 1], grAS[, 2],
		size = unit(5, 'points'), pch = 21, gp = gpar(fill = rgb(1, 0, 0)), default.units = 'points');
	angles = graphAngleSegments(grAS, struct$angle, radius = -.3);
	orient = triangleOrientation(grAS, struct$angle);
	#r = circleSegment(c(0, 0), 1, c(0, pi));
	tri = graphTriangulationSegments(grAS, struct$angle);
	grid.newpage();
	#Grid.draw(list(angles));
	Grid.draw(list(nodes, tri, angles));
}
