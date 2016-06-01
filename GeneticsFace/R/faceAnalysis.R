#
#	faceAnalysis.R
#Fri Jan  8 17:24:31 CET 2016

faceClassifyAndVisualize = function(state, structure, symmetries,
	features = c('coordinate', 'distance', 'area', 'angle'), component = 'feature', Nrepeat = 1) {

	# <p> classification
	fts = extractFeaturesArray(state$coords, features, structure = struct, symmetries = symms);
	data = dataComponents(fts, component);
	rClass = classifyFaceFeatures(state$group, data$data, Nrepeat = Nrepeat);
	rClassT = cbind(as.character(state$group), unlist(rClass$prediction));
	classification = list(
		accuracy = mean(matrix.same(rClassT)),
		rClass = rClass
	);

	# <p> averages
	state = prepareAveraging(state, state$coords);
	state = averageGroups(state, state$group, Niter = 20);

	# <p> importance
	importancePlots(
		state$coords, state$group, glmnetModel(rClass$model), data$desc,
		output = state$output, averageInput = state$outputAverages
	);
	# <p> feature counts
	graph = symmetrizeGraph(meanGraph(state$coords));
	countFeatures = featureCounts(graph, classification, features, component, struct, symms)

	r = c(state, list(classification = classification,
		features = features, components = component,
		countFeatures = countFeatures
	));
	r
}
