#
#	faceFeatures.R
#Thu Sep 24 16:03:45 2015

featureStructureCoordinate = function(graph, direction) {
	nodes = as.matrix(1:nrow(graph));
	list(node = nodes)
}

featureStructureDistance = function(graph, direction) {
	# <p> distances
	distances = do.call(rbind, lapply(as.list(set_combn(nodes, 2L)), unlist));
	list(distance = distances)
}

featureStructureArea = function(graph, direction) {
	triangles = delaunaySymm(graph, direction = direction)$symm;
	list(area = triangles)
}

featureStructureAngle = function(graph, direction) {
	triangles = delaunaySymm(graph, direction = direction)$symm;
	list(angle = triangles)
}

featureStructure = function(graph, direction = 1, features = c('coordinate', 'distance', 'area', 'angle')) {
	fs = nlapply(features, function(feature) {
		feature = get(Sprintf('featureStructure%{feature}u'));
		feature(graph, direction = direction);
	});
}



extractFeatures = function(graph, features) {
	fs = nlapply(features, function(feature) {
		extractor = get(Sprintf('extract%{feature}u'));
		extractor(graph, features[[feature]]);
	});
}

extractFeaturesArray = function(coords, features = list()) {
	apply(coords, 3, function(graph) {
		extractFeatures(graph, features);
	});
}
