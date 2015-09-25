#
#	faceFeatures.R
#Thu Sep 24 16:03:45 2015

featureStructureCoordinate = function(graph, direction) {
	nodes = as.matrix(1:nrow(graph));
	nodes
}

featureStructureDistance = function(graph, direction) {
	# <p> distances
	distances = do.call(rbind, lapply(as.list(set_combn(nodes, 2L)), unlist));
	distances
}

featureStructureArea = function(graph, direction) {
	triangles = delaunaySymm(graph, direction = direction)$symm;
	triangles
}

featureStructureAngle = function(graph, direction) {
	triangles = delaunaySymm(graph, direction = direction)$symm;
	triangles
}

featureStructure = function(graph, direction = 1, features = c('coordinate', 'distance', 'area', 'angle')) {
	fs = nlapply(features, function(feature) {
		feature = get(Sprintf('featureStructure%{feature}u'));
		feature(graph, direction = direction);
	});
}

extractCoordinate = function(graph, structure, symmetries = NULL) {
	coords = apply(structure, 1, function(i)graph[i, ]);
	if (is.null(symmetries)) return(list(feature = ds));
	symm = symmetries$node;
	# <p> symmetrical features
	coordsS = apply(symm$pairs, 1, function(pair) {
		mnSelf = mn = apply(coords[, pair], 1, mean);
		mnSelf[symm$direction] = symm$selfSymmetryRef;
		asymm = coords[, pair[1]] - (if (pair[1] == pair[2]) mnSelf else mn);
		c(mn, asymm)
	});
	# <p> non-symmetrical features
	nonpaired = setdiff(1:nrow(structure), unique(as.vector(symm$pairs)));
	# <p> result
	Ndim = nrow(coordsS)/2;
	r = list(
		feature = c(as.vector(coordsS[1:Ndim, ]), as.vector(coords[, nonpaired])),
		asymm = as.vector(coordsS[(Ndim + 1):nrow(coordsS), ])
	);
	r
}

extractDistance = function(graph, structure, symmetries = NULL) {
	ds = apply(structure, 1, function(pair) {
		vnormv(graph[pair[1], ] - graph[pair[2], ])
	});
	if (is.null(symmetries)) return(list(feature = ds));
	symm = symmetries$distance;
	# <p> symmetrical features
	dsS = apply(symm$pairs, 1, function(pair) {
		mn = mean(ds[pair]);
		asymm = ds[pair[1]] - (if (pair[1] == pair[2]) symm$selfSymmetryRef else mn);
		c(mn, asymm)
	});
	# <p> non-symmetrical features
	nonpaired = setdiff(1:nrow(structure), unique(as.vector(symm$pairs)));
	# <p> result
	r = list(feature = c(dsS[1, ], ds[nonpaired]), asymm = dsS[2, ]);
	r
}


extractFeatures = function(graph, features, structure, symmetries) {
	fs = nlapply(features, function(feature) {
		extractor = get(Sprintf('extract%{feature}u'));
		extractor(graph, structure[[feature]], symmetries);
	});
}

extractFeaturesArray = function(coords, features = 'distance', structure = NULL, symmetries = NULL) {
	# <p> compute structure if not passed as argument
	if (is.null(structure)) structure = featureStructure(coords, features = features);
	# <p> subset structure to needed extraction
	structure = structure[features];
	# <p> extract features
	r = apply(coords, 3, function(graph) {
		extractFeatures(graph, features, structure, symmetries);
	});
	# <p> re-packages results
	fts = lapply(features, function(feature)do.call(rbind, list.kp(r, Sprintf('%{feature}s$feature'))));
	ftsM = do.call(cbind, fts);
	if (is.null(symmetries)) return(list(feature = ftsM));
	# <p> asymmetries
	ftsA = lapply(features, function(feature)do.call(rbind, list.kp(r, Sprintf('%{feature}s$asymm'))));
	ftsAM = do.call(cbind, ftsA);

	r = list(feature = ftsM, asymm = ftsAM);
	r
}
