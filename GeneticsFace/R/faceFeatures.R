#
#	faceFeatures.R
#Thu Sep 24 16:03:45 2015

#
#	<p> graph based helper functions
#

meanGraph = function(coords)apply(coords, 1:2, mean);
graphDimensions = function(graph) {
	mn = apply(graph, 2, min);
	mx = apply(graph, 2, max);
	extend = mx - mn;
	rbind(mn, mx, extend);
}

#
#	<p> extract features
#

featureStructureCoordinate = function(graph, direction) {
	nodes = as.matrix(1:nrow(graph));
	nodes
}

featureStructureDistance = function(graph, direction) {
	nodes = as.matrix(1:nrow(graph));
	# <p> distances
	distances = do.call(rbind, lapply(as.list(set_combn(nodes, 2L)), unlist));
	distances
}

featureStructureArea = function(graph, direction) {
	triangles = delaunaySymm(graph, direction = direction)$triangles;
	triangles
}

featureStructureAngle = function(graph, direction) {
	triangles = delaunaySymm(graph, direction = direction)$triangles;
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

symmetrizeFeature = function(feature, structure, symm) {
	# <p> symmetrical features
	featureS = apply(symm$pairs, 1, function(pair) {
		mn = mean(feature[pair]);
		asymm = feature[pair[1]] - (if (pair[1] == pair[2]) symm$selfSymmetryRef else mn);
		c(mn, asymm)
	});
	# <p> non-symmetrical features
	nonpaired = setdiff(1:nrow(structure), unique(as.vector(symm$pairs)));
	# <p> result
	r = list(feature = c(featureS[1, ], feature[nonpaired]), asymm = featureS[2, ]);
	r
}

extractDistance = function(graph, structure, symmetries = NULL) {
	ds = apply(structure, 1, function(pair) {
		vnormv(graph[pair[1], ] - graph[pair[2], ])
	});
	if (is.null(symmetries)) return(list(feature = ds));
	r = symmetrizeFeature(ds, structure, symmetries$distance);
	r
}

#
#	<p> Triangle notation
#
#		  * C
#		/   \
#	   /	 \
#	A *-------* B
#	A := tri[1], B := tri[2], C := tri[3]
#	a:= BC = d[1, ], b:= AC = d[2, ], c:= AB = d[3, ]
#	al = alpha := ∠BAC, be = beta := ∠ABC, ga = gamma := alpha := ∠ACB
#
#	External Angles are relative to the first coordinate axes after choosing point Z as origin
#		   * C
#		  /  \
#		 /    \
#		/	 --* B
#	   /  --/	
#	A *--/------------------> X
#	ext(A) = alpha external := ∠XAB, ext(B) := ∠XBC, ext(C) := alpha := ∠XCA
#

triangleDistances = function(graph, triangles) {
	d = apply(triangles, 1, function(tri)c(
		vnormv(graph[tri[3], ] - graph[tri[2], ]),	# a
		vnormv(graph[tri[3], ] - graph[tri[1], ]),	# b
		vnormv(graph[tri[2], ] - graph[tri[1], ])	# c
	));
	dimnames(d) = list(c('a', 'b', 'c'), NULL);
	d
}

triangleOrientation = function(graph, triangles) {
	apply(triangles, 1, function(tri)sign(
		vectorCross(graph[tri[['b']], ] - graph[tri[['a']], ], graph[tri[['c']], ] - graph[tri[['a']], ])
	[3]))
}

# structure: triangulation as computed by delaunaySymm
extractArea = function(graph, structure, symmetries = NULL) {
	d = triangleDistances(graph, structure);
	s = apply(d, 2, sum)/2;
	area = sqrt(s * (s - d['a', ]) * (s - d['b', ]) * (s - d['c', ]));
	if (is.null(symmetries)) return(list(feature = area));
	r = symmetrizeFeature(area, structure, symmetries$triangle);
	r
}

extractAngleExternalX = function(graph, structure, symmetries = NULL) {
	d = triangleDistances(graph, structure);
	sta = structure[, 'a'];
	stb = structure[, 'b'];
	stc = structure[, 'c'];
	g1 = graph[, 1];
	g2 = graph[, 2];
	extAr = acos((g1[stb] - g1[sta]) / d['c', ]);
	extA = ifelse((g2[stb] - g2[sta]) < 0, 2*pi - extAr, extAr);
	extBr = acos((g1[stc] - g1[stb]) / d['a', ]);
	extB = ifelse((g2[stc] - g2[stb]) < 0, 2*pi - extBr, extBr);
	extCr = acos((g1[sta] - g1[stc]) / d['b', ]);
	extC = ifelse((g2[sta] - g2[stc]) < 0, 2*pi - extCr, extCr);
	angle = vector.intercalate(extA, extB, extC);

	if (is.null(symmetries)) return(list(feature = angle));
	r = symmetrizeFeature(angle, structure, graphSymmetryExpand(symmetries$triangle, N = 3));
	r
}

extractAngleMinDist = function(graph, structure, symmetries = NULL) {
	d = triangleDistances(graph, structure);
	distA = apply(cbind(d['b', ], d['c', ]), 1, min);
	distB = apply(cbind(d['a', ], d['c', ]), 1, min);
	distC = apply(cbind(d['a', ], d['b', ]), 1, min);
	distMin = vector.intercalate(distA, distB, distC);

	if (is.null(symmetries)) return(list(feature = distMin));
	r = symmetrizeFeature(distMin, structure, graphSymmetryExpand(symmetries$triangle, N = 3));
	r
}

extractAngle = function(graph, structure, symmetries = NULL) {
	d = triangleDistances(graph, structure);
	dq = d^2;
	al = acos((dq['b', ] + dq['c', ] - dq['a', ]) / (2 * d['b', ] * d['c', ]));
	be = acos((dq['a', ] + dq['c', ] - dq['b', ]) / (2 * d['a', ] * d['c', ]));
	ga = acos((dq['a', ] + dq['b', ] - dq['c', ]) / (2 * d['a', ] * d['b', ]));
	angle = vector.intercalate(al, be, ga);

	if (is.null(symmetries)) return(list(feature = angle));
	r = symmetrizeFeature(angle, structure, graphSymmetryExpand(symmetries$triangle, N = 3));
	r
}

extractFeatures = function(graph, features, structure, symmetries) {
	fs = nlapply(features, function(feature) {
		extractor = get(Sprintf('extract%{feature}u'));
		extractor(graph, structure[[feature]], symmetries);
	});
}

featuresDesc = function(types, structure, symmetries) {
	list(features = types, structure = structure, symmetries = symmetries);
}

extractFeaturesArray = function(coords, features = 'distance', structure, symmetries = NULL) {
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

	r = list(feature = ftsM,
		asymm = ftsAM,
		desc = c(featuresDesc(features, structure, symmetries), list(
			indeces = list(
				feature = vectorNamed(sapply(fts, ncol), features),
				asymm = vectorNamed(sapply(ftsA, ncol), features)
			)
		))
	);
	r
}

#
#	<p> create data including meta-data
#

# create full matrix and descriptor based on components
dataComponents = function(fts, components) {
	data = do.call(cbind, fts[components]);
	desc = fts$desc;
	desc$indeces = desc$indeces[components];
	r = list(data = data, desc = desc);
	r
}

structureForType = function(structure, symmetries, type = 'feature') {
	struct = if (!is.null(symmetries)) {
		symmF = symmetries$pairs[, 1];
		# asymmetry features do not include unpaired features
		nonpaired = if (type == 'asymm') c() else
			setdiff(1:nrow(structure), unique(as.vector(symmetries$pairs)));
		structure[c(symmF, nonpaired), , drop = FALSE];
	} else structure;
	struct
}

symmetryMap = list(coordinate = 'node', distance = 'distance', area = 'triangle', angle = 'triangle');
extractStructureFromDesc = function(feature = 'distance', type = 'feature', desc) {
	structureForType(desc$structure[[feature]], desc$symmetries[[symmetryMap[[feature]]]], type);
}

# type can be feature, asymm, for combined data set assume asymm after feature
extractFeatureCoefficients = function(model, feature = 'distance', type = 'feature', desc) {
	cs = cumsumR(desc$indeces);
	csFeature = cumsumI(desc$indeces[[type]], offset = 0);
	csI = which(names(csFeature) == feature);
	# <N> assume intercept included into model
	cfs = model[1 + cs[[type]] + (csFeature[csI]:(csFeature[csI+1] - 1)) ,];
	struct = extractStructureFromDesc(feature, type, desc);
	r = list(coefficients = cfs, structure = struct);
	r
}


