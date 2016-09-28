#
#	faceSymmetry.R
#Mon Sep 21 14:57:18 CEST 2015

# Euclidian vector norm
vnormv = function(v)as.vector(sqrt(v %*% v))
# Euclidian vector norm computed row-wise on matrix/array
vnorm = function(v) {
	if ((is.matrix(v) || is.array(v)) && dim(v) > 1)
		apply(v, 1:(length(dim(v)) - 1), vnormv) else
		vnormv(as.vector(v));
}

# align vector v to vref, so that distances are minimized
#	total distance should be 0 if vref and v are symmetric
# v is a matrix of nodes (nodes per row)
vector_align = function(vref, v) {
	Ns = 1:nrow(vref);
	is = sapply(Ns, function(i)which.min(vnorm(t(t(vref) - v[i, ]))));
	# <p> rectify degenerate cases <A> not optimal assignment in this case;
	if (any(any(duplicated(is)))) is[duplicated(is)] = setdiff(Ns, unique(is));
	if (any(any(duplicated(is)))) stop('did not produced reordering of vector');
	v[is, , drop = F]
}

graphmid = function(graph, direction = 1)mean(graph[, direction])

# graph: graph as an array (dim1: nodes, dim2: dimensions)
# elements: groups of coordinates to be checked as 2-dim matrix (dim1: group, dim2: nodes)
# direction: which symmetry to consider (direction == 1: x-flip, direction == 2, y-flip, ...)
# directionPerElement: try to flip around the mean per element to detect additional symmetries
#		(direction == 1: x-flip, direction == 2, y-flip, ...)

graphSymmetry = function(graph, elements, direction = 1, eps = 1e-5, selfSymmetry = FALSE) {
	midline = graphmid(graph, direction);
	Els = 1:nrow(elements);
	epairs = do.call(rbind, lapply(as.list(set_combn(Els, 2L)), unlist));
	# self-symmetry
	if (selfSymmetry) epairs = rbind(cbind(Els, Els), epairs);
	# <p> establish canonical order
	epairs = epairs[order(epairs[, 1]), ];
	#epairs = t(sapply(as.list(set_combn(1:nrow(elements), 2L)), avu));
	# determine symmetry
	s = apply(epairs, 1, function(p) {
		v1 = graph[elements[p[1], ], , drop = F];
		v2s = v2 = graph[elements[p[2], ], , drop = F];
		#v2[, direction] = midline - (v2[, direction] - midline);
		v2s[, direction] = 2*midline - v2s[, direction];
		d = vnorm(vector_align(v1, v2s) - v1);
		all(d < eps)
	});
	spairs = epairs[s, ];
	dimnames(spairs)[[2]] = paste('symm', 1:2, sep = '');
	nonpaired = setdiff(1:nrow(elements), unique(as.vector(spairs)));
	if (!selfSymmetry) return(list(pairs = spairs, direction = direction, nonpaired = nonpaired));
	r = list(pairs = spairs, direction = direction, selfSymmetryRef = midline,
		dim = dim(graph)[2], nonpaired = nonpaired);
	r
}

#
#	<p> areas/triangulation
#
#	In order to get full symmetry in the triangulation the so called lower paritiion (i.e. the left side)
#	of the face is triangulated. This corresponds one-by-one to a triangulation fo the upper part
#	creating a fully symmetric triangulation

# get a graph segment based on a direction and a percentage of extend to be included
graphPartition = function(graph, direction = 1, threshold = .5, eps = 1e-5) {
	rng = range(graph[, direction]);
	thresholdAbs = rng[1] + (rng[2] - rng[1]) * threshold;
	graphPart = which(graph[, direction] <= thresholdAbs + eps);
	graph[graphPart, ]
}

matrix2dict = function(m)listKeyValue(m[, 1], m[, 2])
matrix2dictSymm = function(m) {
	m0 = unique(rbind(m, cbind(m[, 2], m[, 1])));
	listKeyValue(m0[, 1], m0[, 2])
}

# graph assumed to be symmterized
# use two dimensions to compute triangulation
delaunaySymm = function(graph, direction = 1, useDimensions = 1:2) {
	graph = graph[, useDimensions, drop = F];
	nodeSymm = graphSymmetry(graph, as.matrix(1:nrow(graph)), selfSymmetry = TRUE)$pairs;
	nodeSymmNs = matrix(rownames(graph)[nodeSymm], ncol = 2);	# by name
	nodeSymmDict = matrix2dictSymm(nodeSymmNs);
	# <p> using graphSymmetry has no canonical preference with respect to position (e.g. left/right)
	# therefore graphPartition is used
	# part = graph[nodeSymm[, 1], ];
	part = graphPartition(graph, direction);
	# <p> triangulate "left" part (i.e. everything below .5 range)
	triangles = delaunayn(part);
	trianglesNs = matrix(rownames(part)[triangles], ncol = 3);
	# <p> add symmetric triangles (in terms of node names)
	# <!> nodes w/o symmetry
	# <A> uniqueness
	trianglesSymmNs = nodeSymmDict[trianglesNs];
	# <p> convert back to indeces
	trianglesUnique = matrix(which.indeces(as.character(trianglesNs), rownames(graph)), ncol = 3);
	trianglesSymm = matrix(which.indeces(as.character(trianglesSymmNs), rownames(graph)), ncol = 3);
	trianglesAll = rbind(trianglesUnique, trianglesSymm);
	dimnames(trianglesAll) = list(NULL, c('a', 'b', 'c'));
	r = list(
		triangles = trianglesAll,
		symm = cbind(1:nrow(trianglesUnique), nrow(trianglesUnique) + (1:nrow(trianglesUnique)))
	);
	r
}
graphSymmetryTriangles = function(graph, direction = 1)list(pairs = delaunaySymm(graph, direction = 1)$symm)

#
#	<p> abstract function
#

#
#	<p> graph pre-processing
#

# should work for 2D, 3D, and kD (untested) data
# assume: coord 1 is X-coordinate
symmetrizeGraph = function(graph, nodeSymmetries = symmetry_visigenStd) {
	# <p> indeces and global midline
	nodes = dimnames(graph)[[1]];
	# indeces of nodes on the midline
	midline = which.indeces(
		nodeSymmetries[which(apply(nodeSymmetries, 1, function(nds)nds[1] == nds[2])), 1], nodes
	);
	# pairs of indeces for symmetric landmarks as matrix
	symm = matrix(which.indeces(
		as.vector(nodeSymmetries[which(apply(nodeSymmetries, 1, function(nds)nds[1] != nds[2])), ]), nodes
	), ncol = 2);
	# mean of midline landmarks
	graph_midline = mean(graph[midline, 1]);
	graphs = graph;	# symmetrized graph

	# <p> adjust midline to "midline-center"
	graphs[midline, 1] = graph_midline;

	# <p> adjust symmetric nodes
	graphsn = apply(symm, 1, function(sn) {
		symmprep = rbind(
			c(graph_midline - graph[sn[1], 1], graph[sn[1], -1]),
			c(graph[sn[2], 1] - graph_midline, graph[sn[2], -1])
		);
		means = apply(symmprep, 2, mean);
		symmetrized = c(graph_midline - means[1], means[-1], means[1] + graph_midline, means[-1]);
		symmetrized
	});
	graphs[as.vector(t(symm)), ] = matrix(as.vector(graphsn), byrow = T, ncol = dim(graph)[2]);

	r = list(graphs = graphs, grapha = graph - graphs, midline_x = graph_midline);
	r
}

# if a certain structure (e.g. triangles) produces several features a corresponding symmetry can be computed
#	automatically <i>
graphSymmetryExpand = function(symm, N = 3) {
	p = symm$pair;
	pairsNew = matrix.intercalate((p - 1) * N + 1, (p - 1) * N + 2, (p - 1) * N + 3);
	c(list(pairs = pairsNew), List_(symm, min_ = 'pairs'))
}

#
#	<p> exposed interface functions
#

graphSymmetries = function(graph, direction = 1) {
	# <p> nodes
	nodes = as.matrix(1:nrow(graph));
	nodeSymm = graphSymmetry(graph, nodes, selfSymmetry = TRUE);
	# <p> distances
	distances = do.call(rbind, lapply(as.list(set_combn(nodes, 2L)), unlist));
	distSymm = graphSymmetry(graph, distances);
	# <p> triangles
	triangles = graphSymmetryTriangles(graph, direction = direction);
	r = list(node = nodeSymm, distance = distSymm, triangle = triangles);
	r
}

#
#	<p> general helper functions
#

symmetrizedAverageGraph = function(coords, flip = FALSE, extend = NULL, ...) {
	grMn = meanGraph(coords);
	grS = symmetrizeGraph(grMn, ...)$graphs;
	if (flip) grS[, 2] = extend - grS[, 2];
	grS
}