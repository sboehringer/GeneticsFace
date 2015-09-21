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
	midline = graphmid(graph);
	Els = 1:nrow(elements);
	epairs = do.call(rbind, lapply(as.list(set_combn(Els, 2L)), unlist));
	# self-symmetry
	if (selfSymmetry) epairs = rbind(cbind(Els, Els), epairs);
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
	spairs
}

graphSymmetries = function(graph, direction = 1) {
	coords = 1:nrow(graph);
	
	tri = delaunayn(graphs);
}

