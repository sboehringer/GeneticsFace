#
#	facePlotting.R
#Tue Sep 15 17:41:24 2015

library('ggplot2');

triangulation2segments = function(graph, tri, extend = 512, doFlipY = TRUE) {
	g = graph;
	triSegments = apply(tri, 1, function(triangle) {
		s = sapply(seq_along(triangle), function(i)c(g[triangle[i], ], g[triangle[(i %% 3) + 1], ]));
		s
	});
	dS = Df_(matrix(as.vector(triSegments), byrow = T, ncol = 4), names = c('x', 'y', 'xend', 'yend'));
	if (doFlipY) dS[, c('y', 'yend')] = extend - dS[, c('y', 'yend')];
	dS
}

# annotate("text", x = 2:5, y = 25, label = "Some text")
plotGraph = function(coords, labels = FALSE, extend = 512, flip = TRUE, labelSize = 4, segments = NULL) {
	if (flip) coords[, 2] = extend - coords[, 2];
	p = ggplot(Df_(coords, names = c('x', 'y')), aes(x = x, y = y)) +
		geom_point() +
		scale_x_continuous(limits = c(0, extend)) + scale_y_continuous(limits = c(0, extend)) +
		theme_bw();
	if (!is.null(segments))
		p = p + geom_segment(data = dS, size = .5, color = 'red', alpha = .2,
			mapping = aes(x = x, y = y, xend = xend, yend = yend)
	);
	if (labels) {
		nodes = dimnames(coords)[[1]];
		p = p + annotate('text', x = coords[, 1], y = coords[, 2], label = nodes,
			colour = 'blue', size = labelSize);
	}
	p
}

meanGraph = function(coords)apply(coords, 1:2, mean);
plotMeanGraph = function(coords, labels = FALSE, extend = 512, flip = TRUE) {
	plotGraph(meanGraph(coords), labels, extend, flip);
}

#
#	caricature functions
#

# cross product of two 2d vectors embeded into 3d
vectorCross = function(v1, v2)c(0, 0, det(cbind(v1, v2)))

# vectorized for angle-argument
circlePoint = function(center, radius, angle)
	t(center + matrix(radius * c(cos(angle), sin(angle)), nrow = 2, byrow = T))
circlePoints = function(center, radius, angles, Nlines = 100*(abs(angles[2] - angles[1]) / pi))
	circlePoint(center, radius, seq(angles[1], angles[2], length.out = Nlines))
circleSegment = function(center, radius, angles, Nlines = 100*(abs(angles[2] - angles[1]) / pi), ...) {
	#pts = circlePoints(center, radius, angles, Nlines = 3);
	#xsplineGrob(pts[, 1], pts[, 2], shape = -1);
	pts = circlePoints(center, radius, angles, Nlines);
	linesGrob(pts[, 1], pts[, 2], ...);
}

Grid.draw = function(e, ...) {
	if (!any(class(e) == 'grob')) lapply(e, Grid.draw, ...) else grid.draw(e, ...)
}

graphAngleSegments = function(graph, triangles, radius = -.5, symmetries = NULL) {
	as = matrix(extractAngle(graph, triangles, symmetries)$feature, ncol = 3, byrow = T);
	asE = matrix(extractAngleExternalX(graph, triangles, symmetries)$feature, ncol = 3, byrow = T);
	as = as * triangleOrientation(graph, triangles);
	dist = matrix(extractAngleMinDist(graph, triangles, symmetries)$feature, ncol = 3, byrow = T);
	r = lapply(seq(nrow(as)), function(i) {
		lapply(1:3, function(j) {
			cs = circleSegment(graph[triangles[i, j], ],
				if (radius < 0) dist[i, j] * (-radius) else radius, cumsum(c(asE[i, j], as[i, j])),
				default.units = 'points');
			cs
		})
	});
	r = unlist.n(r, 1);
	r
}
graphTriangulationSegments = function(graph, triangles, width = 1) {
	pts = graph[as.vector(t(cbind(triangles, triangles[, 1]))), ];
	id = rep.each(1:nrow(triangles), ncol(triangles) + 1);
	polylineGrob(pts[, 1], pts[, 2], id = id, default.units = 'points')
}
