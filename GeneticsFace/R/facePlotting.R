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
