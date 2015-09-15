#
#	facePlotting.R
#Tue Sep 15 17:41:24 2015


# annotate("text", x = 2:5, y = 25, label = "Some text")
plotGraph = function(coords, nodes = NULL, extend = 512, flip = TRUE) {
	if (flip) coords[, 2] = extend - coords[, 2];
	p = ggplot(Df_(coords, names = c('x', 'y')), aes(x = x, y = y)) +
		geom_point() +
		scale_x_continuous(limits = c(0, extend)) + scale_y_continuous(limits = c(0, extend)) +
		theme_bw();
	if (!is.null(nodes)) {
		p = p + annotate('text', x = coords[, 1], y = coords[, 2], label = nodes, colour = 'blue', size = 4);
	}
	p
}

plotMeanGraph = function(coords, nodes = NULL, extend = 512, flip = TRUE) {
	coordsMean = apply(coords, 1:2, mean);
	#tri = delaunayn(coordsMean);
	plotGraph(coordsMean, nodes, extend, flip);
}
