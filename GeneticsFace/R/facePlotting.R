#
#	facePlotting.R
#Tue Sep 15 17:41:24 2015

library('ggplot2');

# annotate("text", x = 2:5, y = 25, label = "Some text")
plotGraph = function(coords, labels = FALSE, extend = 512, flip = TRUE, labelSize = 4) {
	if (flip) coords[, 2] = extend - coords[, 2];
	p = ggplot(Df_(coords, names = c('x', 'y')), aes(x = x, y = y)) +
		geom_point() +
		scale_x_continuous(limits = c(0, extend)) + scale_y_continuous(limits = c(0, extend)) +
		theme_bw();
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
