#
#	facePreprocess.R
#Mon Aug 31 12:18:52 CEST 2015

convertImage = function(input, output, quality = 80L) {
	cmd = Sprintf("convert -quality %{quality}s %{input}Q %{output}Q");
	System(cmd, 3);
}

convertImages = function(path, output, extension = 'jpg', quality = 80L) {
	files = list_files_with_exts(path, extension);
	Dir.create(output, recursive = TRUE);
	r = lapply(files, function(path){
		file = splitPath(path)$file;
		convertImage(path, Sprintf("%{output}s/%{file}s"), quality);
	});
	r
}

defaultExts = c('jpg', 'jpeg', 'JPG', 'JPEG');
listImages = function(path, extension = defaultExts)list_files_with_exts(path, extension);
readImages = function(path, extension = defaultExts, ...)lapply(listImages(path, extension), function(file)readImage(file));
images2grey = function(is)lapply(is, function(i)channel(i, 'gray'))

readCoordinatesFromFile = function(path, path2metaRegex = NULL, reader = readCoordinateFile_visigen) {
	d0 = reader(path);
	if (!('id' %in% names(d0))) {
		name = splitPath(path)$base;
		d0$id = name;
	}
	if (is.null(path2metaRegex)) {
		d0$type = NA;
	} else {
		d0$type = as.character(fetchRegexpr(regex, name, captures = T));
		if (length(type) == 0) type = NA;
	}
	d0
}

dataCheckNodeAlignment = function(data, stopIfUnaligned = TRUE) {
	nodes = sapply(unique(data$id), function(id)data$node[data$id == id]);
	unaligned = any(!apply(nodes, 1, function(e)all(e == e[1])));
	if (stopIfUnaligned && unaligned) stop("Nodes in data set not in same order");
	!unaligned
}

dataAlignNodes = function(data) {
	d = data;
	ids = unique(d$id);
	nodesRef = d$node[d$id == ids[1]];
	r1 = lapply(ids, function(id) {
		nodesThis = d$node[d$id == id];
		d[d$id == id, ][order_align(nodesRef, nodesThis), ];
	});
	d = do.call(rbind, r1);
	dataCheckNodeAlignment(d);
	list(coords = d, nodes = nodesRef)
}

#' @arg path Path to folder from which files with extension `pos` are parsed
#' @arg path2metaRegex Provide regular expression that is applied to the path and caputres of which
#'   are used as meta-information returned in the type column of the result data frame
readCoordinateData = function(path, path2metaRegex = NULL, reader = readCoordinateFile_visigen) {
	files = list_files_with_exts(path, 'pos');
	r0 = lapply(files, readCoordinatesFromFile, reader = reader, path2metaRegex = path2metaRegex);
	r1 = dataAlignNodes(do.call(rbind, r0));
	# <p> return
	d = data.frame.types(r1$coords, factor = c('id', 'type', 'node'));
	r = list(coord = d, nodes = r1$nodes);
	r
}


table2array = function(d, id = 'id', coords = c('x', 'y'), node = 'node') {
	ids = levels(d[[id]]);	# individuals
	nodes = as.character(d[[node]])[1:(nrow(d)/length(ids))];
	# expected order of nodes
	nodeRows = rep(rep.each(nodes, length(coords)), length(ids));
	d1 = reshape.long(d, coords);
	if (any(d1$node != nodeRows)) stop('order of nodes inconsistent bewteen individuals');
	a = array(d1$value, c(length(coords), length(nodes), length(ids)));
	a = aperm(a, c(2, 1, 3));
	dimnames(a)[[1]] = nodes;
	dimnames(a)[[3]] = levels(d$id);
	a
}

array2table = function(a, ids, nodes, coords = c('x', 'y')) {
	Nr = dim(a)[3];

	dfs = lapply(1:Nr, function(i) {
		dfCoords = data.frame.types(a[,,i], names = coords);
		d1 = reshape.wide(data.frame(id = ids[i], node = nodes, dfCoords), ids = 'id', vars = 'node',
			blockVars = T, reverseNames = T);
	});
	d1 = rbindDataFrames(dfs);
	d1
}


coords2distances = function(d, vertices, patterns = c('%s.x', '%s.y')) {
	# <p> find columns
	ns = names(d);
	coords = nlapply(vertices, function(n) {
		which.indeces(c(sprintf(patterns[1], n), sprintf(patterns[2], n)), ns)
	});

	# <A> identify based on x-coordinate, assume y-coordinate to exist
	modelList = list(x = vertices, y = vertices);
	r = iterateModels(modelList, .constraint = function(x, y)(x < y), function(i, x, y) {
		#Log(sprintf('Vertices (%s, %s)', x, y), 4);
		dims = length(coords[[x]]);
		dist = sapply(1:dims, function(dim)((d[, coords[[x]][dim]] - d[, coords[[y]][dim]])^2));
		distEuclidian = apply(dist, 1, function(r)sqrt(sum(r)));
		distEuclidian
	}, lapply__ = lapply);
	distanceNames = apply(r$models_symbolic, 1, function(r)sprintf('d_%s', join(r, sep = '_')));
	distances = data.frame.types(sapply(r$results, identity), names = distanceNames);
	distances
}

#
#	<p> procrustes
#

performProcrustes = function(a) {
	pa = procGPA(a, affine = T);
	#dr = pa$rotated;	# rotated data
	pa
}

procrustesTransformations = function(pa) {
}

#
#	<p> transformations
#

# euclidian norm
vector.norm = function(v)sqrt(v %*% v)
# homogeneous coordinates
hom = function(d)cbind(d, 1)
# homogeneous coordinates -> euclidian coordinates
homI = function(d)d[, 1:2]

# X: data matrix
# Xt: transformed data
# X' T = Xr' (priming -> homogenous coordinates) => T = X'^- Xr'
affineBetween = function(X, Xt) {
	ginv(hom(X)) %*% hom(Xt)
}

# http://math.stackexchange.com/questions/237369/given-this-transformation-matrix-how-do-i-decompose-it-into-translation-rotati
# http://math.stackexchange.com/questions/13150/extracting-rotation-scale-values-from-2d-transformation-matrix/13165#13165
# @result list to be read as follows: rot -> scale -> translate or scale -> rot -> translate

affine2components = function(m) {
	m = t(m);
	sx = vector.norm(m[1, 1:2]) * sign(m[1, 1]);
	sy = vector.norm(m[2, 1:2]) * sign(m[2, 2]);
	phi = atan2(-m[1, 2], m[1, 1]);
	phi1 = atan2(m[2, 1], m[2, 2]);
	aff = list(translate = m[1:2, 3], scale = c(sx, sy), rot = c(phi, phi1));
	aff
}
components2affine = function(aff) with(aff, {
	m = matrix(c(
		scale[1] * cos(rot[1]), -scale[1] * sin(rot[1]), translate[1],
		scale[2] * sin(rot[1]), scale[2] * cos(rot[1]), translate[2],
		0, 0, 1
	), nrow = 3, ncol = 3, byrow = T);
	t(m)
})
affineCreate = function(scale = c(1, 1), rot = c(0, 0), translate = c(0, 0)) {
	components2affine(list(translate = translate, scale = scale, rot = rot));
}

affRotate = function(data, aff) with(aff, {
	m = matrix(c(
		cos(rot[1]), -sin(rot[1]),
		sin(rot[1]), cos(rot[1])
	), nrow = 2, ncol = 2, byrow = T);
	data %*% t(m)
	
})

affScale = function(data, aff) with(aff, {
	m = matrix(c(
		scale[1], 0,
		0, scale[2]
	), nrow = 2, ncol = 2, byrow = T);
	data %*% t(m)
})

affTranslate = function(data, aff) with(aff, {
	cbind(data[, 1, drop = F] + translate[1], data[, 2, drop = F] + translate[2])
})

# same as data %*% components2affine(aff)
affPerform = function(data, aff) {
	affTranslate(affScale(affRotate(data, aff), aff), aff);
}

boundingBox = function(a) {
	list(
		xlim = c(min(a[, 1, ]), max(a[, 1, ])),
		ylim = c(min(a[, 2, ]), max(a[, 2, ]))
	)
}
boundingBoxSquared = function(bb, dim = 256, border = .1) {
	# do not center
	ext = max(abs(bb$xlim), abs(bb$ylim));
	scale = dim * (1 - border) / ext;
	scale
}


#
#	GraphicsMagick
#

pictureDimensions = function(path) {
	cmd = Sprintf("identify -format '%%w %%h' %{path}Q");
	as.integer(splitString(' ', System(cmd, 5, return.output = T)$output))
}

#
#	transform image to reflect procrustes rotation
#	assumption: graph coordinates in pixels
#	apply transform -> centered coordinates
#	add scaling to for target coordinates
#	add translation to center to target coordinates
#	compute borders to result in a square image
#	generate imagemagick/graphicsmagick command
#

pictureMap2output = function(aff, outputDim, scale = 1) {
	# rescale to fit target size
	aff = aff %*% affineCreate(scale = rep(scale, 2));
	aff[3, 1:2] = aff[3, 1:2] + outputDim/2;
	aff
}

#' @arg recenter: assume aff is centered around the origin, recenter to middle coordinate of the picture
pictureTransform = function(path, aff, outputDir, outputDim = pictureDimensions(path),
	extraArgs = '-type Grayscale', recenter = T, scale = 1, quality = .9) {
	Dir.create(outputDir, recursive = TRUE);
	file = splitPath(path)$file;

	aff = pictureMap2output(aff, outputDim, scale);
	a = sapply(as.vector(t(aff[, 1:2])), function(v)Sprintf('%{v}.2f'));
	aS = join(a, ',');
	draw = Sprintf("affine %{aS}s image over 0,0 0,0 %{path}Q");
	cmd = Sprintf('convert -background grey -size %{w}sx%{h}s xc:snow3 -draw %{draw}Q %{extraArgs}s -quality %{q}.0f %{outputDir}Q/%{file}Q',
		w = outputDim[1], h = outputDim[2], q = 1e2*quality);
	System(cmd, 2);
}

pictureAnnotate = function(path, graph, outputDir, scale = 1, borders = rep(0, 4)) {
	Dir.create(outputDir);
	file = splitPath(path)$file;
	graphS = paste(
		c('', apply(graph, 1, function(r)join(
				c(join(sprintf('%.0f', r), ','), join(sprintf('%.0f', r + 2), ','))
			)
		))
	, collapse = ' circle ');
	cmd = Sprintf('convert %{path}Q -draw %{graphS}Q %{outputDir}Q/%{file}Q');
	System(cmd, 2);
}

#
#	<p> align batch of pictures to common image positions
#


boundingBox = function(a) {
	list(
		xlim = c(min(a[, 1, ]), max(a[, 1, ])),
		ylim = c(min(a[, 2, ]), max(a[, 2, ]))
	)
}
maxAbs = function(e, ...)max(abs(e), ...)
ext = function(e)(e[2] - e[1])
scaleFromBoundingBox = function(bb, dimTarget, margin = .05) {
	dimMax = sapply(bb, maxAbs) * 2;
	scaleD = (dimTarget * (1 - margin)) / dimMax;
	# rescale to uniformly fit into target size
	scale = min(scaleD);
	scale
}

procrustesTransformImages = function(path, output,
	outputDirRaw = con(output, '/00_raw'),
	outputDirAnnotation = con(output, '/01_annotated'),
	outputDirTransform = con(output, '/02_transformed'),
	outputDirAnnotationTransf = con(output, '/03_transformed_annotated'),
	readCoordinates = readCoordinateData,
	annotate = T, dimTarget = rep(512, 2), margin = .05, ...) {

	convertImages(path, outputDirRaw);
	d = readCoordinates(path, ...);
	a = table2array(d$coord);
	pa = performProcrustes(a);
	scale = scaleFromBoundingBox(boundingBox(pa$rotated), dimTarget, margin);

	ip = listImages(path);
	coords = ilapply(ip, function(path, i) {
		file = splitPath(path)$file;
		pdim = pictureDimensions(path);
		if (annotate) pictureAnnotate(path, a[, , i], outputDirAnnotation);

		aff = affineBetween(a[, , i], pa$rotated[, , i]);
		rot = affine2components(aff);
		#rot$rot[1] = 20/360 * (2 * pi);
		aff = components2affine(rot);
		pictureTransform(path, aff, outputDirTransform, outputDim = dimTarget, scale = scale);

		affImg = pictureMap2output(aff, dimTarget, scale);
		paCoords = homI(hom(a[, , i]) %*% affImg);
		if (annotate) pictureAnnotate(
			Sprintf('%{outputDirTransform}s/%{file}s'), paCoords, outputDirAnnotationTransf);
		paCoords
	});
	coords = array(unlist(coords), dim = c(dim(coords[[1]]), length(coords)));
	dimnames(coords) = dimnames(a);
	ids = sapply(ip, function(e)splitPath(e)$base);
	r = list(id = ids,
		coords = coords,
		# paths
		images = ip, input = path, output = output,
		outputDirRaw = outputDirRaw, outputDirAnnotation = outputDirAnnotation,
		outputDirTransform = outputDirTransform, outputDirAnnotationTransf = outputDirAnnotationTransf
	);
	r
}

prepareAveraging = function(collection, dataArray,
	preparer = prepareAveraging_ini, output = con(collection$output, '/04_averaging_prepare'), ...) {
	myoutput = output;

	with(collection, {
	ilapply(id, function(id, i)
		preparer(Sprintf('%{outputDirTransform}s/%{id}s.jpg'), dataArray[, , i], output = myoutput, ...));
	r = c(collection, list(outputDirAveragingPrepare = myoutput));
	r
})}

averageGroups = function(collection, groups,
	output = con(collection$output, '/05_averages'), averager = averageGraphs_ini) {
	ls = levels(as.factor(groups));
	sapply(ls, function(level) {
		sel = collection$id[which(groups == level)];
		averager(collection, level, sel, output);
	});
	r = c(collection, list(outputAverages = output));
	r
}

#
#	<p> graph pre-processing
#

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
			c(graph_midline - graph[sn[1], 1], graph[sn[1], 2]),
			c(graph[sn[2], 1] - graph_midline, graph[sn[2], 2])
		);
		means = apply(symmprep, 2, mean);
		symmetrized = c(graph_midline - means[1], means[2], means[1] + graph_midline, means[2]);
		symmetrized		
	});
	graphs[as.vector(t(symm)), ] = matrix(as.vector(graphsn), byrow = T, ncol = 2);

	r = list(graphs = graphs, grapha = graph - graphs, midline_x = graph_midline);
	r
}
