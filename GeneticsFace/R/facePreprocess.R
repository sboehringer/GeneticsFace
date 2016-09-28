#
#	facePreprocess.R
#Mon Aug 31 12:18:52 CEST 2015

#
#	<p> Coordinates
#
#	Coordinates are handled in the first quadrant. Reader should therefore convert to 1st quadrant
#	coordinates. EBimage handles images in the 4th quadrant, requiring coordinate transformation.

defaultExts = c('jpg', 'jpeg', 'JPG', 'JPEG', 'png', 'PNG');

removeExclusions = function(files, removeIds = NULL) {
	if (!is.null(removeIds)) {
		remove = sapply(files, function(e)splitPath(e)$base) %in% removeIds;
		Logs("Removed %{N}d data sets.", N = sum(remove), logLevel = 5);
		files = files[!remove];
	}
	files
}

applyPathesFiles = function(pathes, f, ..., extension, noMerge = FALSE, removeIds = NULL) {
	r = lapply(pathes, function(p) {
		files = removeExclusions(list_files_with_exts(p$path, extension), removeIds);
		lapply(files, function(path)f(path, ...))
	});
	if (noMerge) return(r);
	r1 = merge.lists(r, listOfLists = T, concat = T, useIndeces = T);
	r1
}

applyPathes = function(pathes, f, ..., noMerge = FALSE) {
	r = lapply(pathes, function(p)f(p$path, ...))
	if (noMerge) return(r);
	r1 = merge.lists(r, listOfLists = T, concat = T, useIndeces = T);
	r1
}

convertImageRaw = function(input, output, quality = 80L) {
	cmd = Sprintf("convert -quality %{quality}s %{input}Q %{output}Q");
	System(cmd, 5);
}

convertImage = function(path, output, outputFormat = 'jpg', quality = 95L) {
	base = splitPath(path)$base;
	outputPath = Sprintf("%{output}s/%{base}s.%{outputFormat}s");
	convertImageRaw(path, outputPath, quality);
	outputPath
}

convertImages = function(path, output, extension = defaultExts, outputFormat = 'jpg', quality = 95L) {
	files = list_files_with_exts(path, extension);
	Dir.create(output, recursive = TRUE);
	r = lapply(files, convertImage, output = output, outputFormat = outputFormat, quality = quality);
	r
}

listImages = function(path, extension = defaultExts)list_files_with_exts(path, extension);
readImages = function(path, extension = defaultExts, ...)lapply(listImages(path, extension), function(file) readImage(file));
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
		type = unlist(matchRegex(path2metaRegex, name)$capture);
		#d0$type = as.character(fetchRegexpr(regex, name, captures = T));
		d0$type = if (length(type) == 0) NA else type;
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

#' @arg path Path to folder from which files with extension coordinateExt are parsed
#' @arg path2metaRegex Provide regular expression that is applied to the path and caputres of which
#'   are used as meta-information returned in the type column of the result data frame
#' @arg coordinateExt file extension of coordinate files. Defaults to 'pos'.
readCoordinateData = function(path, path2metaRegex = NULL, reader = readCoordinateFile_visigen,
	imageExts = defaultExts, type = NULL, removeIds = NULL, coordinateExt = 'pos') {
	files = list_files_with_exts(path, coordinateExt);
	if (!all(
		sort(sapply(list_files_with_exts(path, imageExts), function(f)splitPath(f)$base)) ==
		sort(sapply(files, function(f)splitPath(f)$base)))) {
		stop(Sprintf('Mismatch coordinate, image files in: %{path}s'));
	}
	files = removeExclusions(files, removeIds);
	r0 = lapply(files, readCoordinatesFromFile, reader = reader, path2metaRegex = path2metaRegex);
	r1 = dataAlignNodes(do.call(rbind, r0));
	# <p> return
	d = Df_(r1$coords, as_factor = c('id', 'type', 'node'));
	if (!is.null(type)) d$type = type;
	r = list(coord = d, nodes = r1$nodes);
	r
}
readCoordinatesDefaults = list(flip = FALSE, flipExtend = 512, reader = 'readCoordinateFile_visigen');
readCoordinateDataFromPathes = function(pathList,
	imageExts = defaultExts, removeIds = NULL, reader = 'readCoordinateFile_visigen', coordinateExt = 'pos') {
	defaults = merge.lists(readCoordinatesDefaults, list(coordinateExt = coordinateExt, reader = reader));
	d = lapply(pathList, function(p) {
		p = merge.lists(defaults, p);
		if (is.character(p$reader)) p$reader = get(p$reader);
		d = readCoordinateData(p$path, type = p$type,
			path2metaRegex = p$groupRegex, imageExts = imageExts, removeIds = removeIds,
			reader = p$reader, coordinateExt = p$coordinateExt);
		if (p$flip) d$coord$y = p$flipExtend - d$coord$y;
		d
	});
	nodes = do.call(cbind, lapply(list.kp(d, 'nodes'), sort));
	if (!all(apply(nodes, 1, same.vector))) {
		print(nodes);
		stop('Sets of data with incompatible nodes');
	}
	r1 = dataAlignNodes(do.call(rbind, list.kp(d, 'coord')));
	# <p> return
	d1 = Df_(r1$coords, as_factor = c('id', 'type', 'node'));
	r = list(coord = d1, nodes = r1$nodes);
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
	t(ginv(hom(X)) %*% hom(Xt))
}

# http://math.stackexchange.com/questions/237369/given-this-transformation-matrix-how-do-i-decompose-it-into-translation-rotati
# http://math.stackexchange.com/questions/13150/extracting-rotation-scale-values-from-2d-transformation-matrix/13165#13165
# @result list to be read as follows: rot -> scale -> translate or scale -> rot -> translate

# assume coords to be row-wise
affApply = function(aff, coords)homI(t(aff %*% t(hom(coords))))

affine2components = function(m) {
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
	m
})
affineCreate = function(scale = c(1, 1), rot = c(0, 0), translate = c(0, 0)) {
	components2affine(list(translate = translate, scale = scale, rot = rot));
}

affRotate = function(data, aff) with(aff, {
	m = matrix(c(
		cos(rot[1]), -sin(rot[1]),
		sin(rot[1]), cos(rot[1])
	), nrow = 2, ncol = 2, byrow = T);
	affApply(m, data)
})

affScale = function(data, aff) with(aff, {
	m = matrix(c(
		scale[1], 0,
		0, scale[2]
	), nrow = 2, ncol = 2, byrow = T);
	affApply(m, data)
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
#	multiplication of is right to left (right-most transformation applied first)
#

pictureMap2output = function(aff, outputDim, scale = 1) {
	# rescale to fit target size
	affS = affineCreate(scale = rep(scale, 2));
	# final translation is on original scale
	#affT = affineCreate(translate = rep(outputDim/2, 2)) %*% ginv(affS);
	affT = affineCreate(translate = rep(outputDim/2, 2));
	#aff[3, 1:2] = aff[3, 1:2] + outputDim/2;
	#aff[3, 1:2] = affC$translate + outputDim/2;
	affR = affT %*% affS %*% aff;
	affR
}

#' @arg recenter: assume aff is centered around the origin, recenter to middle coordinate of the picture
pictureTransformImageMagick = function(path, aff, outputDir, outputDim = pictureDimensions(path),
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

# aff: 3x3 matrix for right multiplication of column vectors
# picture coordinate system is in 4th quadrant -> conjugate with flipping operation
pictureTransform = function(image, aff, output = NULL, outputDim = 512, quality = 1, background = gray(.6)) {
	image = channel(image, 'gray');
	affC = affineCreate(scale = c(1, -1), translate = c(0, outputDim));
	affP = ginv(affC) %*% aff %*% affC;
	#afft = t(aff[1:2, ])
	afft = t(affP[1:2, ])
	imageTrans = affine(image, afft, output.dim = outputDim, bg.col = background);
	if (!is.null(output)) {
		Dir.create(splitPath(output)$dir, recursive = TRUE);
		writeImage(imageTrans, output, quality = quality * 1e2);
	}
	imageTrans
}

pictureAnnotateImageMagick = function(path, graph, outputDir, scale = 1, borders = rep(0, 4)) {
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

# <!> png makes white stripes
#install_github('oldGridExtra', 'ttriche')
# flipY depracted, should be handled when reading the file
pictureAnnotate = function(image, graph, output,
	colorNode = rgb(1, 0, 0), gridUnit = 'points', flipY = FALSE, nodeSize = 4, draw = FALSE,
	quality = 95L) {
	dimI = dim(image)[1:2];
	if (flipY) graph[, 2] = dimI[2] - graph[, 2];
	imgGrob = ebimageGrob(image, x = unit(dimI[1]/2, gridUnit), y = unit(dimI[2]/2, gridUnit));
	#imgGrob = ebimageGrob(image);
	nodes = pointsGrob(graph[, 1], graph[, 2],
			size = unit(nodeSize, gridUnit), pch = 21, gp = gpar(fill = colorNode), default.units = gridUnit);
	imageAnnotated = gList(imgGrob, nodes);
	if (draw) grid.draw(imageAnnotated);
	if (!is.null(output)) {
		Dir.create(splitPath(output)$dir, recursive = TRUE);
		plot_save(grid.draw(imageAnnotated), plot_path = output,
			width = valueU(dimI[1], 'points'), height = valueU(dimI[2], 'points'),
			unit = 'points', unit_out = 'points', quality = quality);
	}
	imageAnnotated
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
	readCoordinates = readCoordinateDataFromPathes,
	annotate = T, dimTarget = rep(512, 2), margin = .05, ..., backGround, removeIds = NULL) {

	Dir.create(output);
	if (!is.list(path)) path = list(list(path = path));
	d = readCoordinates(path, ..., removeIds = removeIds);
	Dir.create(outputDirRaw);
	pathesRaw = unlist(applyPathesFiles(path, convertImage, removeIds = removeIds, noMerge = T,
			output = outputDirRaw, extension = defaultExts));
	if (length(levels(d$coord$id)) != length(pathesRaw)) {
		stop(Sprintf("Coordinate data (N = %{Ncoord}d) does not match image files (N = %{Nfiles}d)",
			Ncoord = length(levels(d$coord$id)), Nfiles = length(pathesRaw)));
	}
	a = table2array(d$coord);
	pa = performProcrustes(a);
	scale = scaleFromBoundingBox(boundingBox(pa$rotated), dimTarget, margin);

	coords = ilapply(pathesRaw, function(path, i) {
		image = readImage(path);
		base = splitPath(path)$base;
		outputAnnot = Sprintf('%{outputDirAnnotation}s/%{base}s.jpg');
		if (annotate) pica = pictureAnnotate(image, a[, , i], outputAnnot);

		aff = affineBetween(a[, , i], pa$rotated[, , i]);
		#rot = affine2components(aff);
		#rot$rot[1] = 20/360 * (2 * pi);
		#aff = components2affine(rot);
		affImg = pictureMap2output(aff, dimTarget, scale);
		#affImg = affineCreate(translate = c(10, 10));
		outputTrans = Sprintf('%{outputDirTransform}s/%{base}s.jpg');
		imageTrans = pictureTransform(image, affImg, output = outputTrans, outputDim = dimTarget);

		paCoords = affApply(affImg, a[, , i]);
		outputTransAnnot = Sprintf('%{outputDirAnnotationTransf}s/%{base}s.jpg');
		if (annotate) picat = pictureAnnotate(imageTrans, paCoords, output = outputTransAnnot);
		paCoords
	});
	coords = array(unlist(coords), dim = c(dim(coords[[1]]), length(coords)));
	dimnames(coords) = dimnames(a);
	ids = sapply(pathesRaw, function(e)splitPath(e)$base);
	r = list(id = ids, group = unique(d$coord[, c('id', 'type')])$type, coords = coords,
		outputExtend = dimTarget,
		# paths
		images = pathesRaw, input = path, output = output,
		outputDirRaw = outputDirRaw, outputDirAnnotation = outputDirAnnotation,
		outputDirTransform = outputDirTransform, outputDirAnnotationTransf = outputDirAnnotationTransf
	);
	r
}

collectionRemoveIds = function(collection, ids) {
	i = which.indeces(ids, collection$id);
	r = merge.lists(collection, list(
		id = collection$id[-i], group = collection$group[-i], coords = collection$coords[,, -i],
		images = collection$images[-i]
	));
	r
}

prepareAveraging = function(collection, dataArray,
	preparer = prepareAveraging_ini, output = con(collection$output, '/04_averaging_prepare'), ...) {
	myoutput = output;

	with(collection, {
	ilapply(id, function(id, i)
		preparer(Sprintf('%{outputDirTransform}s/%{id}s.jpg'), dataArray[, , i], output = myoutput, ...,
			outputExtend = outputExtend[2]));
	r = c(collection, list(outputDirAveragingPrepare = myoutput));
	r
})}

averageGroups = function(collection, groups,
	output = con(collection$output, '/05_averages'), averager = averageGraphs_ini, ...) {
	ls = levels(as.factor(groups));
	sapply(ls, function(level) {
		sel = collection$id[which(groups == level)];
		averager(collection, level, sel, output, ...);
	});
	r = c(collection, list(outputAverages = output));
	r
}
