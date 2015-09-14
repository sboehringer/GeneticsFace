#
#	faceData.R
#Thu Sep  3 14:33:52 2015

#
#	<p> visigen file format
#

readCoordinateFile_visigen = function(path) {
	readTable(sprintf('[SEP=T,HEADER=F,NAMES=node;x;y]:%s', path));
}

#
#	<p> ini xml files
#

#
# adapted form stackoverflow
#

require('XML');
list2xmlRaw = function(node, l) {
	for(i in seq_along(l)) {
		child = newXMLNode(names(l)[i], parent = node);
		if (typeof(l[[i]]) == "list") list2xmlRaw(child, l[[i]]) else xmlValue(child) = l[[i]];
	} 
}
list2xml = function(l, rootName = 'root', asDocument = FALSE) {
	document = if (asDocument) newXMLDoc() else NULL;
	root = newXMLNode(rootName, doc = document);
	list2xmlRaw(root, l);
	r = if (asDocument) document else root;
	r
}

xmlTemplateModelGraph = list(
	GraphType = 'MODEL_GRAPH',
	NumberOfNodes = NA,
	NumberOfEdges = 0,
	Trafo = list(
		TrafoType = 'GABOR',
		NumberOfLevels = 5,
		NumberOfDirections = 8,
		Maximum = 1.5708,
		Factor = 0.707107,
		Sigma = 6.28319,
		DcFree = 1,
		SpaceAccuracy = 3
	),
	#Topology = list(Edge = list(SourceNode = 1, TargetNode = 2))
	Topology = list(),
	#EdgeList = list(Edge = list(Index = 1, EdgeInfoType = 'MODEL_EDGE_INFO'))
	EdgeList = list()
);

#' @par: graph array containing coordinates
coordinateStructure_ini = function(graph, rootName = 'ModelGraph', asDocument = TRUE) {
	graphStructureL = xmlTemplateModelGraph;
	graphL = ilapply(apply(graph, 1, as.list), function(coords, i) {
		list(
			Index = i - 1, ToBeUsed = 1, NodeInfoType = 'MODEL_NODE_INFO',
			Position = list(Element = avu(coords[1]), Element = avu(coords[2]))
		)
	});
	names(graphL) = rep('Node', length(graphL));
	graphStructureL$NumberOfNodes = nrow(graph);
	graphStructureL$NodeList = graphL;
	r = list2xml(graphStructureL, rootName = rootName, asDocument = asDocument);
	r
}

writeCoordinateFile_ini = function(graph, output = NULL) {
	header = '<?xml version="1.0"?>';
	ini = coordinateStructure_ini(graph);
	if (!is.null(output)) Dir.create(output, recursive = TRUE, treatPathAsFile = TRUE);
	saveXML(ini, output);
}

readCoordinateFile_ini = function(path, input) {
	xml = xmlParseDoc(file);
	coords = as.integer(
		xpathSApply(xml,'//ModelGraph/NodeList/Node/Position/Element', function(i)xmlValue(i))
	);
}

prepareAveraging_ini = function(pathImg, coords, output,
	pathImgConverted = Sprintf('%{output}s/%{base}s.tif', base = splitPath(pathImg)$base),
	pathGraph = Sprintf('%{output}s/%{base}s.xml', base = splitPath(pathImg)$base)) {

	writeCoordinateFile_ini(coords, output = pathGraph);
	convertImage(pathImg, pathImgConverted);

	cmd = Sprintf(con('ini_X64 --extract-jets ',
		'--input-graph %{pathGraph}Q ',
		'--input-image %{pathImgConverted}Q ',
		'--output-graph %{pathGraph}Q'));
	Log(cmd, 2);
	System(cmd, 2);
}