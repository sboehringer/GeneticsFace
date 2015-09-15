#
#	faceData.R
#Thu Sep  3 14:33:52 2015

#
#	<p> visigen file format
#

symmetry_visigenStd = matrix(c(
	'k1', 'k1',
	'k2', 'k3',
	'k4', 'k5',
	'k6', 'k10',
	'k7', 'k9',
	'k8', 'k8',
	'k11', 'k14',
	'k12', 'k13',
	'k16', 'k20',
	'k22', 'k29',
	'k24', 'k27',
	'k17', 'k19',
	'k23', 'k28',
	'k15', 'k21',
	'k25', 'k26',
	'k30', 'k33',
	'k18', 'k18',
	'k31', 'k32',
	'k34', 'k38',
	'k35', 'k37',
	'k36', 'k36',
	'k39', 'k42',
	'k40', 'k41',
	'k43', 'k47',
	'k45', 'k45',
	'k44', 'k46',
	'k48', 'k48',
	'k51', 'k51',
	'k49', 'k53',
	'k50', 'k52',
	'k54', 'k54'
), byrow = T, ncol = 2);
Nsymmetry_visigenStd = length(unique(as.vector(symmetry_visigenStd)));

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

	cmd = Sprintf(con('ini_X64 --extract-jets --gabor-trafo-type Gray ',
		'--input-graph %{pathGraph}Q ',
		'--input-image %{pathImgConverted}Q ',
		'--output-graph %{pathGraph}Q'));
	System(cmd, 2);
}

# size_x/.75
averageGraphs_ini = function(collection, name, sel, output, size_x = 512, size_y = size_x, Niter = 0L) {
	Dir.create(output, recursive = T);
	# <p> average
	averagingInput = collection$outputDirAveragingPrepare;
	averageGraph = Sprintf('%{output}Q/%{name}Q.xml');
	averageGraphs = join(sapply(sel, function(id)Sprintf('%{averagingInput}Q/%{id}Q.xml')), "\n");
	averageGraphList = tempfile();
	writeFile(averageGraphList, averageGraphs);
	cmd = Sprintf(con('ini_X64 --average-graphs ',
		'--input-graph-list %{averageGraphList}Q ',
		'--output-graph %{averageGraph}Q'));
	System(cmd, 2);

	# <p> reconstruct
	cmd = Sprintf(con('ini_X64 --reconstruct --reconstruction-iterations %{Niter}s ',
		'--reconstructed-image-size %{size_x}s,%{size_y}s ',
		'--input-graph %{averageGraph}Q ',
		'--output-image %{output}Q/%{name}Q.tif '));
	System(cmd, 2);
	convertImage(Sprintf('%{output}Q/%{name}Q.tif'), Sprintf('%{output}Q/%{name}Q.jpg'));
}
