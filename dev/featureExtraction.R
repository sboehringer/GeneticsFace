#
#	expample illustrating feature extraction from classification model
#Tue Oct  6 14:16:24 CEST 2015

library('tools');
library('shapes');
library('EBImage');
library('geometry');
library('sets');
library('calibrate')
source('../GeneticsFace/R/Rdata.R');
source('../GeneticsFace/R/Rsystem.R');
source('../GeneticsFace/R/Rgraphics.R');
source('../GeneticsFace/R/facePreprocess.R');
source('../GeneticsFace/R/faceData.R');
source('../GeneticsFace/R/facePlotting.R');
source('../GeneticsFace/R/faceSymmetry.R');
source('../GeneticsFace/R/faceFeatures.R');
source('../GeneticsFace/R/faceClassification.R');
source('../GeneticsFace/R/faceVisualization.R');

if (0) {
	load('classifier.RData');
	importancePlot(meanGraph, model, modelDesc, pars = list(), output, average);
	#importancePlot(meanGraph=gr$graphs, model=rClass$model, modelDesc=dataFeature$desc, pars = list(), output, average);	
}

if (1) {
	load('classifier-data.RData');
	r = lapply(r, function(e)if (is.character(e)) gsub('data/pictures_converted', '.', e) else e);
	importancePlots(
		r$coords, r$group, model, dataFeature$desc,
		output = r$output, averageInput = r$outputAverages
	);
}
