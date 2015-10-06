#
#	expample illustrating feature extraction from classification model
#Tue Oct  6 14:16:24 CEST 2015

library('tools');
library('shapes');
library('EBImage');
library('geometry');
library('sets');
source('../GeneticsFace/R/Rdata.R');
source('../GeneticsFace/R/facePreprocess.R');
source('../GeneticsFace/R/faceData.R');
source('../GeneticsFace/R/facePlotting.R');
source('../GeneticsFace/R/faceSymmetry.R');
source('../GeneticsFace/R/faceFeatures.R');
source('../GeneticsFace/R/faceClassification.R');

visualizeClassfifier = function(meanGraph, model, modelDesc) {

	cfsL = lapply(c('coordinate', 'distance', 'area', 'angle'), function(feature) {
		print(feature);
		# extract coefficients corresponding to feature from model
		# cfs is a list with components
		#	coefficients: these are regression coefficients from the LASSO model
		#	structure: this describes the corresponding feature further
		#		feature == 'coordinate': matrix with one column, each entry inicates a node corresponding to
		#			*two* coefficients (x, y)
		#		feature == 'distance': matrix with two columns, each row inicates the pair of nodes 
		#			corresponding to the distance
		#		feature == 'area': matrix with three columns, each row inicates the triple of nodes 
		#			corresponding to the triangle with its area
		#		feature == 'angle': matrix with three columns, each row inicates the triple of nodes 
		#			corresponding to the triangle for which *three* coefficients corresond to the angles
		#			in this triangle, the exact layout is described in an ASCII-art in faceFeature.R
		# type == 'feature' means that actual values of the features are extracted
		# type == 'asymm' means that corresponding assymmetry values are extracted
		# for the moment, type == 'feature' can be assumed
		cfs = extractFeatureCoefficients(model, feature, type = 'feature', modelDesc);

		# do something to visualize here
	});
}

if (1) {
	load('classifier.RData');
	visualizeClassfifier(gr$graphs, rClass$model, dataFeature$desc);
}
