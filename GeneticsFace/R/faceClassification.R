#
#	faceClassification.R
#Tue Sep 29 14:51:50 2015

require('glmnet');

glmnet_test_binomial = function(model, data, threshold = .5) {
	r = as.vector(as.matrix(data) %*% model);
	pred = as.integer(plogis(r) > threshold);
	pred
}

glmnet_train = function(data, nfolds = -1, family = 'binomial', type.measure = 'class', alpha = .5) {
	# <p> crossvalidation partitions
	N = nrow(data);
	perm = Sample(1:N, N);
	fold = if (nfolds > 0) nfolds else as.integer(N/-nfolds);
	folds = index2listPosition(splitListEls(perm, fold, returnElements = T));
	# <p> glmnet
	X = as.matrix(Df(data, min_ = 'response'));
	Y = data$response;
	r = cv.glmnet(X, Y, foldid = folds, type.measure = type.measure,
		grouped = FALSE, alpha = alpha, family = family);
	r0 = glmnet(X, Y, lambda = r$lambda.min, alpha = alpha, family = family);
	r1 = coefficients(r0);
	r1
}

#@arg cv_fold By default leave-one-out crossvalidation is performed: cv_fold == -1

classifyFaceFeatures = function(response, feature, cv_fold = -1, parallel = FALSE,
	trainer = glmnet_train, tester = glmnet_test_binomial, Nrepeat = 1) {
	feature = t(apply(feature, 1, standardize));
	data = data.frame(response = response, feature = feature);
	crossvalidate(trainer, tester, data = data, parallel = parallel, cv_repeats = Nrepeat);
}

