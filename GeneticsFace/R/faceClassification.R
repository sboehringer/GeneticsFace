#
#	faceClassification.R
#Tue Sep 29 14:51:50 2015

require('glmnet');

crossvalidationPartitions = function(N, Nfolds = -1) {
	perm = Sample(1:N, N);
	fold = if (Nfolds > 0) Nfolds else as.integer(N/-Nfolds);
	folds = index2listPosition(splitListEls(perm, fold, returnElements = T));
	folds
}

glmnet_lambda_min = function(X, Y,..., Nfolds = -1, Nrepeats = if (Nfolds == -1) 1 else 10) {
	sapply(1:Nrepeats, function(i) {
		folds = crossvalidationPartitions(nrow(X), Nfolds);
		cv.glmnet(X, Y, ..., foldid = folds)$lambda.min
	})
}

glmnet_test_binomial = function(model, data, threshold = .5) {
	r = as.vector(as.matrix(data) %*% model);
	pred = as.integer(plogis(r) > threshold);
	pred
}

glmnet_train = function(data, cv_fold = -1, family = 'binomial', type.measure = 'class', alpha = .5) {
	# <p> glmnet
	X = as.matrix(Df(data, min_ = 'response'));
	Y = data$response;
	r = glmnet_lambda_min(X, Y, Nfolds = cv_fold,
		type.measure = type.measure, grouped = FALSE, alpha = alpha, family = family);
	r0 = glmnet(X, Y, lambda = median(r), alpha = alpha, family = family);
	r1 = coefficients(r0);
	r1
}

#@arg cv_fold By default leave-one-out crossvalidation is performed: cv_fold == -1

classifyFaceFeatures = function(response, feature, cv_fold = -1, parallel = FALSE,
	trainer = glmnet_train, tester = glmnet_test_binomial, Nrepeat = 1, cv_fold_inner = -1,
	alpha = .5, family = 'binomial') {
	feature = t(apply(feature, 1, standardize));
	data = data.frame(response = response, feature = feature);
	modelData = trainer(data, cv_fold = cv_fold_inner,
		family = family, type.measure = 'class', alpha = alpha);
	pred = crossvalidate(trainer, tester, data = data, parallel = parallel, cv_repeats = Nrepeat);
	r = list(model = modelData, prediction = pred);
	r
}

