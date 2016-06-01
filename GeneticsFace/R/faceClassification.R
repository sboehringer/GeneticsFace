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

glmnet_test_multinomial = function(model, data, threshold = .5) {
	X = as.matrix(Df(data, min_ = 'response'));
	pred = predict(model, X, type = 'response');
	predClass = apply(pred, 1, which.max);
	predClass
}

glmnet_train = function(data, cv_fold = -1, family = 'multinomial', type.measure = 'class', alpha = .5) {
	# <p> glmnet
	X = as.matrix(Df(data, min_ = 'response'));
	Y = data$response;
	r = glmnet_lambda_min(X, Y, Nfolds = cv_fold,
		type.measure = type.measure, grouped = FALSE, alpha = alpha, family = family);
	r0 = glmnet(X, Y, lambda = median(r), alpha = alpha, family = family);
	r0
}

# create matrix with coefficients
glmnetModel = function(r0) {
	r1 = do.call(cbind, lapply(as.list(coefficients(r0)), as.vector));
	dimnames(r1) = list(row.names(coefficients(r0)[[1]]), names(coefficients(r0)));
	r1
}

aPrioriClass = function(group) {
	freqs = table(group)/length(group);
	r = (freqs %*% freqs)[1, 1];
	r
}

testClassification = function(group, prediction) {
	truthPred = cbind(as.character(group), unlist(prediction));
	accuracy = mean(matrix.same(truthPred));
	aP = aPrioriClass(group);
	N = length(group);
	bt = binom.test(round(accuracy * N, 0), N, p = aP, alternative = 'greater');
	r = list(N = N, accuracyApriori = aP, accuracy = accuracy, p.value = bt$p.value, groups = table(group));
	r
}

#@arg cv_fold By default leave-one-out crossvalidation is performed: cv_fold == -1

classifyFaceFeaturesRaw = function(response, feature, cv_fold = -1, parallel = FALSE,
	trainer = glmnet_train, tester = glmnet_test_multinomial, Nrepeat = 1, cv_fold_inner = -1,
	alpha = .5, family = 'multinomial') {
	feature = t(apply(feature, 1, standardize));
	data = data.frame(response = response, feature = feature);
	respClasses = if (is.factor(response)) levels(response) else unique(response);
	modelData = trainer(data,
		family = family, type.measure = 'class', alpha = alpha);
	pred = crossvalidate(trainer, tester, data = data, parallel = parallel,
		cv_fold = cv_fold_inner, cv_repeats = Nrepeat);
	# <p> reinstate original response-labels
	pred1 = lapply(pred, function(e)sapply(e, function(e) {
		r = respClasses[e];
		names(r) = names(e);
		r
	}));
	r = list(model = modelData, alpha = alpha, coefficients = glmnetModel(modelData),
		prediction = pred1,
		test = testClassification(response, pred1)
	);
	r
}

classifyFaceFeatures = function(response, feature, cv_fold = -1, parallel = FALSE,
	trainer = glmnet_train, tester = glmnet_test_multinomial, Nrepeat = 1, cv_fold_inner = -1,
	alphas = seq(.05, .5, length.out = 10), family = 'multinomial') {
	rAlphas = lapply(alphas, function(alpha)
		classifyFaceFeaturesRaw(response, feature, cv_fold, parallel,
			trainer, tester, Nrepeat, cv_fold_inner, alpha, family));
	iBest = which.min(list.kp(rAlphas, 'test$accuracy', do.unlist = T));
	r = c(rAlphas[[iBest]], list(classificationAlphas = rAlphas));
	r
}

