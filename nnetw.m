clear all, clc
% Training a Deep Neural Network for mitochondria three-class
% classification

% Load the training data into memory

% XTrainImages are features, tTrain is labels
[xTrainImages, tTrain] = traindata;

rng('default')

% Set the size of the hidden layer for the autoencoder. 
hiddenSize1 = 10;

autoenc1 = trainAutoencoder(xTrainImages,hiddenSize1, ...
    'MaxEpochs',400, ...
    'L2WeightRegularization',0.004, ...
    'SparsityRegularization',4, ...
    'SparsityProportion',0.15, ...
    'ScaleData', false);


view(autoenc1)

%plotWeights(autoenc1);

feat1 = encode(autoenc1,xTrainImages);

% Training the second autoencoder

hiddenSize2 = 5;
autoenc2 = trainAutoencoder(feat1,hiddenSize2, ...
    'MaxEpochs',100, ...
    'L2WeightRegularization',0.002, ...
    'SparsityRegularization',4, ...
    'SparsityProportion',0.1, ...
    'ScaleData', false);

view(autoenc2)

feat2 = encode(autoenc2,feat1);


% Training the final softmax layer

softnet = trainSoftmaxLayer(feat2,tTrain,'MaxEpochs',400);

view(softnet)

% Forming a stacked neural network

view(autoenc1)
view(autoenc2)
view(softnet)

deepnet = stack(autoenc1,autoenc2,softnet);

view(deepnet)

% Actually just 16 features for each image
imageWidth = 1;
imageHeight = 16;
inputSize = imageWidth*imageHeight;

% Load the test data
[xTestImages, tTest] = testdata;
xTest = zeros(inputSize,numel(xTestImages));
for i = 1:numel(xTestImages)
    xTest(:,i) = xTestImages{i}(:);
end

% You can visualize the results with a confusion matrix. 

y = deepnet(xTest);
plotconfusion(tTest,y);

xTrain = zeros(inputSize,numel(xTrainImages));
for i = 1:numel(xTrainImages)
    xTrain(:,i) = xTrainImages{i}(:);
end

% Perform fine tuning
deepnet = train(deepnet,xTrain,tTrain);


y = deepnet(xTest);
plotconfusion(tTest,y);


displayEndOfDemoMessage(mfilename)