% This is the program using random forest method to predict
% cell moving direction 
clear all, clc

% load feature matrix and labels for cell moving direction
fastslow = load('labeldirection.mat');
NormalizedFeatureMatrix = load('CellFeatureMatrixNewNonAbs.mat');
fastslow = fastslow.labeldirection;
NormalizedFeatureMatrix = NormalizedFeatureMatrix.CellFeatureMatrix;

% randomize the train, validation, test sequences
non = 1:61;
seq = randperm(1358);

T = int16(fastslow);
% change labels to one-hot form
for i = 1:1358
    if T(i) == 0
        nex(:,i) = [1;0];
    else
        nex(:,i) = [0;1];
    end
end
T = nex;

% choose number of internal nodes for each layer
node1 = 21;
node2 = 7;

prcerror1 = 0;
x = NormalizedFeatureMatrix;
t = T;
pe = 1;
percentErrors = 1;


trainFcn = 'trainscg';

% Create a Pattern Recognition Network
hiddenLayerSize = [node1,node2];
net = patternnet(hiddenLayerSize);

net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};
%net.Layers{1}.transferFcn = 'logsig';
%net.Layers{2}.transferFcn = 'logsig';

% Change Network Parameter Setting
net.trainParam.max_fail = 100;
net.trainParam.epochs = 1000;
net.trainParam.min_grad = 0.0000001;


    net.divideFcn = 'divideind';
    net.divideParam.trainInd = seq(1:950);
    net.divideParam.valInd = seq(951:1154);
    net.divideParam.testInd = seq(1155:1358);


net.performFcn = 'mse';

net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
    'plotconfusion', 'plotroc'};

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y)
tind = vec2ind(t);
yind = vec2ind(y);
tind1 = tind(net.divideParam.testInd);
yind1 = yind(net.divideParam.testInd);
percentErrors = sum(tind1 ~= yind1)/numel(tind1);
pe = sum(tind ~= yind)/numel(tind);

% Recalculate Training, Validation and Test Performance
trainTargets = t .* tr.trainMask{1};
valTargets = t .* tr.valMask{1};
testTargets = t .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,y)
valPerformance = perform(net,valTargets,y)
testPerformance = perform(net,testTargets,y)


if (false)
    
    genFunction(net,'myNeuralNetworkFunction');
    y = myNeuralNetworkFunction(x);
end
if (false)
    genFunction(net,'myNeuralNetworkFunction','MatrixOnly','yes');
    y = myNeuralNetworkFunction(x);
end
if (false)
    gensim(net);
end
