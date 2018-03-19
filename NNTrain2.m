% This program has two separate parts. 
% The first one deals with 20/80 motile/nonmotile cells
% The second one deals with 10/90 motile/nonmotile cells
clear all, clc
% load whiten cell feature matrix and their labels
fastslow = load('fastslowbinary.mat');
NormalizedFeatureMatrix = load('WhitenData.mat');
fastslow = fastslow.fastslow;
NormalizedFeatureMatrix = NormalizedFeatureMatrix.a;
NormalizedFeatureMatrix = NormalizedFeatureMatrix';
non = 1:61;
load('label2080');
load('gp4_1')
load('gp4_2')

T = int16(fastslow);

% change labels to one-hot form
for i = 1:1358
    if T(i) == 1
        nex(:,i) = [1;0];
    else
        nex(:,i) = [0;1];
    end
end
T = nex;

% Initialize prediction accuracy matrix for feature reduction
prcerror = zeros(61,60);

% Set number of nodes for each layer
node1 = 21;
node2 = 7;

prcerror1 = 0;
percentErrors = 1;
x = NormalizedFeatureMatrix;
t = T;
minerror = 1;

% Run the original feature matrix (without feature reduction)
% Record the average and the best result
for count1 = 1:400  
    
    
trainFcn = 'trainscg';

hiddenLayerSize = [node1,node2];
net = patternnet(hiddenLayerSize);

net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};


net.trainParam.max_fail = 20;
net.trainParam.epochs = 1000;
net.trainParam.min_grad = 1e-8;

    net.divideFcn = 'divideind';
    net.divideParam.trainInd = [seqtrain(1:950), seqtest(81:170)];
    net.divideParam.valInd = seqtrain(951:1108);
    net.divideParam.testInd = [seqtest(1:80), seqtest(171:250)];
    

net.performFcn = 'mse';

net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
    'plotconfusion', 'plotroc'};

[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y)
tind = vec2ind(t);
yind = vec2ind(y);
y1 = y(net.divideParam.testInd);
tind1 = tind(net.divideParam.testInd);
yind1 = yind(net.divideParam.testInd);
[XN,YN,TN, AUC] = perfcurve(tind1,yind1,2);
percentErrors = sum(tind1 ~= yind1)/numel(tind1);

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
prcerror1 = prcerror1+percentErrors/400;
if percentErrors < minerror
    minerror = percentErrors;
end
end


xlswrite('Feature2080_1s.xlsx', prcerror1,4,strcat('A',num2str(1)));
xlswrite('Feature2080_1s.xlsx', minerror,6,strcat('A',num2str(1)));

% Begin feature reduction
% Delete one feature at a time
% Record the prediction accuracy matrix in XLSX file for each deletion
for im = 1:60
matrixsize = size(NormalizedFeatureMatrix);
minerror = 1;
for ik = 1:matrixsize(1)
x = NormalizedFeatureMatrix;
x(ik,:) = zeros(1,1358);
t = T;
for count1 = 1:20

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
% trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.
trainFcn = 'trainscg';

% Create a Pattern Recognition Network
hiddenLayerSize = [node1,node2];
net = patternnet(hiddenLayerSize);

% Choose Input and Output Pre/Post-Processing Functions
net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};

% Change Network Parameter Setting
net.trainParam.max_fail = 20;
net.trainParam.epochs = 1000;
net.trainParam.min_grad = 1e-8;

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide


    net.divideFcn = 'divideind';
    net.divideParam.trainInd = [seqtrain(321:950), seqtest(1:170)];
    net.divideParam.valInd = seqtrain(951:1108);
    net.divideParam.testInd = [seqtrain(1:320), seqtest(171:250)];



% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
% net.performFcn = 'crossentropy';  % Cross-Entropy
net.performFcn = 'mse';

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
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

% Recalculate Training, Validation and Test Performance
trainTargets = t .* tr.trainMask{1};
valTargets = t .* tr.valMask{1};
testTargets = t .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,y)
valPerformance = perform(net,valTargets,y)
testPerformance = perform(net,testTargets,y)

% View the Network
% view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotconfusion(t,y)
%figure, plotroc(t,y)

% Deployment
% Change the (false) values to (true) to enable the following code blocks.
% See the help for each generation function for more information.
if (false)
    % Generate MATLAB function for neural network for application
    % deployment in MATLAB scripts or with MATLAB Compiler and Builder
    % tools, or simply to examine the calculations your trained neural
    % network performs.
    genFunction(net,'myNeuralNetworkFunction');
    y = myNeuralNetworkFunction(x);
end
if (false)
    % Generate a matrix-only MATLAB function for neural network code
    % generation with MATLAB Coder tools.
    genFunction(net,'myNeuralNetworkFunction','MatrixOnly','yes');
    y = myNeuralNetworkFunction(x);
end
if (false)
    % Generate a Simulink diagram for simulation or deployment with.
    % Simulink Coder tools.
    gensim(net);
end
prcerror(ik,im) = prcerror(ik,im)+percentErrors/20;
if percentErrors < minerror
    minerror = percentErrors;
end
end
end
maxerror = min(prcerror(1:(62-im),im));
delete = find(prcerror(:,im) == maxerror);
delete = delete(1);
deleteno = non(delete);
NormalizedFeatureMatrix(delete,:)=[];
non(delete) = [];
xlswrite('Feature2080_1s.xlsx', prcerror(:,im)',strcat('A',num2str(ik),':CK',num2str(ik)));
xlswrite('Feature2080_1s.xlsx', deleteno, 2, strcat('A',num2str(ik)));
xlswrite('Feature2080_1s.xlsx', minerror, 5, strcat('A',num2str(ik)));
xlswrite('Feature2080_1s.xlsx', maxerror, 3, strcat('A',num2str(ik)));
end


% load whiten cell feature matrix and their labels 
clear all, clc
fastslow = load('fastslowbinary.mat');
NormalizedFeatureMatrix = load('WhitenData.mat');
fastslow = fastslow.fastslow;
NormalizedFeatureMatrix = NormalizedFeatureMatrix.a;
NormalizedFeatureMatrix = NormalizedFeatureMatrix';
non = 1:61;

load('1090gp2_test')
load('1090gp2_train')
load('label1090')


T = int16(fastslow);
% change labels to one-hot form
for i = 1:1358
    if label1090(i) == 1
        nex(:,i) = [1;0];
    else
        nex(:,i) = [0;1];
    end
end
T = nex;

% Initialize prediction accuracy matrix for feature reduction
prcerror = zeros(61,60);

% Set number of nodes for each layer
node1 = 21;
node2 = 7;

prcerror1 = 0;
percentErrors = 1;
x = NormalizedFeatureMatrix;
t = T;
minerror = 1;

% Run the original feature matrix (without feature reduction)
% Record the average and the best result
for count1 = 1:400
    
trainFcn = 'trainscg';

hiddenLayerSize = [node1,node2];
net = patternnet(hiddenLayerSize);

net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};


net.trainParam.max_fail = 20;
net.trainParam.epochs = 1000;
net.trainParam.min_grad = 1e-8;


    net.divideFcn = 'divideind';
    net.divideParam.trainInd = [seqtrain(361:1050), seqtest(1:94)];
    net.divideParam.valInd = seqtrain(1051:1224);
    net.divideParam.testInd = [seqtrain(1:360), seqtest(95:134)];


net.performFcn = 'mse';

net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
    'plotconfusion', 'plotroc'};

[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y)
tind = vec2ind(t);
yind = vec2ind(y);
y1 = y(net.divideParam.testInd);
tind1 = tind(net.divideParam.testInd);
yind1 = yind(net.divideParam.testInd);
%y2 = [-ones(1,125),ones(1,125)] ;
[XN,YN,TN, AUC] = perfcurve(tind1,yind1,2);
percentErrors = sum(tind1 ~= yind1)/numel(tind1);

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
prcerror1 = prcerror1+percentErrors/400;
if percentErrors < minerror
    minerror = percentErrors;
end
end


xlswrite('Feature1090_2s.xlsx', prcerror1,4,strcat('A',num2str(1)));
xlswrite('Feature1090_2s.xlsx', minerror,6,strcat('A',num2str(1)));


% Begin feature reduction
% Delete one feature at a time
% Record the prediction accuracy matrix in XLSX file for each deletion
for im = 1:60
matrixsize = size(NormalizedFeatureMatrix);
minerror = 1;
for ik = 1:matrixsize(1)
x = NormalizedFeatureMatrix;
x(ik,:) = zeros(1,1358);
t = T;
for count1 = 1:20

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
% trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.
trainFcn = 'trainscg';

% Create a Pattern Recognition Network
hiddenLayerSize = [node1,node2];
net = patternnet(hiddenLayerSize);

% Choose Input and Output Pre/Post-Processing Functions
net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};

% Change Network Parameter Setting
net.trainParam.max_fail = 20;
net.trainParam.epochs = 1000;
net.trainParam.min_grad = 1e-8;

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide

    net.divideFcn = 'divideind';
    net.divideParam.trainInd = [seqtrain(361:1050), seqtest(1:94)];
    net.divideParam.valInd = seqtrain(1051:1224);
    net.divideParam.testInd = [seqtrain(1:360), seqtest(95:134)];

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
% net.performFcn = 'crossentropy';  % Cross-Entropy
net.performFcn = 'mse';

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
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

% Recalculate Training, Validation and Test Performance
trainTargets = t .* tr.trainMask{1};
valTargets = t .* tr.valMask{1};
testTargets = t .* tr.testMask{1};
trainPerformance = perform(net,trainTargets,y)
valPerformance = perform(net,valTargets,y)
testPerformance = perform(net,testTargets,y)

% View the Network
% view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotconfusion(t,y)
%figure, plotroc(t,y)

% Deployment
% Change the (false) values to (true) to enable the following code blocks.
% See the help for each generation function for more information.
if (false)
    % Generate MATLAB function for neural network for application
    % deployment in MATLAB scripts or with MATLAB Compiler and Builder
    % tools, or simply to examine the calculations your trained neural
    % network performs.
    genFunction(net,'myNeuralNetworkFunction');
    y = myNeuralNetworkFunction(x);
end
if (false)
    % Generate a matrix-only MATLAB function for neural network code
    % generation with MATLAB Coder tools.
    genFunction(net,'myNeuralNetworkFunction','MatrixOnly','yes');
    y = myNeuralNetworkFunction(x);
end
if (false)
    % Generate a Simulink diagram for simulation or deployment with.
    % Simulink Coder tools.
    gensim(net);
end
prcerror(ik,im) = prcerror(ik,im)+percentErrors/20;
if percentErrors < minerror
    minerror = percentErrors;
end
end
end
maxerror = min(prcerror(1:(62-im),im));
delete = find(prcerror(:,im) == maxerror);
delete = delete(1);
deleteno = non(delete);
NormalizedFeatureMatrix(delete,:)=[];
non(delete) = [];
xlswrite('Feature1090_2s.xlsx', prcerror(:,im)',strcat('A',num2str(ik),':CK',num2str(ik)));
xlswrite('Feature1090_2s.xlsx', deleteno, 2, strcat('A',num2str(ik)));
xlswrite('Feature1090_2s.xlsx', minerror, 5, strcat('A',num2str(ik)));
xlswrite('Feature1090_2s.xlsx', maxerror, 3, strcat('A',num2str(ik)));
end