% This is the program using random forest method to predict
% cell moving direction and sort important features
clear all,clc
% load feature matrix and labels for cell moving direction
fastslow = load('labeldirection.mat');
NormalizedFeatureMatrix = load('CellFeatureMatrixNewNonAbs.mat');
train1 = load('celltrain');
test1 = load('celltest');
fastslow = fastslow.labeldirection;
NormalizedFeatureMatrix = NormalizedFeatureMatrix.CellFeatureMatrix;
train1 = train1.train;
test1 = test1.test;

FeatureName = load('FeatureName.mat');
FeatureName = FeatureName.name;

seq = randperm(1358);
X = NormalizedFeatureMatrix(:,seq);
T = int16(fastslow(seq));


for i = 1:1358
    if T(i) == 0
        nex{i} = 'up';
    else
        nex{i} = 'down';
    end
end

T = nex;

inputs = X';
targets = T;

rng(1); % For reproducibility
Mdl = TreeBagger(1000,inputs,targets,'OOBPrediction','On','OOBPredictorImportance', 'On','Surrogate','on',...
    'Method','classification');

oobErrorBaggedEnsemble = oobError(Mdl);
figure; plot(oobErrorBaggedEnsemble)
xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');


disp('Sorting importance into descending order')
weights=Mdl.OOBPermutedVarDeltaError;
figure;bar(weights);
title('Feature Importance for Cell Movement');
ylabel('Predictor importance estimates');
xlabel('Predictors');
h = gca;
h.FontSize = 12;
set(gca,'XTick',1:61);
h.XTickLabel = FeatureName;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
hold on 

%{
bar(weights.*(weights>0),'c');
bar(weights.*(weights>0.08),'g');
bar(weights.*(weights>0.13),'y');
bar(weights.*(weights>0.25),'r');
%}


bar(weights.*(weights>0),'c');
bar(weights.*(weights>0.2),'g');
bar(weights.*(weights>0.4),'y');
bar(weights.*(weights>0.8),'r');



[B,iranked] = sort(weights,'descend');
disp(['Plotting a horizontal bar graph of sorted labeled weights.']) 
figure
barh(weights(iranked),'g');
xlabel('Variable Importance','FontSize',30,'Interpreter','latex');
ylabel('Variable Rank','FontSize',30,'Interpreter','latex');
title(...
    ['Relative Importance of Inputs'],...
    'FontSize',17,'Interpreter','latex'...
    );
hold on
barh(weights(iranked(1:10)),'y');
barh(weights(iranked(1:5)),'r');



%{
y = flip(y);
x = flip(x);
figure;barh(x,'y');
hold on
%barh(x.*(x>0.5),'g');
%hold on
barh(x.*(x>0.7),'r');
h = gca;
set(h,'fontWeight','bold','fontsize',12);
h.YTick = 1:numel(y);
h.YTickLabel = y(1:end);
title('Top 10 Important Features for Direction Prediction');
xlabel('Predictor importance estimates');
ylabel('Predictors');
hold on 
%}
