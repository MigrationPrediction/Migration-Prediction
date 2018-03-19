% This program runs mitochondria classification to three categories:
% Filamentous, intermediate, and dots
% This program uses random forest method to do classification
clear,clc
[inputs, targets] = train3;
[inpu, targ] = test3;
rng(1); % For reproducibility
Mdl = TreeBagger(50,inputs,targets,'OOBPrediction','On','OOBPredictorImportance', 'On','Surrogate','on',...
    'Method','classification');
%{
for i = 1 : 2052
    if(strcmp(targets(i),'fiber'))
        pre(i)=1;
    elseif(strcmp(targets(i),'dot'))
        pre(i)=2;
    else
        pre(i)=3;
    end
end
Mdl = fitensemble(inputs,pre,'LSBoost',100,'Tree');            
imp = predictorImportance(Mdl);
%}

predict_label = predict(Mdl, inpu)';
temp=size(targ);
% Change labels to one-hot form
for i=1:temp(2)
    if(strcmp(predict_label{i},'fiber'))
        outp(:,i)=[1;0;0];
    elseif(strcmp(predict_label{i},'dot'))
        outp(:,i)=[0;1;0];
    else
        outp(:,i)=[0;0;1];
    end
    if(strcmp(targ{i},'fiber'))
        tar(:,i)=[1;0;0];
    elseif(strcmp(targ{i},'dot'))
        tar(:,i)=[0;1;0];
    else
        tar(:,i)=[0;0;1];
    end
end
wrong=[];
for i=1:temp(2)
    if(((outp(:,i)')*tar(:,i))==0)
        if(outp(:,i)==[1;0;0])
            o=1;
        elseif(outp(:,i)==[0;1;0])
            o=2;
        else
            o=3;
        end
        if(tar(:,i)==[1;0;0])
            t=1;
        elseif(tar(:,i)==[0;1;0])
            t=2;
        else
            t=3;
        end
        wrong=[wrong,[i;o;t]];
    end
end
figure;plotconfusion(tar, outp);

oobErrorBaggedEnsemble = oobError(Mdl);
figure; plot(oobErrorBaggedEnsemble)
xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');

%{
disp('Sorting importance into descending order')
weights=Mdl.OOBPermutedVarDeltaError;
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
%}
