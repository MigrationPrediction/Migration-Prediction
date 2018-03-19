% Training a SVM for mitochondria three-class
% classification
clear,clc
[inputs, targets] = trainnew;
[inpu, targ] = testnew;
rng(1); % For reproducibility
SVMModels = cell(3,1);
classes = unique(targets);

for j = 1:numel(classes);
    indx = strcmp(targets,classes(j)); % Create binary classes for each classifier
    SVMModels{j} = fitcsvm(inputs,indx,'ClassNames',[false true],'Standardize',true,...
        'KernelFunction','rbf','BoxConstraint',1);
end
predict_label1 = predict(SVMModels{1}, inpu)';
predict_label2 = predict(SVMModels{2}, inpu)';
predict_label3 = predict(SVMModels{3}, inpu)';
temp=size(targ);

for i=1:temp(2)
    if(predict_label1(i)==1)
        outp(:,i)=[0;1;0];
    elseif(predict_label2(i)==1)
        outp(:,i)=[1;0;0];
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
plotconfusion(tar, outp);
