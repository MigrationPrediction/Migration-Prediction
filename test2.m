function [inputs,targets] = test2
% digittrain_dataset   Synthetic handwritten digit dataset for training
%   [inputs, targets] = digittrain_dataset will load a dataset of synthetic
%   handwritten digits, where inputs will be a 1-by-5000 cell array, with
%   each cell containing a 28-by-28 matrix representing a synthetic image
%   of a handwritten digit, and targets will be a 10-by-5000 matrix
%   containing the labels for the images in 1-of-K format.
%
%   Example:
%       Train a sparse autoencoder on images of handwritten digits and
%       visualize the learned results.
%
%       x = digittrain_dataset;
%       hiddenSize = 25;
%       autoenc = trainAutoencoder(x, hiddenSize, ...
%           'L2WeightRegularization', 0.004, ...
%           'SparsityRegularization', 4, ...
%           'SparsityProportion', 0.2);
%       plotWeights(autoenc);
%
%   See also digittest_dataset, digitsmall_dataset

% Copyright 2014-2015 The MathWorks, Inc.
B{1}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure12_T1\feature');
B{2}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure12_T2\feature');
B{3}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure13_T1\feature');
B{4}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure14_T1\feature');
B{5}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure14_T2\feature');
B{6}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure15_T1\feature');
B{7}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure16_T1\feature');
B{8}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure17_T1\feature');

A{1}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure12_T1\result');
A{2}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure12_T2\result');
A{3}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure13_T1\result');
A{4}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure14_T1\result');
A{5}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure14_T2\result');
A{6}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure15_T1\result');
A{7}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure16_T1\result');
A{8}=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure17_T1\result');


l = 1;
for i = 1:8;
    F = size(A{i}.B);
    for j = 1:F(2)
        K(:,l) = B{i}.FEA(:,j);
        if(A{i}.D(1,j)==1)
            R{l}='fiber';
        elseif(A{i}.D(2,j)==1)
            R{l}='dot';
        else
            R{l}='inter';
        end
        l = l+1;
    end
end



inputs = K';
targets = R;
end