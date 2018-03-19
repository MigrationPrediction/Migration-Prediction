load('label1090');
a = find(label1090 == 0);
b = find(label1090 == 1);
c = find(label1090 == 2);
asize = size(a);
bsize = size(b);
no = min(asize(2), bsize(2))/2;
asize = asize(2);
bsize = bsize(2);
an = randperm(asize);
bn = randperm(bsize);
seqtest = [a(an(1:no)) b(bn(1:no))];
seqtrain = [c a(an(no+1:end)) b(bn(no+1:end))];

clear all, clc
load('1090gp1_test')
load('1090gp1_train')