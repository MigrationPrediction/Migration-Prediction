clear all,clc
for count1 = [37]

c = xlsread('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Nucleus_Label_Check\Data.xlsx', 1, strcat('E', num2str(count1)));
r = xlsread('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Nucleus_Label_Check\Data.xlsx', 1, strcat('F', num2str(count1)));

c1 = xlsread('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Nucleus_Label_Check\Correct.xlsx', 1, strcat('E', num2str(count1)));
r1 = xlsread('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Nucleus_Label_Check\Correct.xlsx', 1, strcat('F', num2str(count1)));

Error(count1) = sqrt((c-c1)*(c-c1)+(r-r1)*(r-r1))/sqrt(2);
end