clear all; close all; clc;
Filepath = '\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells';

CenterMove = load('CenterMove.mat');
Data = load('Data.mat');
FiberRatio = load('Dataprocessimg');
FiberRatio = FiberRatio.data;
CenterMove = CenterMove.data';
Data = Data.data;

threshold = 15;
threshold1 = 150;
threshold2 = 1000;
boundarythres = 10;
widthmax=300;
image_minsize = 10000;
widthred = 50;
heightred = 100;

%[filenameM, pathnameM, filterindexM] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'},'Please select the multi-channel image', Filepath);
for count1=23:40
    for count2=1:2

        L = Data((count1+1952)*2-2+count2, 3);
        R = Data((count1+1952)*2-2+count2, 4);
        T = Data((count1+1952)*2-2+count2, 5);
        B = Data((count1+1952)*2-2+count2, 6);
        
        Center1 = Data((count1+1952)*2-2+count2, 1);
        Center2 = Data((count1+1952)*2-2+count2, 2);
        
        

temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cell Movement 0702\ChamberImg', num2str(count1*2+80*(count2-1), '%04i'));
    ImgMLocation = strcat(temp, '.tif');
    
Image_Original = imread(ImgMLocation);
%Imageperctile = prctile(Image_Original(:), 98);
Image_G = Image_Original(1:end,1:end,1);


threshold = 7;

if count1 == 22
    threshold = 6;
end



smooth = Image_G;
smooth = smooth .* uint8(smooth > threshold);
AreaFilter = bwareaopen(smooth, 10, 4);
Image_G = smooth .* uint8(AreaFilter);
ProcessedImg = Image_G;
bw1 = ProcessedImg;
bw1 = bwareaopen(bw1, 100, 4);

bw2 = bwselect(bw1, Center1, Center2, 4);
X_Sum = sum(bw2);
Y_Sum = sum(bw2');
L_Bound = min(find(X_Sum > 0));
R_Bound = max(find(X_Sum > 0));
T_Bound = min(find(Y_Sum > 0));
B_Bound = max(find(Y_Sum > 0));


Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
Cropped_Image = imcrop(bw1&bw2, Rect);


[L, num] = bwlabel(Cropped_Image,4);
maxarea = 1;
maxno = 1;

    for i = 1:num
        DetectRegion = (L == i);
        temparea = sum(sum(DetectRegion));
        if temparea > maxarea
            maxno = i;
            maxarea = temparea;
        end
    end
    
    tem1 = L_Bound;
    tem2 = T_Bound;
    
    DetectRegion = (L == maxno);
    X_Sum = sum(DetectRegion);
    Y_Sum = sum(DetectRegion');
    L_Bound = min(find(X_Sum > 0))+tem1;
    R_Bound = max(find(X_Sum > 0))+tem1;
    T_Bound = min(find(Y_Sum > 0))+tem2;
    B_Bound = max(find(Y_Sum > 0))+tem2;
    Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
    Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
    Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
    Cropped_Image = imcrop(bw1&bw2, Rect);

Weighted_Center = CenterMass(uint8(bw1));
%figure; imshow(Image_Original);
%hold on; plot(Weighted_Center(2), Weighted_Center(1), 'r*');
%hold on; plot(Graph_Center(1), Graph_Center(2),'bd');

Weighted_Center2 = CenterMass(uint8(Cropped_Image));
BW_path = strcat(Filepath,'\Process_Img\BW', num2str(count1+1952), '_', num2str(count2), '.jpg');
%figure; imshow(Cropped_Image);
Cropped_Image = imcrop(ProcessedImg,Rect);
Cropped_Image = Cropped_Image.*uint8((Cropped_Image>6));
if count2 == 1
    Weighted_Center3 = CenterMass(uint8(Cropped_Image));
    L1 = Weighted_Center3(2)+L_Bound;
    T1 = Weighted_Center3(1)+T_Bound;
else
    Weighted_Center4 = CenterMass(uint8(Cropped_Image));
    L2 = Weighted_Center4(2)+L_Bound;
    T2 = Weighted_Center4(1)+T_Bound;
    xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\Center.xlsx', T2-T1, 1, strcat('A', num2str(count1)));
    xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\Center.xlsx', L2-L1, 1, strcat('B', num2str(count1)));
end

%imwrite(Cropped_Image, BW_path);
%hold on; plot(Weighted_Center2(2), Weighted_Center2(1), 'r*');
%hold on; plot(0.5*(R_Bound - L_Bound), 0.5*(B_Bound - T_Bound), 'bd');


Cro_Orig = imcrop(Image_Original, Rect);
Cro_bw = imcrop(bw1&bw2, Rect);
FL_path = strcat(Filepath,'\Process_FL\FL', num2str(count1+1952), '_', num2str(count2), '.jpg');
%figure; imshow(Cro_Orig);
%imwrite(Cro_Orig, FL_path);
%hold on; plot(Weighted_Center2(2), Weighted_Center2(1), 'r*');
%hold on; plot(0.5*(R_Bound - L_Bound), 0.5*(B_Bound - T_Bound), 'bd');


%figure;
%imshow(bw1);
bw2 = bw1;
%{
temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cell Movement 0702\ChamberImg', num2str(count1*2-1+80*(count2-1), '%04i'));
    ImgMLocation = strcat(temp, '.tif');
    graythre = 0.1;
    
Image_Original = imread(ImgMLocation);
Imageperctile = prctile(Image_Original(:), 99);
Image_G = Image_Original(1:end,1:end,2);

[counts,binLocations] = imhist(Image_G);
diffcounts = counts(1:end-1)-counts(2:end);
threshold = find(diffcounts == max(diffcounts(4:end)));

if (Imageperctile > 10 && threshold < 20)
    threshold = Imageperctile;
end

if Imageperctile < 4
    threshold = 4;
end

if threshold >15
    threshold = 15;
end

smooth = Image_G;
smooth = smooth .* uint8(smooth > threshold);
AreaFilter = bwareaopen(smooth, 10, 4);
Image_G = smooth .* uint8(AreaFilter);

Background = imopen(Image_G,strel('disk',3));
ProcessedImg = imadjust(Image_G-Background);
K = wiener2(ProcessedImg,[5 5]);

Centroid = regionprops(K,'Centroid');

R_path = strcat(Filepath,'\Process_R\R', num2str(count1+1952), '_', num2str(count2), '.jpg');
imwrite(Cro_Orig, R_path);
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\Data.xlsx', L_Bound, 1, strcat('A', num2str((count1+1253-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\Data.xlsx', R_Bound, 1, strcat('B', num2str((count1+1253-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\Data.xlsx', T_Bound, 1, strcat('C', num2str((count1+1253-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\Data.xlsx', B_Bound, 1, strcat('D', num2str((count1+1253-1)*2+count2)));

level1 = graythresh(K)*graythre;
bw1 = im2bw(ProcessedImg, level1);

 %   Cro_Orig(:,:,1) = 100.*uint8(Cro_Orig(:,:,1)>=8);

Cro_Orig(:,:,1) = 100.*Cro_bw;

bw1 = bwareaopen(bw1, 15, 4);
Cropped_Image2 = imcrop(bw1, Rect);
Cro_Orig(:,:,2) = Cropped_Image2*100;

figure;
imshow(Cro_Orig);
RG_path = strcat(Filepath,'\Process_RG\RG', num2str(count1+1952), '_', num2str(count2), '.jpg');
imwrite(Cro_Orig, RG_path);
%}
    end
end