clear all; close all; clc;
Filepath = '\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells';

threshold = 15;
threshold1 = 150;
threshold2 = 1000;
boundarythres = 10;
graythre = 0.7;
widthmax=400;
image_minsize = 800;

%[filenameM, pathnameM, filterindexM] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'},'Please select the multi-channel image', Filepath);
for count1=8:8
    for count2=1:1
        
%ImgMLocation = strcat(pathnameM, filenameM);
temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Nucleus_Label_Check\ChamberImg', num2str(count1*2-2, '%04i'));
    ImgMLocation = strcat(temp, '.tif');
    graythre = 0.1;
    
Image_Original = imread(ImgMLocation);
Imageperctile = prctile(Image_Original(:), 99);
%Image_G = Image_Original(1:end,1:end,2);
Image_G = Image_Original;

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

if threshold == 8
    threshold = 7;
end

if count1 == 3
    threshold = 7;
end

if count1 == 4
    threshold = 11;
end



smooth = Image_G;
smooth = smooth .* uint8(smooth > threshold);
AreaFilter = bwareaopen(smooth, 10, 4);
Image_G = smooth .* uint8(AreaFilter);

Background = imopen(Image_G,strel('disk',3));
ProcessedImg = imadjust(Image_G-Background);
K = wiener2(ProcessedImg,[5 5]);

Centroid = regionprops(K,'Centroid');

level1 = graythresh(K)*graythre;
bw1 = im2bw(ProcessedImg, level1);
bw1 = bwareaopen(bw1, 15, 4);

X_Sum = sum(bw1);
Y_Sum = sum(bw1');
L_Bound = min(find(X_Sum > 0));


temp=size(bw1);
if(L_Bound+widthmax<temp(2))
    bw11=bw1(1:end, 1:L_Bound+widthmax);
    XX_Sum=sum(bw11);
end


if(L_Bound<15)
    L_Bound=15;
end


R_Bound = max(find(XX_Sum > 0));

T_Bound = min(find(Y_Sum > 0));


B_Bound = max(find(Y_Sum > 0));


Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
Cropped_Image = imcrop(bw1, Rect);

if sum(sum(Cropped_Image))<image_minsize
    R_Bound = max(find(X_Sum > 0));
    if(R_Bound-widthmax>0)
        bw11=bw1(1:end, (R_Bound)-widthmax:end);
        XX_Sum=sum(bw11);
        L_Bound = min(find(XX_Sum>0))+R_Bound-widthmax;
    end
end

Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
Cropped_Image = imcrop(bw1, Rect);

if sum(sum(Cropped_Image))<image_minsize
        bw11=bw1(1:end, (min(find(X_Sum>0))+widthmax):(R_Bound)-widthmax);
        XX_Sum=sum(bw11);
        R_Bound=max(find(XX_Sum>0))+(min(find(X_Sum>0))+widthmax);
        L_Bound=min(find(XX_Sum>0))+(min(find(X_Sum>0))+widthmax);
end

Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
Cropped_Image = imcrop(bw1, Rect);


while sum(sum(Cropped_Image(1:end, 1:20)))<boundarythres
        bw11=bw1(T_Bound:B_Bound, (L_Bound+20):(L_Bound+20+widthmax));
        XX_Sum=sum(bw11);
        R_Bound=max(find(XX_Sum>0))+20+L_Bound;
        L_Bound=min(find(XX_Sum>0))+20+L_Bound;
        Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
        Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
        Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
        Cropped_Image = imcrop(bw1, Rect);
end

while sum(sum(Cropped_Image(1:end, end-19:end)))<boundarythres
        bw11=bw1(T_Bound:B_Bound, L_Bound:(R_Bound-20));
        XX_Sum=sum(bw11);
        R_Bound=max(find(XX_Sum>0))+L_Bound;
        Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
        Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
        Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
        Cropped_Image = imcrop(bw1, Rect);
end

while sum(sum(Cropped_Image(1:20, 1:end)))<boundarythres
        bw11=bw1((T_Bound+20):B_Bound, L_Bound:R_Bound);
        YY_Sum=sum(bw11');
        T_Bound=min(find(YY_Sum>0))+T_Bound+20;
        Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
        Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
        Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
        Cropped_Image = imcrop(bw1, Rect);
end


while sum(sum(Cropped_Image(end-19:end, 1:end)))<boundarythres
        bw11=bw1(T_Bound:(B_Bound-20), L_Bound:R_Bound);
        YY_Sum=sum(bw11');
        B_Bound=max(find(YY_Sum>0))+T_Bound;
        Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
        Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
        Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
        Cropped_Image = imcrop(bw1, Rect);
end


Weighted_Center = CenterMass(uint8(bw1));
%figure; imshow(Image_Original);
%hold on; plot(Weighted_Center(2), Weighted_Center(1), 'r*');
%hold on; plot(Graph_Center(1), Graph_Center(2),'bd');

Weighted_Center2 = CenterMass(uint8(Cropped_Image));
BW_path = strcat(Filepath,'\Process_ImgDrug\BWFast_', num2str(count1), '_', num2str(count2), '.jpg');
figure; imshow(Cropped_Image);

ima(:,:,1) = zeros(size(Cropped_Image));
ima(:,:,2) = 100.*Cropped_Image;
ima(:,:,3) = zeros(size(Cropped_Image));

figure; imshow(ima);
%imwrite(Cropped_Image, BW_path);
%hold on; plot(Weighted_Center2(2), Weighted_Center2(1), 'r*');
%hold on; plot(0.5*(R_Bound - L_Bound), 0.5*(B_Bound - T_Bound), 'bd');


%{
Cro_Orig = imcrop(Image_Original, Rect);
FL_path = strcat(Filepath,'\Process_FLDrug\FLFast_', num2str(count1), '_', num2str(count2), '.jpg');
figure; imshow(Cro_Orig);
imwrite(Cro_Orig, FL_path);
%hold on; plot(Weighted_Center2(2), Weighted_Center2(1), 'r*');
%hold on; plot(0.5*(R_Bound - L_Bound), 0.5*(B_Bound - T_Bound), 'bd');

BoundaryUp = xlsread('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', 1, strcat('B', num2str((count1)*2+count2-2)));
Upper_Img = Cropped_Image(1:(BoundaryUp-T_Bound), 1:end);
Lower_Img = Cropped_Image((BoundaryUp-T_Bound):end, 1:end);
U_path = strcat(Filepath,'\Process_UpperDrug\UFast_', num2str(count1), '_', num2str(count2), '.jpg');
L_path = strcat(Filepath,'\Process_LowerDrug\LFast_', num2str(count1), '_', num2str(count2), '.jpg');
UF_Img = Cro_Orig(1:(BoundaryUp-T_Bound), 1:end);
LF_Img = Cro_Orig((BoundaryUp-T_Bound):end, 1:end);
UF_path = strcat(Filepath,'\Process_Upper_FLDrug\UFFast_', num2str(count1), '_', num2str(count2), '.jpg');
LF_path = strcat(Filepath,'\Process_Lower_FLDrug\LFFast_', num2str(count1), '_', num2str(count2), '.jpg');
imwrite(Upper_Img, U_path);
imwrite(Lower_Img, L_path);
imwrite(UF_Img, UF_path);
imwrite(LF_Img, LF_path);
%}

figure;
imshow(bw1);

imag(:,:,1) = zeros(size(bw1));
imag(:,:,2) = 100.*bw1;
imag(:,:,3) = zeros(size(bw1));

figure; imshow(imag);
imwrite(imag, '\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\SampleImages\Label1.tif');
%{
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', L_Bound, 1, strcat('C', num2str((count1-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', R_Bound, 1, strcat('D', num2str((count1-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', T_Bound, 1, strcat('E', num2str((count1-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', B_Bound, 1, strcat('F', num2str((count1-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', Graph_Center(1), 1, strcat('G', num2str((count1-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', Graph_Center(2), 1, strcat('H', num2str((count1-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', Weighted_Center2(2)+L_Bound, 1, strcat('I', num2str((count1-1)*2+count2)));
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\FastData.xlsx', Weighted_Center2(1)+T_Bound, 1, strcat('J', num2str((count1-1)*2+count2)));
%}
    end
end