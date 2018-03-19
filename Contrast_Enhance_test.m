clear all; close all; clc;
Filepath = 'P:\Research\20170510_SUM159_Migration_Lili';


threshold1 = 150;
threshold2 = 1000;

[filenameM, pathnameM, filterindexM] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'},'Please select the multi-channel image', Filepath);
ImgMLocation = strcat(pathnameM, filenameM);
Image_Original = imread(ImgMLocation);
Image_G = Image_Original(1:end,1:end,2);

[counts,binLocations] = imhist(Image_G);
diffcounts = counts(1:end-1)-counts(2:end);
threshold = find(diffcounts == max(diffcounts));

%threshold = 7;
smooth = Image_G;
smooth = smooth .* uint8(smooth > threshold);
AreaFilter = bwareaopen(smooth, 10, 8);
Image_G = smooth .* uint8(AreaFilter);

Background = imopen(Image_G,strel('disk',3));
ProcessedImg = imadjust(Image_G-Background);
K = wiener2(ProcessedImg,[5 5]);

Centroid = regionprops(K,'Centroid');

level1 = graythresh(K)*0.7;
bw1 = im2bw(ProcessedImg, level1);
bw1 = bwareaopen(bw1, 15, 8);

X_Sum = sum(bw1);
Y_Sum = sum(bw1');
L_Bound = min(find(X_Sum > 0));
R_Bound = max(find(X_Sum > 0));
T_Bound = min(find(Y_Sum > 0));
B_Bound = max(find(Y_Sum > 0));
Bounding_Box = [L_Bound, R_Bound, T_Bound, B_Bound];
Rect =[L_Bound, T_Bound, R_Bound - L_Bound, B_Bound - T_Bound];
Graph_Center = [L_Bound + 0.5*(R_Bound - L_Bound), T_Bound + 0.5*(B_Bound - T_Bound)];
Cropped_Image = imcrop(bw1, Rect);


%{
[m,n] = size(Image_G);
Final_image = zeros(m,n,3);
Final_image(:,:,1) = Final_image(:,:,1) + Fiber_Area .* 100;
Final_image(:,:,2) = Final_image(:,:,2) + Intermediate_Area .* 100;
Final_image(:,:,3) = Final_image(:,:,3) + Dot_Area .* 100;
%%savefig('D:\Research\20170504_SUM159_Mito_Migration\Process_Img\4_T1_Center.fig');
figure; imshow(Final_image);
%%imwrite(Final_image, 'D:\Research\20170504_SUM159_Mito_Migration\Process_Img\4_T1_MitoClass.jpg');
%}
figure; imshow(Image_Original);

figure;
subplot(2,2,1);imshow(Image_G);title('Original Image');
subplot(2,2,2);imshow(ProcessedImg);title('Background Noise subtration');
subplot(2,2,3);imshow(K);title('Noise removal and contrast enhancement');
subplot(2,2,4);imshow(bw1);title('self-adjust thresholding');
%{
figure;
subplot(2,2,1);surf(double(Image_Original(1:8:end,1:8:end)),'EdgeColor','none');
subplot(2,2,2);surf(double(K(1:8:end,1:8:end)),'EdgeColor','none');
subplot(2,2,3);surf(double(bw1(1:8:end,1:8:end)),'EdgeColor','none');


Dot_sum = sum(sum(Dot_Area))
Inter_sum = sum(sum(Intermediate_Area))
Fiber_sum = sum(sum(Fiber_Area))
%}