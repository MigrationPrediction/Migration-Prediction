clear all; close all; clc;
Filepath = '\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk';

threshold = 15;
threshold1 = 150;
threshold2 = 1000;

[filenameM, pathnameM, filterindexM] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';'*.*','All Files'},'Please select the multi-channel image', Filepath);
ImgMLocation = strcat(pathnameM, filenameM);
Image_Original = imread(ImgMLocation);
Imageperctile = prctile(Image_Original(:), 99);
if(Imageperctile < 5)
    threshold=9;
end
Image_G = Image_Original(100:end,1:end,2);

smooth = Image_G;
smooth = smooth .* uint8(smooth > threshold);
AreaFilter = bwareaopen(smooth, 10, 4);
Image_G = smooth .* uint8(AreaFilter);

Background = imopen(Image_G,strel('disk',2));
ProcessedImg = imadjust(Image_G-Background);
K = wiener2(ProcessedImg,[5 5]);

Centroid = regionprops(K,'Centroid');

level1 = graythresh(K)*0.7;
bw1 = im2bw(ProcessedImg, level1);
bw1 = bwareaopen(bw1, 15, 4);

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

Weighted_Center = CenterMass(uint8(bw1));
figure; imshow(Image_Original);
hold on; plot(Weighted_Center(2), Weighted_Center(1), 'r*');
hold on; plot(Graph_Center(1), Graph_Center(2),'bd');

Weighted_Center2 = CenterMass(uint8(Cropped_Image));
figure; imshow(Cropped_Image);
hold on; plot(Weighted_Center2(2), Weighted_Center2(1), 'r*');
hold on; plot(0.5*(R_Bound - L_Bound), 0.5*(B_Bound - T_Bound), 'bd');

[L, num] = bwlabel(bw1,4);
Dot_Area = zeros(size(bw1));
Intermediate_Area = zeros(size(bw1));
Fiber_Area = zeros(size(bw1));

figure; Mito_Num = 0;
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure19_T2\Label_Data.xlsx', 'Mito_Num', 1, 'A1');
xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure19_T2\Label_Data.xlsx', 'Mito_Label(1 as Elongate, 2 as Dotted, 3 as not sure)', 1, 'B1');
for i = 1:num
    flag = 0;
    DetectRegion = (L == i);
    s = sum(sum(DetectRegion));
    imshow(DetectRegion);
    Mito_Num = Mito_Num + 1;
    Mito_Label = 0;
    % Construct a question dialog with three options
    Choice = questdlg('What type of mitochondria is this?', ...
	'Mitochondria Image', ...
	'Elongated', 'Dotted','Not Sure', 'Not Sure');

    % Handle response
    switch Choice
        case 'Elongated'
            disp([Choice])
            Mito_Label = 1;
        case 'Dotted'
            disp([Choice])
            Mito_Label = 2;
        case 'Not Sure'
           flag = 1;
    end
    
    if flag == 1
         Choice2 = questdlg('If not surem, what type of mitochondria is this?', ...
            'Mitochondria Image', ...
            'Intermediate', 'Net', 'Exit','Intermediate');
            switch Choice2
                case 'Intermediate'
                    disp([Choice2])
                    Mito_Label = 3;
                case 'Net'
                    disp([Choice2])
                    Mito_Label = 4;
                case 'Exit'
                    pause;
            end
    end
    SaveLocation = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure19_T2\', num2str(Mito_Num), '.jpg');
    imwrite(DetectRegion, SaveLocation);
    xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure19_T2\Label_Data.xlsx', Mito_Num, 1, strcat('A', num2str(Mito_Num + 1)));
    xlswrite('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure19_T2\Label_Data.xlsx', Mito_Label, 1, strcat('B', num2str(Mito_Num + 1)));
end

%{
[m,n] = size(Image_G);
Final_image = zeros(m,n,3);
Final_image(:,:,1) = Final_image(:,:,1) + Fiber_Area .* 100;
Final_image(:,:,2) = Final_image(:,:,2) + Intermediate_Area .* 100;
Final_image(:,:,3) = Final_image(:,:,3) + Dot_Area .* 100;
savefig('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure19_T2\4_T1_Center.fig');
figure; imshow(Final_image);
imwrite(Final_image, '\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure19_T2\4_T1_MitoClass.jpg');


figure;
subplot(2,2,1);imshow(Image_G);title('Original Image');
subplot(2,2,2);imshow(ProcessedImg);title('Background Noise subtration');
subplot(2,2,3);imshow(K);title('Noise removal and contrast enhancement');
subplot(2,2,4);imshow(bw1);title('self-adjust thresholding');

figure;
subplot(2,2,1);surf(double(Image_Original(1:8:end,1:8:end)),'EdgeColor','none');
subplot(2,2,2);surf(double(K(1:8:end,1:8:end)),'EdgeColor','none');
subplot(2,2,3);surf(double(bw1(1:8:end,1:8:end)),'EdgeColor','none');
%}

Dot_sum = sum(sum(Dot_Area))
Inter_sum = sum(sum(Intermediate_Area))
Fiber_sum = sum(sum(Fiber_Area))