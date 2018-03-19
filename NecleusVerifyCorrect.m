%clear all, clc
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
temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Nucleus_Label_Check\ChamberImg', num2str(count1*2-1, '%04i'));
    ImgMLocation = strcat(temp, '.tif');
    graythre = 0.1;
    
Image_Original = imread(ImgMLocation);
Imageperctile = prctile(Image_Original(:), 99);
%Image_G = Image_Original(1:end,1:end,2);
Image_G = Image_Original;

[counts,binLocations] = imhist(Image_G);
diffcounts = counts(1:end-1)-counts(2:end);
threshold = find(diffcounts == max(diffcounts(4:end)));

threshold = 20;

Image_New = im2bw(Image_G, 0.2);

Nucleus = bwselect(Image_New, 243, 217, 4);

figure;imshow(Image_G);

image1 = imread('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\SampleImages\Label1.tif');

image1(:,:,3) = 180.*Image_New;

figure;imshow(Image_New);

figure;imshow(image1);

image2 = imcrop(image1, Rect);

figure; imshow(image2);
imwrite(image2, '\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\SampleImages\NucleusLabelling.tif');
Weighted_Center = CenterMass(uint8(Nucleus));
    end
end