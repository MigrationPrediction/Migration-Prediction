function [Morph_Feature] = Morphology_Feature(Image_Fluo) 

%% Select regions using the detected center from cell counting function
%Image_Fluo_BW = bwselect(Image_Fluo, CellCenterX, CellCenterY, 8);
flag = 0;
Image_Fluo_BW = im2bw(Image_Fluo, 0.5);
Morph_Feature = struct('CenterShift', [], 'MajorAxis', [], 'MinorAxis', [], ...
'AspectRatio', [], 'Area', [], 'Extent', [], 'Perimeter', [], 'Orientation', [], 'Eccentricity', [], 'EulerNumber', [], 'EquivDiameter', [], ...
'ConvexArea', [], 'Solidity', [], 'BrightFieldEntropy', [], 'FluoEntropy', [], 'FluoContrast', []);

if nnz(Image_Fluo_BW) == 0
    flag = 1;
    return;
end

%% Extract morphology features
MajorAxis = regionprops(Image_Fluo_BW, 'MajorAxisLength');
MinorAxis = regionprops(Image_Fluo_BW, 'MinorAxisLength');
if flag == 1
    a = 1;
end
AspectRatio = MajorAxis.MajorAxisLength ./ MinorAxis.MinorAxisLength;
Area = regionprops(Image_Fluo_BW, 'Area');
Extent = regionprops(Image_Fluo_BW, 'Extent');
Perimeter = regionprops(Image_Fluo_BW, 'Perimeter');
Orientation = regionprops(Image_Fluo_BW, 'Orientation');
Eccentricity = regionprops(Image_Fluo_BW, 'Eccentricity');
EulerNumber = regionprops(Image_Fluo_BW, 'EulerNumber');
EquivDiameter = regionprops(Image_Fluo_BW, 'EquivDiameter');
ConvexArea = regionprops(Image_Fluo_BW, 'ConvexArea');
Solidity = regionprops(Image_Fluo_BW, 'Solidity');
%MeanIntensity = regionprops(Image_Fluo_BW, Image_Fluo, 'MeanIntensity');
%MaxIntensity = regionprops(Image_Fluo_BW, Image_Fluo, 'MaxIntensity');
PixelIntensity = regionprops(Image_Fluo_BW, Image_Fluo, 'PixelValues');
%TotalIntensity = sum(PixelIntensity.PixelValues);
        
% Center Difference between weighted center and graphical center
WeightedCenter = regionprops(Image_Fluo_BW, Image_Fluo, 'WeightedCentroid');
RegionCenter = regionprops(Image_Fluo_BW,'Centroid');
Center1 = RegionCenter.Centroid;
Center2 = WeightedCenter.WeightedCentroid;
CenterShift = sqrt((Center1(1)-Center2(1)) .^2 + (Center1(2)-Center2(2)) .^2);


%% Extract Texture features
%Cell_Gray_Image = rgb2gray(Image_BF) .* uint8(Image_Fluo_BW);
Cell_Gray_Image = uint8(Image_Fluo_BW);

Entropy_BF = entropy(Cell_Gray_Image);
Entropy_Fluo = entropy(Image_Fluo);
GLCM_BF = graycomatrix(Cell_Gray_Image);
GLCM_Fluo = graycomatrix(Image_Fluo);
Texture_Stats_BF = graycoprops(GLCM_BF);
Texture_Stats_Fluo = graycoprops(GLCM_Fluo);


Morph_Feature = struct('CenterShift', CenterShift, 'MajorAxis', MajorAxis.MajorAxisLength, 'MinorAxis', MinorAxis.MinorAxisLength, ...
'AspectRatio', AspectRatio, 'Area', Area.Area, 'Extent', Extent.Extent, 'Perimeter', Perimeter.Perimeter, 'Orientation', Orientation.Orientation, 'Eccentricity', ...
Eccentricity.Eccentricity, 'EulerNumber', EulerNumber.EulerNumber, 'EquivDiameter', EquivDiameter.EquivDiameter, 'ConvexArea', ConvexArea.ConvexArea, ...
'Solidity', Solidity.Solidity, 'BrightFieldEntropy', Entropy_BF, ...
'FluoEntropy', Entropy_Fluo, 'FluoContrast', Texture_Stats_Fluo.Contrast);


%{

Boundary = bwboundaries(Image_Fluo_BW);
Aspect_Ratio = zeros(length(MajorAxis), 1);
CellNum = length(RegionCenter);
Total_Intensity = zeros(CellNum,1);
Num_Marker = cell(CellNum,1);
Cell_Position = zeros(CellNum,2);
Cell_Area = zeros(CellNum,1);
Counter = 0;
for i = 1:CellNum
    if Area(i).Area > 1000
        Counter = Counter + 1;
        Total_Intensity(Counter) = sum(PixelIntensity(i).PixelValues);
        Aspect_Ratio(Counter) = MajorAxis(i).MajorAxisLength ./ MinorAxis(i).MinorAxisLength;
        Cell_Position(Counter,:) = RegionCenter(i).Centroid;
        Num_Marker{Counter} = num2str(i);
        Cell_Area(Counter) = Area(i).Area;
        Mean_Intensity(Counter) = MeanIntensity(i).MeanIntensity;
    end
end

CellNum = Counter;
Cell_Position(Counter + 1 : end,:) = [];
Num_Marker(Counter + 1 : end) = [];
Aspect_Ratio(Counter + 1 : end) = [];
Total_Intensity(Counter + 1 : end) = [];
Cell_Area(Counter + 1 : end) = [];
Mean_Intensity(Counter + 1 : end) = [];



I = insertText(ProcessedImg,Cell_Position,Num_Marker,'TextColor','blue','FontSize',18,'BoxOpacity',0,'AnchorPoint','Center');
figure;
imshow(I);hold on;


%figure;
%imshow(imresize(ProcessedImg,AmpFactor));
%hold on;
for i = 1:length(Boundary)
    plot(Boundary{i}(:,2),Boundary{i}(:,1),'-');
end

I2 = insertText(ProcessedImgFl_G,Cell_Position,Num_Marker,'TextColor','blue','FontSize',18,'BoxOpacity',0,'AnchorPoint','Center');
figure;
imshow(I2);hold on;
%}

