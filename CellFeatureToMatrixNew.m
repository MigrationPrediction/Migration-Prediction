clear all, clc;
%CellFeature = load('FeatureAbs.mat');
%CellFeature = load('CellFeatureMatrixNew.mat');%61 Cell use
CellFeature = load('NewCellFeature.mat');
%CellFeature = load('CellFeatureRed.mat');
CellFeature = CellFeature.stru;
%move = load('absvalidmove.mat');
%move = move.move;
nooo = 1;

for i = 1:1358
    CellFeatureMatrix(:,i) = [CellFeature{i}.centershift, CellFeature{i}.fiberupdownratio, CellFeature{i}.aspectratiogreen, CellFeature{i}.aspectratiored, ...
        CellFeature{i}.redtotalarea, CellFeature{i}.greentotalarea, CellFeature{i}.greenredarearatio, CellFeature{i}.fiberratio, CellFeature{i}.dotratio, ...
        CellFeature{i}.prctile99, CellFeature{i}.prctile50, CellFeature{i}.greenprctile50, CellFeature{i}.prctile99prctile50ratio, CellFeature{i}.totalintensityred, ...
        CellFeature{i}.totalintensitygreen, CellFeature{i}.totalintensityredgreenratio, CellFeature{i}.centerthickness, CellFeature{i}.convexarea, ...
        CellFeature{i}.solidity, CellFeature{i}.eccentricity, CellFeature{i}.equivdiameter, CellFeature{i}.extent, CellFeature{i}.majoraxis, ...
        CellFeature{i}.minoraxis, CellFeature{i}.orientation, CellFeature{i}.perimeter, CellFeature{i}.maxwidth, CellFeature{i}.maxwidthsum, ...
        CellFeature{i}.totalarearatio, CellFeature{i}.totalperimeterratio, CellFeature{i}.topdownconvexxarearatio, CellFeature{i}.topdownsolidityratio, ...
        CellFeature{i}.topdowneccentricityratio, CellFeature{i}.topdownequivdiameterratio, CellFeature{i}.topdownextentratio, CellFeature{i}.topdownmajoraxislengthratio, ...
        CellFeature{i}.topdownminoraxislengthratio, CellFeature{i}.topdownprctile99, CellFeature{i}.topdownprctile50, CellFeature{i}.topdownprctile99prctile50ratio, ...
        CellFeature{i}.topdowntotalintensityratio, CellFeature{i}.topdownperimeterarearatio, CellFeature{i}.topdownperimeterboundingboxratio, ...
        CellFeature{i}.topdownmeanintensityratio, CellFeature{i}.majorminoraxisratio, CellFeature{i}.perimeterarearatio, CellFeature{i}.perimeterboundingboxratio, CellFeature{i}.meanintensity, ...
        CellFeature{i}.perimetersquarearearatio, CellFeature{i}.totalperimeternewratio, CellFeature{i}.averagewidth, CellFeature{i}.leftrightbound, ...
        CellFeature{i}.leftrightboundratio, CellFeature{i}.headaveragewidth, CellFeature{i}.headaveragewidthratio, CellFeature{i}.centerred, ...
        CellFeature{i}.redshift, CellFeature{i}.redfoot, CellFeature{i}.redfootupdown, CellFeature{i}.redgreendist, CellFeature{i}.totalareanewratio];    
    
end

%{
for i = 1:1358
    if CellFeatureMatrix(2,i) > 4
        CellFeatureMatrix(2,i) = 4;
    elseif CellFeatureMatrix(2,i) < -3
        CellFeatureMatrix(2,i) = -3;
    end
    
    if CellFeatureMatrix(3,i) > 10
        CellFeatureMatrix(3,i) = 10;
    end
    
    if CellFeatureMatrix(4,i) > 7.5
        CellFeatureMatrix(4,i) = 7.5;
    end
    
    if CellFeatureMatrix(5,i) > 102504
        CellFeatureMatrix(5,i) = 102504;
    end
    
    if CellFeatureMatrix(6,i) > 17569
        CellFeatureMatrix(6,i) = 17569;
    end
    
    if CellFeatureMatrix(7,i) > 30
        CellFeatureMatrix(7,i) = 30;
    end
    
    if CellFeatureMatrix(10,i) > 120
        CellFeatureMatrix(10,i) = 120;
    end
    
    if CellFeatureMatrix(12,i) > 50
        CellFeatureMatrix(12,i) = 50;
    end
    
    if CellFeatureMatrix(15,i) > 543925
        CellFeatureMatrix(15,i) = 543925;
    end
    
    if CellFeatureMatrix(16,i) > 0.60345
        CellFeatureMatrix(16,i) = 0.60345;
    end
    
    if CellFeatureMatrix(25,i) < 57
        CellFeatureMatrix(25,i) = 57;
    end
    
    if CellFeatureMatrix(33,i) > 1
        CellFeatureMatrix(33,i) = 1;
    elseif CellFeatureMatrix(33,i) < -1
        CellFeatureMatrix(33,i) = -1;
    end
    
    if CellFeatureMatrix(35,i) > 0.8
        CellFeatureMatrix(35,i) = 0.8;
    end
    
    if CellFeatureMatrix(45,i) > 10.32
        CellFeatureMatrix(45,i) = 10.32;
    end
end
%}

no = 61;
%{
for i = 1:no
%for i = 60:60
    figure;boxplot(CellFeatureMatrix(i,:));
end

for i = 1:no
%for i = 60:60
    A = CellFeatureMatrix(i,:);
    min1 = min(A);
    max1 = max(A);
    A = 1+9.*((A-min1)./(max1-min1));
    NormalizedFeatureMatrix(i,:) = A;
    figure;boxplot(A);
end
%}

serumcorrect = load('serumcorrect.mat');
serumcorrect = serumcorrect.serumcorrect;

load('absspeed.mat');

centermove = load('CenterMove.mat');
centermove = centermove.data;
noo = 1;
for i = 1:3342
    if serumcorrect(i) == 1 
        %figure; imshow(strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_RG\RG',num2str(i),'_1.jpg'));
        move(noo) = centermove(i);
        mov(noo) = absspeed(i);
        noo = noo +1;
    end
end
absmove = abs(move);
ratio = CellFeatureMatrix(42,:);
%figure; scatter(ratio,move);


for i = 1:no
    a(i,:)=prctile(CellFeatureMatrix(i,:),[33 50 67]);
end