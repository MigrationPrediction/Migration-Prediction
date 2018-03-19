clear, clc
%{
CenterMove = load('CenterMoveRed.mat');
Dataup = load('Dataprocessupper.mat');
Datalower = load('Dataprocesslower.mat');
FiberRatio = load('Dataprocessimg');
Data = load('Data.mat');
Reddata = load('reddata.mat');
Serumcorrect = load('serumcorrectred.mat');
%Serumcorrect = load('serumwrong.mat');
updown = load('updownred.mat');
%updown = load('updownwrong.mat');
CenterMove = CenterMove.data';
CenterMove = CenterMove(1,:);
Dataup = Dataup.data;
Datalower = Datalower.data;
FiberRatio = FiberRatio.data;
Data = Data.data;
Reddata = Reddata.data;
Serumcorrect = Serumcorrect.serumcorrect;
updown = updown.updown;
%Serumcorrect = Serumcorrect.serumwrong;
%updown = updown.updownwrong;
%}
CenterMove = load('CenterMove.mat');
Dataup = load('Dataprocessupper.mat');
Datalower = load('Dataprocesslower.mat');
FiberRatio = load('Dataprocessimg');
Data = load('Data.mat');
Reddata = load('reddata.mat');
Serumcorrect = load('serumcorrect.mat');
%Serumcorrect = load('serumwrong.mat');
updown = load('updown.mat');
%updown = load('updownwrong.mat');
CenterMove = CenterMove.data';
CenterMove = CenterMove(1,:);
Dataup = Dataup.data;
Datalower = Datalower.data;
FiberRatio = FiberRatio.data;
Data = Data.data;
Reddata = Reddata.data;
Serumcorrect = Serumcorrect.serumcorrect;
updown = updown.updown;
%Serumcorrect = Serumcorrect.serumwrong;
%updown = updown.updownwrong;

no = 1;
rate = 0.7;
for i=[1011:3342]
%for i = 1362
    if(Serumcorrect(i) == 1)
        center = Data(2*i-1, 2);
        up = Data(2*i-1, 5);
        low = Data(2*i-1, 6);
        if updown(i) == 1
            stru{no}.centershift = log((center-up)/(low-center));
        elseif updown(i) == 2
            stru{no}.centershift = -log((center-up)/(low-center));
        end
        fiber1 = Dataup(i,2); 
        fiber2 = Datalower(i,2); 
        
        if fiber1 == 0
            fiber1 = 15;
        end
        if fiber2 == 0
            fiber2 = 15;
        end
        if updown(i) == 1
            stru{no}.fiberupdownratio = log(fiber1/fiber2);
        elseif updown(i) == 2
            stru{no}.fiberupdownratio = log(fiber2/fiber1);
        end
        width1 = Data(2*i-1,4)-Data(2*i-1,3);
        height1 = Data(2*i-1,6)-Data(2*i-1,5);
        stru{no}.aspectratiogreen = height1/width1;
        if i < 670
            width2 = Reddata(2*(i-609)-1,2)-Reddata(2*(i-609)-1,1);
            height2 = Reddata(2*(i-609)-1,4)-Reddata(2*(i-609)-1,3);
        else
            width2 = Reddata(2*(i-699)-1,2)-Reddata(2*(i-699)-1,1);
            height2 = Reddata(2*(i-699)-1,4)-Reddata(2*(i-699)-1,3);
        end
        BoundingBoxPerimeter = 2*(width2+height2);
        stru{no}.aspectratiored = height2/width2;
        temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_RG\RG', num2str(i),'_1');
        ImgMLocation = strcat(temp, '.jpg');
        Image_Original = imread(ImgMLocation);
        GreenBWLocation = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_Img\BW', num2str(i),'_1.jpg');
        ImgGreenBW = imread(GreenBWLocation);
        ImgGreenBW = im2bw(ImgGreenBW, 0.2);
        Image_R = Image_Original(1:end,1:end,1);
        Image_G = Image_Original(1:end,1:end,2);
        Image_R = Image_R .* uint8(Image_R > 7);
        Image_G = Image_G .* uint8(Image_G > 7);
        Image_R_bw = im2bw(Image_R, 0.07);
        %imshow(Image_R_bw);
        Image_G_bw = im2bw(Image_G, 0.07);
        redhighest = min(find(sum(Image_R_bw')>0))+Reddata(2*(i-699)-1,3);
        redlowest = max(find(sum(Image_R_bw')>0))+Reddata(2*(i-699)-1,3);
        greenhighest = min(find(sum(ImgGreenBW')>0))+Data(2*i-1,5);
        greenlowest = max(find(sum(ImgGreenBW')>0))+Data(2*i-1,5);
        stru{no}.redtotalarea = sum(sum(Image_R_bw));
        stru{no}.greentotalarea = sum(sum(ImgGreenBW));
        stru{no}.greenredarearatio = stru{no}.redtotalarea/stru{no}.greentotalarea;
        stru{no}.fiberratio = FiberRatio(i,2)/(FiberRatio(i,2)+FiberRatio(i,4)+FiberRatio(i,6));
        stru{no}.dotratio = FiberRatio(i,6)/(FiberRatio(i,2)+FiberRatio(i,4)+FiberRatio(i,6));
        
        [L, num] = bwlabel(Image_R_bw,4);
        maxarea = 1;
        maxno = 1;
        for j = 1:num
            DetectRegion = (L == j);
            temparea = sum(sum(DetectRegion));
            if temparea > maxarea
              maxno = j;
              maxarea = temparea;
            end
        end
        DetectRegion = (L == maxno);
        temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\R', num2str(i),'_1');
        ImgMLocationRed = strcat(temp, '.jpg');
        Image_RedOriginal = imread(ImgMLocationRed);
        ImageRed = Image_RedOriginal(:,:,1).*uint8(DetectRegion);
        GreenLocation = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_FL\FL', num2str(i),'_1.jpg');
        GreenImg = imread(GreenLocation);
        GreenImg = GreenImg(:,:,2).*uint8(ImgGreenBW);
        ImageRed(find(ImageRed(:)==0))=[];
        GreenImg(find(GreenImg(:)==0))=[];
        stru{no}.prctile99 = double(prctile(ImageRed, 99));
        stru{no}.prctile50 = double(prctile(ImageRed, 50));
        stru{no}.greenprctile50 = double(prctile(GreenImg, 50));
        stru{no}.prctile99prctile50ratio = stru{no}.prctile99/stru{no}.prctile50;
        stru{no}.totalintensityred = sum(ImageRed);
        stru{no}.totalintensitygreen = sum(GreenImg);
        stru{no}.totalintensityredgreenratio = stru{no}.totalintensitygreen/stru{no}.totalintensityred;
        if i < 670
            width3 = Reddata(2*(i-609)-1,1);
            height3 = Reddata(2*(i-609)-1,3);
        else
            width3 = Reddata(2*(i-699)-1,1);
            height3 = Reddata(2*(i-699)-1,3);
        end
        relativecenter1 = Data(2*i-1,1)-width3;
        relativecenter2 = Data(2*i-1,2)-height3;
        
        temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\R', num2str(i),'_1');
        ImgMLocationRed = strcat(temp, '.jpg');
        Image_RedOriginal = imread(ImgMLocationRed);
        ImageRed = Image_RedOriginal(:,:,1).*uint8(DetectRegion);
        centerintensityaver = double(sum(sum(ImageRed((relativecenter2-7):(relativecenter2+7), (relativecenter1-7):(relativecenter1+7))))/225);
        ImageRed = im2bw(ImageRed, 0.02);
        
        [L, num] = bwlabel(ImageRed,4);
        maxarea = 1;
        maxno = 1;
        for j = 1:num
            DetectRegionR = (L == j);
            temparea = sum(sum(DetectRegionR));
            if temparea > maxarea
              maxno = j;
              maxarea = temparea;
            end
        end
        ImageRed = (L == maxno);
        %figure; imshow(DetectRegion);
        %figure; imshow(ImageRed22);
        stru{no}.centerthickness = centerintensityaver/double(stru{no}.prctile50);
        stru{no}.convexarea = regionprops(ImageRed, 'ConvexArea');
        stru{no}.convexarea = stru{no}.convexarea.ConvexArea;
        stru{no}.solidity = regionprops(ImageRed, 'Solidity');
        stru{no}.solidity = stru{no}.solidity.Solidity;
        stru{no}.eccentricity = regionprops(ImageRed, 'Eccentricity');
        stru{no}.eccentricity = stru{no}.eccentricity.Eccentricity;
        stru{no}.equivdiameter = regionprops(ImageRed, 'Equivdiameter');
        stru{no}.equivdiameter = stru{no}.equivdiameter.EquivDiameter;
        stru{no}.extent = regionprops(ImageRed, 'Extent');
        stru{no}.extent = stru{no}.extent.Extent;
        stru{no}.majoraxis = regionprops(ImageRed, 'MajorAxisLength');
        stru{no}.majoraxis = stru{no}.majoraxis.MajorAxisLength;
        stru{no}.minoraxis = regionprops(ImageRed, 'MinorAxisLength');
        stru{no}.minoraxis = stru{no}.minoraxis.MinorAxisLength;
        stru{no}.orientation = regionprops(ImageRed, 'Orientation');
        stru{no}.orientation = abs(stru{no}.orientation.Orientation);
        stru{no}.perimeter = regionprops(ImageRed, 'Perimeter');
        stru{no}.perimeter = stru{no}.perimeter.Perimeter;
        center = Data(2*i-1, 2);
        if i < 670
            up = Reddata(2*(i-609)-1, 3);
            low = Reddata(2*(i-609)-1, 4);
        else
            up = Reddata(2*(i-699)-1, 3);
            low = Reddata(2*(i-699)-1, 4);
        end
        uplevel = int16(rate*up + (1-rate)*center-up+1);
        lowlevel = int16(rate*low + (1-rate)*center-up+1);
        
        
        
        %Problematic
        
        temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_R\R', num2str(i),'_1');
        ImgMLocationRed = strcat(temp, '.jpg');
        Image_RedOriginal = imread(ImgMLocationRed);
        Image_Red = Image_RedOriginal(:,:,1).*uint8(DetectRegion);
        
        %{
        temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_RG\RG', num2str(i),'_1');
        ImgMLocation = strcat(temp, '.jpg');
        Image_Original = imread(ImgMLocation);
        Image_Red = Image_Original(:,:,1);
        %}
        Bound_Up = max(find(Image_Red(1:uplevel,:) > 0))-min(find(Image_Red(1:uplevel,:) > 0));  
        Bound_Low = max(find(Image_Red(lowlevel:end,:) > 0))-min(find(Image_Red(lowlevel:end,:) > 0));
        Bound_Up = max(Bound_Up);
        Bound_Low = max(Bound_Low);
        if updown(i) == 1
            stru{no}.maxwidth = log(Bound_Up/Bound_Low);
        elseif updown(i) == 2
            stru{no}.maxwidth = -log(Bound_Up/Bound_Low);
        end
        Up = sum(Image_Red(1:uplevel,:)');   
        BoundingBox1 = double(2*(uplevel+width2));
        BoundingBox2 = double(2*(height2+width2-lowlevel));
        Low = sum(Image_Red(lowlevel:end,:)');
        Up = max(Up);
        Low = max(Low);
        if updown(i) == 1
            stru{no}.maxwidthsum = log(Up/Low);
        elseif updown(i) == 2
            stru{no}.maxwidthsum = -log(Up/Low);
        end
        
        %imshow(ImageRed);
        Area1 = sum(sum(ImageRed(1:uplevel,:)));
        Area2 = sum(sum(ImageRed(lowlevel:end,:)));
        if updown(i) == 1
            stru{no}.totalarearatio = log(Area1/Area2);
        elseif updown(i) == 2
            stru{no}.totalarearatio = -log(Area1/Area2);
        end
        [L, num] = bwlabel(Image_Red(1:uplevel,:),4);
        maxarea = 1;
        maxno = 1;
        for j = 1:num
            DetectRegion1 = (L == j);
            temparea = sum(sum(DetectRegion1));
            if temparea > maxarea
              maxno = j;
              maxarea = temparea;
            end
        end
        DetectRegion1 = (L == maxno);
        Perimeter1 = regionprops(DetectRegion1,'Perimeter');
        Perimeter1 = Perimeter1.Perimeter;
        [L, num] = bwlabel(Image_Red(lowlevel:end,:),4);
        maxarea = 1;
        maxno = 1;
        for j = 1:num
            DetectRegion2 = (L == j);
            temparea = sum(sum(DetectRegion2));
            if temparea > maxarea
              maxno = j;
              maxarea = temparea;
            end
        end
        DetectRegion2 = (L == maxno);
        Perimeter2 = regionprops(DetectRegion2,'Perimeter');
        Perimeter2 = Perimeter2.Perimeter;
        if updown(i) == 1
            stru{no}.totalperimeterratio = log(Perimeter1/Perimeter2);
        elseif updown(i) == 2
            stru{no}.totalperimeterratio = -log(Perimeter1/Perimeter2);
        end
        
        
        Convexarea1 = regionprops(DetectRegion1, 'ConvexArea');
        Convexarea2 = regionprops(DetectRegion2, 'ConvexArea');
        Convexarea1 = Convexarea1.ConvexArea;
        Convexarea2 = Convexarea2.ConvexArea;
        if updown(i) == 1
            stru{no}.topdownconvexxarearatio = log(Convexarea1/Convexarea2);
        elseif updown(i) == 2
            stru{no}.topdownconvexxarearatio = -log(Convexarea1/Convexarea2);
        end
        Solidity1 = regionprops(DetectRegion1, 'Solidity');
        Solidity2 = regionprops(DetectRegion2, 'Solidity');
        Solidity1 = Solidity1.Solidity;
        Solidity2 = Solidity2.Solidity;
        if updown(i) == 1
            stru{no}.topdownsolidityratio = log(Solidity1/Solidity2);
        elseif updown(i) == 2
            stru{no}.topdownsolidityratio = -log(Solidity1/Solidity2);
        end
        Eccentricity1 = regionprops(DetectRegion1, 'Eccentricity');
        Eccentricity2 = regionprops(DetectRegion2, 'Eccentricity');
        Eccentricity1 = Eccentricity1.Eccentricity;
        Eccentricity2 = Eccentricity2.Eccentricity;
        if updown(i) == 1
            stru{no}.topdowneccentricityratio = log(Eccentricity1/Eccentricity2);
        elseif updown(i) == 2
            stru{no}.topdowneccentricityratio = -log(Eccentricity1/Eccentricity2);
        end
        EquivDiameter1 = regionprops(DetectRegion1, 'EquivDiameter');
        EquivDiameter2 = regionprops(DetectRegion2, 'EquivDiameter');
        EquivDiameter1 = EquivDiameter1.EquivDiameter;
        EquivDiameter2 = EquivDiameter2.EquivDiameter;
        if updown(i) == 1
            stru{no}.topdownequivdiameterratio = log(EquivDiameter1/EquivDiameter2);
        elseif updown(i) == 2
            stru{no}.topdownequivdiameterratio = -log(EquivDiameter1/EquivDiameter2);
        end
        Extent1 = regionprops(DetectRegion1, 'Extent');
        Extent2 = regionprops(DetectRegion2, 'Extent');
        Extent1 = Extent1.Extent;
        Extent2 = Extent2.Extent;
        if updown(i) == 1
            stru{no}.topdownextentratio = log(Extent1/Extent2);
        elseif updown(i) == 2
            stru{no}.topdownextentratio = -log(Extent1/Extent2);
        end
        MajorAxisLength1 = regionprops(DetectRegion1, 'MajorAxisLength');
        MajorAxisLength2 = regionprops(DetectRegion2, 'MajorAxisLength');
        MajorAxisLength1 = MajorAxisLength1.MajorAxisLength;
        MajorAxisLength2 = MajorAxisLength2.MajorAxisLength;
        if updown(i) == 1
            stru{no}.topdownmajoraxislengthratio = log(MajorAxisLength1/MajorAxisLength2);
        elseif updown(i) == 2
            stru{no}.topdownmajoraxislengthratio = -log(MajorAxisLength1/MajorAxisLength2);
        end
        MinorAxisLength1 = regionprops(DetectRegion1, 'MinorAxisLength');
        MinorAxisLength2 = regionprops(DetectRegion2, 'MinorAxisLength');
        MinorAxisLength1 = MinorAxisLength1.MinorAxisLength;
        MinorAxisLength2 = MinorAxisLength2.MinorAxisLength;
        if updown(i) == 1
            stru{no}.topdownminoraxislengthratio = log(MinorAxisLength1/MinorAxisLength2);
        elseif updown(i) == 2
            stru{no}.topdownminoraxislengthratio = -log(MinorAxisLength1/MinorAxisLength2);
        end
        
        
        
        
        
        %imshow(Image_Red);
        
        
        
        
        
        ImageRed1 = Image_Red(1:uplevel,:).*uint8(DetectRegion1);
        ImageRed1(find(ImageRed1(:)==0))=[];
        ImageRed2 = Image_Red(lowlevel:end,:).*uint8(DetectRegion2);
        ImageRed2(find(ImageRed2(:)==0))=[];
        if updown(i) == 1
            stru{no}.topdownprctile99 = log(double(prctile(ImageRed1, 99))/double(prctile(ImageRed2, 99)));
            stru{no}.topdownprctile50 = log(double(prctile(ImageRed1, 50))/double(prctile(ImageRed2, 50)));
            stru{no}.topdownprctile99prctile50ratio = log(double(prctile(ImageRed1, 99))*double(prctile(ImageRed2, 50))/double(prctile(ImageRed2, 99))/double(prctile(ImageRed1, 50)));
            stru{no}.topdowntotalintensityratio = log(sum(sum(ImageRed1))/sum(sum(ImageRed2)));
            stru{no}.topdownperimeterarearatio = log((Perimeter1/Area1)/(Perimeter2/Area2));
            stru{no}.topdownperimeterboundingboxratio = log((Perimeter1/BoundingBox1)/(Perimeter2/BoundingBox2));
            stru{no}.topdownmeanintensityratio = log((sum(sum(ImageRed1))/Area1)/(sum(sum(ImageRed2))/Area2));
        elseif updown(i) == 2
            stru{no}.topdownprctile99 = -log(double(prctile(ImageRed1, 99))/double(prctile(ImageRed2, 99)));
            stru{no}.topdownprctile50 = -log(double(prctile(ImageRed1, 50))/double(prctile(ImageRed2, 50)));
            stru{no}.topdownprctile99prctile50ratio = -log(double(prctile(ImageRed1, 99))*double(prctile(ImageRed2, 50))/double(prctile(ImageRed2, 99))/double(prctile(ImageRed1, 50)));
            stru{no}.topdowntotalintensityratio = -log(sum(sum(ImageRed1))/sum(sum(ImageRed2)));
            stru{no}.topdownperimeterarearatio = -log((Perimeter1/Area1)/(Perimeter2/Area2));
            stru{no}.topdownperimeterboundingboxratio = -log((Perimeter1/BoundingBox1)/(Perimeter2/BoundingBox2));
            stru{no}.topdownmeanintensityratio = -log((sum(sum(ImageRed1))/Area1)/(sum(sum(ImageRed2))/Area2));
        end
        
        stru{no}.majorminoraxisratio = stru{no}.majoraxis/stru{no}.minoraxis;
        stru{no}.perimeterarearatio = stru{no}.perimeter/stru{no}.redtotalarea;
        stru{no}.perimeterboundingboxratio = stru{no}.perimeter/BoundingBoxPerimeter;
        stru{no}.meanintensity = stru{no}.totalintensityred/stru{no}.redtotalarea;
        stru{no}.perimetersquarearearatio = stru{no}.perimeter*stru{no}.perimeter/stru{no}.redtotalarea;
        
        
        
        [L, num] = bwlabel(Image_Red(1:(center-up+1),:),4);
        maxarea = 1;
        maxno = 1;
        for j = 1:num
            DetectRegion1 = (L == j);
            temparea = sum(sum(DetectRegion1));
            if temparea > maxarea
              maxno = j;
              maxarea = temparea;
            end
        end
        DetectRegion1 = (L == maxno);
        Perimeter1 = regionprops(DetectRegion1,'Perimeter');
        Perimeter1 = Perimeter1.Perimeter;
        [L, num] = bwlabel(Image_Red((center-up+1):end,:),4);
        maxarea = 1;
        maxno = 1;
        for j = 1:num
            DetectRegion2 = (L == j);
            temparea = sum(sum(DetectRegion2));
            if temparea > maxarea
              maxno = j;
              maxarea = temparea;
            end
        end
        DetectRegion2 = (L == maxno);
        Perimeter2 = regionprops(DetectRegion2,'Perimeter');
        Perimeter2 = Perimeter2.Perimeter;
        if updown(i) == 1
            stru{no}.totalperimeternewratio = log(Perimeter1/Perimeter2);
        elseif updown(i) == 2
            stru{no}.totalperimeternewratio = -log(Perimeter1/Perimeter2);
        end
        
        newArea1 = sum(sum(DetectRegion1));
        newArea2 = sum(sum(DetectRegion2));
        
        
        imshow(Image_Red);
        redsize = size(ImageRed);
        averagewid = 0;
        for j = 1:redsize(1)
            if max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0))>0
                averagewid = averagewid + max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0));
            end
        end
        stru{no}.averagewidth = averagewid/redsize(1);
        leftbound = sum(ImageRed(:,10));
        rightbound = sum(ImageRed(:,end-9));
        leftboundratio = leftbound/redsize(1);
        rightboundratio = rightbound/redsize(1);
        stru{no}.leftrightbound = max(leftbound, rightbound);
        stru{no}.leftrightboundratio = max(leftboundratio,rightboundratio);
        headaveragewid = 0;
        headcount = int32(redsize(1)/6);
        if updown(i) == 1
            for j = 1:headcount
                if max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0))>0
                    headaveragewid = headaveragewid + double(max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0)))/double(headcount);
                end
            end
        elseif updown(i) == 2
            for j = (redsize(1)-headcount+1):redsize(1)
                if max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0))>0
                    headaveragewid = headaveragewid + double(max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0)))/double(headcount);
                end
            end
        end
        stru{no}.headaveragewidth = headaveragewid;
        stru{no}.headaveragewidthratio = double(stru{no}.headaveragewidth)/double(stru{no}.averagewidth);
        
        maxwidthred = 0;
        maxwidthredno = 0;
        for j = 1:redsize(1)
            if (max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0)))>maxwidthred
                maxwidthred = max(find(ImageRed(j,:) > 0))-min(find(ImageRed(j,:) > 0));
                maxwidthredno = j;
            end
        end
        stru{no}.centerred = double(max(find(ImageRed(relativecenter2,:) > 0))-min(find(ImageRed(relativecenter2,:) > 0)))/double(maxwidthred);
        
        if updown(i) == 1
            stru{no}.redshift = log((maxwidthredno-1)/(redsize(1)-maxwidthredno));
        elseif updown(i) == 2
            stru{no}.redshift = -log((maxwidthredno-1)/(redsize(1)-maxwidthredno));
        end
        
        sumintensity = double(0);
        if updown(i) == 1
            for j = 13:(redsize(2)-12)
                redmin = min(find(ImageRed(:,j) > 0));
                for k = 1:25
                    addred = double(Image_Red(redmin-1+k,j));
                    if addred > 0
                        sumintensity = sumintensity + addred;
                    end
                end
            end
            stru{no}.redfoot = sumintensity/25/(redsize(2)-24)/centerintensityaver;
        elseif updown(i) == 2
            for j = 13:(redsize(2)-12)
                redmin = max(find(ImageRed(:,j) > 0));
                for k = 1:25
                    addred = double(Image_Red(redmin+1-k,j));
                    if addred > 0
                        sumintensity = sumintensity + addred;
                    end
                end
            end
            stru{no}.redfoot = sumintensity/25/(redsize(2)-24)/centerintensityaver;
        end
        
        temp = strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\Cells\Process_Img\BW', num2str(i),'_1');
        ImgMLocationBW = strcat(temp, '.jpg');
        Image_BW = imread(ImgMLocationBW);
        %figure;imshow(Image_BW);
        Image_BW = im2bw(Image_BW, 0.1);
        %figure;imshow(Image_BW);
        cen1 = Data(2*i-1, 2);
        cen2 = Data(2*i-1, 1);
        lef = Data(2*i-1, 3);
        rig = Data(2*i-1,4);
        up = Data(2*i-1, 5);
        down = Data(2*1-1,6);
        relativecenter3 = cen2-lef+1;
        relativecenter4 = cen1-up+1;
        if relativecenter3-20 < 1
            lef1 = 1;
        else
            lef1 = relativecenter3-20;
        end
        if relativecenter3+20 > (rig-lef+1)
            rig1 = rig-lef+1;
        else
            rig1 = relativecenter3+20;
        end
        if relativecenter4-20 < 1
            top1 = 1;
        else
            top1 = relativecenter4-20;
        end
        if relativecenter4+20 > (down-up+1)
            down1 = down-up+1;
        else
            down1 = relativecenter4+20;
        end
        
        mitoarea = sum(sum(Image_BW(top1:down1,lef1:rig1)));
        stru{no}.centermito = double(mitoarea/(rig1-lef1+1)/(down1-top1+1));
        
        if updown(i) == 1
            stru{no}.redgreendist = (greenhighest-redhighest)/height2;
        elseif updown(i) == 2
            stru{no}.redgreendist = (redlowest-greenlowest)/height2;
        end
        
        if updown(i) == 1
            stru{no}.totalareanewratio = log(newArea1/newArea2);
        elseif updown(i) == 2
            stru{no}.totalareanewratio = -log(newArea1/newArea2);
        end
        
        
        
        sumintensity1 = double(0);
        sumintensity2 = double(0);
        if updown(i) == 1
            for j = 13:(redsize(2)-12)
                redmin = min(find(ImageRed(:,j) > 0));
                for k = 1:24
                    addred = double(Image_Red(redmin-1+k,j));
                    if addred > 0
                        sumintensity1 = sumintensity1 + addred;
                    end
                end
            end
            
            for j = 13:(redsize(2)-12)
                redmin = max(find(ImageRed(:,j) > 0));
                for k = 1:24
                    addred = double(Image_Red(redmin+1-k,j));
                    if addred > 0
                        sumintensity2 = sumintensity2 + addred;
                    end
                end
            end
            
            stru{no}.redfootupdown = sumintensity1/sumintensity2;
        elseif updown(i) == 2
            for j = 13:(redsize(2)-12)
                redmin = max(find(ImageRed(:,j) > 0));
                for k = 1:24
                    addred = double(Image_Red(redmin+1-k,j));
                    if addred > 0
                        sumintensity1 = sumintensity1 + addred;
                    end
                end
            end
            
            for j = 13:(redsize(2)-12)
                redmin = min(find(ImageRed(:,j) > 0));
                for k = 1:24
                    addred = double(Image_Red(redmin-1+k,j));
                    if addred > 0
                        sumintensity2 = sumintensity2 + addred;
                    end
                end
            end
            
            stru{no}.redfootupdown = sumintensity1/sumintensity2;
        end
        
        no = no+1;
    end
end