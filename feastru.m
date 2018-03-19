clear, clc
E=load('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure25_T1\result');
F=size(E.B);
for i=1:F(2)
fig=imread(strcat('\\engin-labs.m.storage.umich.edu\chenlili\windat.v2\Desktop\label_smallerdisk\figure25_T1\',num2str(i), '.jpg'));
stru=Morphology_Feature(fig);
FEA{i}=[stru.CenterShift; stru.MajorAxis; stru.MinorAxis; stru.AspectRatio; stru.Area; stru.Extent; stru.Perimeter; ...
        stru.Orientation; stru.Eccentricity; stru.EulerNumber; stru.EquivDiameter; stru.ConvexArea; stru.Solidity; ...
        stru.BrightFieldEntropy; stru.FluoEntropy; stru.FluoContrast];
end