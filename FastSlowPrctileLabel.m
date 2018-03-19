load('absvalidspeed.mat');
prctile(absmove,[10 20 30 40 60 70 80 90]);


for i = 1:1358
    if absmove(i)<=0.0114
        label2080(i) = 0;
    elseif absmove(i)>=0.0519
        label2080(i) = 1;
    else
        label2080(i) =2;
    end
    
    if absmove(i)<=0.0064
        label1090(i) = 0;
    elseif absmove(i)>=0.0665
        label1090(i) = 1;
    else
        label1090(i) =2;
    end
    
    if absmove(i)<=0.0169
        label3070(i) = 0;
    elseif absmove(i)>=0.0407
        label3070(i) = 1;
    else
        label3070(i) =2;
    end
    
    if absmove(i)<=0.0224
        label4060(i) = 0;
    elseif absmove(i)>=0.0333
        label4060(i) = 1;
    else
        label4060(i) =2;
    end
end
