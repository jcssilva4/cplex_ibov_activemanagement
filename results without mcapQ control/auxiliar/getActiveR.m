function activeR = getActiveR(d, P)
    assets = P.initialAssetList;
    activeR = 0;
    for i=1:length(P.initialAssetList)  
        activeR = activeR + d(i)*(assets{1,i}.alphaScore); 
    end
end