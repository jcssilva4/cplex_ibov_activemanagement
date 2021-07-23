function [Ropt] = getOverallR ( w, P)
    somaWR = 0;
    for i=1:length(P.initialAssetList)
        somaWR = somaWR + w(i)*P.initialAssetList{1,i}.r;   
    end
    Ropt = somaWR;
end