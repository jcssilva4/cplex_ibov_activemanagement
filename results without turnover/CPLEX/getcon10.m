function [Acon10, bcon10] = getcon10(x, P)
    assets = P.initialAssetList;
    Acon10 = zeros(2,2*length(assets));
    bcon10 = zeros(2,1);
    %r10.1 
    bcon10(1,1) = 0.4;
    for i = 1:length(assets)
        if(x(i) < assets{1,i}.benchWeight) %min(wi,wb_i)
            Acon10(1,i) = 1;
        else %the min is bench
            bcon10(1,1) = bcon10(1,1) - assets{1,i}.benchWeight;
        end
    end
    %r10.2
    bcon10(2,1) = 0;
    for i = 1:length(assets)
        if(x(i) < assets{1,i}.benchWeight) %min(wi,wb_i)
            Acon10(2,i) = -1;
        else %the min is bench
            bcon10(2,1) = bcon10(2,1) + assets{1,i}.benchWeight;
        end
    end
end