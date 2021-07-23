function f = getobjfnVal(x, P, lambda)  

    f = 0; %initialize
    nAssets = length(P.initialAssetList);
    assets = P.initialAssetList;
    %di = wi-wbenchi
    di = 0;
    dj = 0;
    %d = x(1:nAssets) - assets{1,:}.benchWeight; %compute di
    %(wi-wbench)*covij*(wj-wbench)
    cov = 0;
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); % get covij + covji
        else
            cov = P.covij(t); %get vari
        end
        di = x(P.i_asset(t))-assets{1,P.i_asset(t)}.benchWeight;
        dj = x(P.j_asset(t))-assets{1,P.j_asset(t)}.benchWeight;
        f= f + di*cov*dj;
    end
    %-lambda*d*alpha_score
    for t=1:nAssets  %loop over assets vector
        f= f - (lambda*(x(t) - assets{1,t}.benchWeight)*assets{1,t}.alphaScore);
    end
end