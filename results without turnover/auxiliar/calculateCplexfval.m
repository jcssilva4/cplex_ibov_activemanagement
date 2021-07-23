function f = calculateCplexfval(x, P, lambda)
    f = 0; %initialize
    nAssets = length(P.initialAssetList);
    assets = P.initialAssetList;
    %d = x(1:nAssets) - assets{1,:}.benchWeight; %compute di
    %(wi)*covij*(wj)
    cov = 0;
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); % get covij + covji
        else
            cov = P.covij(t); %get vari
        end
        f = f + x(P.i_asset(t))*cov*x(P.j_asset(t));
        %-covij*wi*wbj - covij*wj*wbi 
        f = f - (cov*(x(P.i_asset(t))*assets{1,P.j_asset(t)}.benchWeight + x(P.j_asset(t))*assets{1,P.i_asset(t)}.benchWeight));
    end
        
    %-lambda*d*alpha_score
    for t=1:nAssets  %loop over assets vector
        f= f - (lambda*(x(t)*assets{1,t}.alphaScore));
    end
    f = f + getobjfnConstants(P, 1);
end