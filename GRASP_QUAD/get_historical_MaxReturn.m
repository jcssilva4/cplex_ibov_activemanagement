function maxReturn = get_historical_MaxReturn(Portfolios, t0, tf)
    %get min return from the last 3 months
    maxReturn = -1e3; %arbitrary large negative number to begin
    for period = t0+1:tf
        ret = 0;
        P = Portfolios{1,period};
        assets = P.initialAssetList; %get assetlist 
        for i = 1:length(assets)
            ret = ret + assets{1,i}.benchWeight*assets{1,i}.r;
        end
        if(ret > maxReturn)
            maxReturn = ret; %if current risk is lower than minRisk, then minRisk = risk
        end
    end
end