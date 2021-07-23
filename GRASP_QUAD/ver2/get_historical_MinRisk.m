function minRisk = get_historical_MinRisk(Portfolios, t0, tf)
    %get min risk from the last 3 months
    minRisk = 1e3; %arbitrary large number to begin
    for period = t0:tf
        var = 0;
        P = Portfolios{1,period};
        assets = P.initialAssetList; %get assetlist 
        for t=1:length(P.covij)  %loop over covij vector
            if(P.i_asset(t) ~= P.j_asset(t))
                cov = 2*P.covij(t); %sum covij and covji
            else
                cov = P.covij(t);
            end 
            xiBench = assets{1,P.i_asset(t)}.benchWeight;
            xjBench = assets{1,P.j_asset(t)}.benchWeight;
            var = var + xiBench*cov*xjBench;
        end
        risk = sqrt(var);
        if(risk < minRisk && risk > 0)
            minRisk = risk; %if current risk is lower than minRisk, then minRisk = risk
        end
    end
end