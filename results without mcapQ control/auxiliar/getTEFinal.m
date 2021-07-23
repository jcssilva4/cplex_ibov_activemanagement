function TE = getTEFinal(x, P) %get tracking error
    %fprintf('Calculating Tracking Error')
    assets = P.initialAssetList;
    cov = 0;
    TE = 0;
    %Calculate TE
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
        end
        di = x(P.i_asset(t));
        dj = x(P.j_asset(t));
        TE = TE + di*cov*dj;
    end
    TE = sqrt(TE);
    %{
    assets = P.initialAssetList;
    cov = 0;
    TE = 0;
    %Calculate TE
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
        end
        di = x(P.i_asset(t))-assets{1,P.i_asset(t)}.benchWeight;
        dj = x(P.j_asset(t))-assets{1,P.j_asset(t)}.benchWeight;
        TE = TE + di*cov*dj;
    end
    TE = sqrt(TE);
    %}
end