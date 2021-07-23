function const = getobjfnConstants(P, lambda)
    const = 0;
    assets = P.initialAssetList;
    %fval = foptfound + const; %because you expanded di^2 and didj
    %const = (wbi*COVij*wbj)->bench var + lambda*wbi*alphai->bench mean return
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
            const = const +  lambda*assets{1,P.i_asset(t)}.benchWeight*assets{1,P.i_asset(t)}.alphaScore;
        end
        const = const + assets{1,P.i_asset(t)}.benchWeight*cov*assets{1,P.j_asset(t)}.benchWeight;
    end
end