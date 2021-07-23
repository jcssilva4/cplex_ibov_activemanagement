function AS = getActiveShare(x, P, S) %returns 0 or 1

    % if [sum = sum(min(wi,wbi))] <= 0.4, return AS = 1
    sum = 0;
    min = 0;
    assets = P.initialAssetList;
    for i = 1:length(P.initialAssetList)
        if(ismember(i,S)) %if stock i is included in the Portfolio S
            if(x(find(ismember(S,i))) < assets{1,i}.benchWeight) %if wi < wbi (min = wi)
                min = x(find(ismember(S,i)));
            else %min = wbi
                min = assets{1,i}.benchWeight;
            end
        else %wi is not included in S, therefore wi = 0 is added without question
            min = 0;
        end
        sum = sum + min;
    end
    AS = sum;
end