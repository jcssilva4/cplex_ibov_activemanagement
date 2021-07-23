function AS = checkActiveShare(x, P, S) %returns 0 or 1

    % if [sum = sum(min(wi,wbi))] <= 0.4, return AS = 1
    sum = 0;
    assets = P.initialAssetList;
    for i = 1:length(P.initialAssetList)
        if(ismember(i,S)) %if stock i is included in the Portfolio S
            if(x(find(ismember(S,i))) < )
            x(i) = x(find(ismember(S,i)));
        else %wi is not included in S, therefore x = 0 , and wbi is added without question
           sum = sum + assets{1,i}.benchWeight;
        end
    end
end