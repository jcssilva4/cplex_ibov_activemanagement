function numStocksperQ = getnumStocksperQ(P, S)
    %count assets in each quintile
    numStocksperQ = zeros(1, 5);
    for q = 1:5
        for i = 1:length(S)
            if(ismember(S(i),P.idxQ{q}))
                numStocksperQ(q) = numStocksperQ(q) + 1;
            end
        end
    end
end