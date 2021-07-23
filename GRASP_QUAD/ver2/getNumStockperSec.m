function numStockperSec = getNumStockperSec(P, S)
    numStockperSec = zeros(1, length(P.sectorSet));
    %count assets in each sector
    for s = 1:length(P.sectorSet) %loop over all sectors
        for i = 1:length(S) %loop over all assets
            if(ismember(S(i),P.idxSec{1,s}))
                numStockperSec(s) = numStockperSec(s) + 1;
            end
        end
    end
end