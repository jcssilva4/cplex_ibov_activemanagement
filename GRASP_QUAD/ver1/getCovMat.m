function [COVIJ] = getCovMat(P, S, COV)
    assets = P.initialAssetList;
    if(nargin < 2)
        %get full COVIJ
        COVIJ = zeros(length(assets), length(assets));
        for t=1:length(P.covij)  %loop over covij vector
            if(P.i_asset(t) ~= P.j_asset(t))
                COVIJ(P.j_asset(t), P.i_asset(t)) = P.covij(t); %get lower triangular
            end
            COVIJ(P.i_asset(t), P.j_asset(t)) = P.covij(t); %get diagonal and upper triangular
        end
    else
        %extract a reduced COVIJ from the original COV matrix
        %S is the sorted(ascend dir) index list of selected stocks
        %COV is the full COV matrix
        %fprintf('get reduced COVIJ')
        COVIJ = zeros(length(S), length(S)); %initialization
        for i=1:length(S)  %loop over covij vector
            for j = 1:length(S)
                COVIJ(i, j) = COV(S(i),S(j)); 
            end
        end
    end
end