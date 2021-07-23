function [A, b] = getLinCon_TE(P)%original
%function [A, lhs, rhs] = getLinCon_Classical(P)%original

%original
    %WITH CARD AND FLOOR-CEILING CONSTRAINTS
    nAssets = length(P.initialAssetList);
    cardSectorSet = length(P.sectorSet); 
    nVars = 3*nAssets;
    %{
    nIneqLinConstraints = (r6)(2*length(P.sectorSet)) + (r9.1) (r9.2&r9.3)2*nAssets;
    %}
    nLinConstr = 2 + (2*nAssets) + (2*length(P.sectorSet));
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    A=zeros(nLinConstr,nVars);
    b=zeros(1,nLinConstr);
    auxc = 1; %auxiliar index that loops through constraints 
    
    %r6.1: -sum(di) <= 0.20, where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
        b(auxc) = 0.05;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = -1;
            end
        end
        auxc = auxc + 1;
    end
    %r6.2: sum(di) <= 0.20, where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
        b(auxc) = 0.2;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = 1;
            end
        end
        auxc = auxc + 1;
    end
    
    % sigma constraints
    %9.1 20 <= sum(sigmai) <= 40
    for i = (2*nAssets)+1:nVars %loop through sigmas 
        A(auxc,i) = -1; %r9.1: -sum(Sigmai) <= -20
        A(auxc + 1,i) = 1; %r9.2: sum(Sigmai) <= 40)
    end
    b(auxc) = -5;
    b(auxc + 1) = 5;
    auxc = auxc + 2; %2 constraints added
    %boundary*Sigma constraints (without asset limit (r9.1)! This is not a NP-HARD problem anymore!)
    %9.3 wi - ub*sigma <= 0
    sub = 1; %w upper bound
    for i = nAssets+1:2*nAssets
        A(auxc,i) = 1;
        A(auxc,i + nAssets) = -sub;
        b(auxc) = 0;
        auxc = auxc + 1;
    end
end