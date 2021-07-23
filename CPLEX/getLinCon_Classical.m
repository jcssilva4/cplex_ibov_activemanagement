function [A, b] = getLinCon_Classical(P)%original
%function [A, lhs, rhs] = getLinCon_Classical(P)%original
%original
    %WITH CARD AND FLOOR-CEILING CONSTRAINTS
    nAssets = length(P.initialAssetList);
    cardSectorSet = length(P.sectorSet); 
    nVars = 2*nAssets;
    %{
    nIneqLinConstraints = (r6)(2*length(P.sectorSet)) + (r9.1) (r9.2&r9.3)2*nAssets;
    %}
    nLinConstr = 2 + (2*nAssets);% + (2*length(P.sectorSet));
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    A=zeros(nLinConstr,nVars);
    b=zeros(1,nLinConstr);
    auxc = 1; %auxiliar index th    at loops through constraints 
    auxw = 1; %auxiliar index that loops through weights vars
    %{
    %r6.1: -sum(xi) <= -0.05, where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
        b(auxc) = -0;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = -1;
            end
            
        end
        auxc = auxc + 1;
    end
    %r6.2: sum(xi) <= 0.15, where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
        b(auxc) = 0.5;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = 1;
            end
        end
        auxc = auxc + 1;
    end
    %}
    % sigma constraints
    %9.1 20 <= sum(sigmai) <= 40
    for i = nAssets+1:nVars %loop through sigmas 
        A(auxc,i) = -1; %r9.1: -sum(Sigmai) <= -20
        A(auxc + 1,i) = 1; %r9.2: sum(Sigmai) <= 40)
    end
    b(auxc) = -5;
    b(auxc + 1) = 5;
    auxc = auxc + 2; %2 constraints added
    %boundary*Sigma constraints (without asset limit (r9.1)! This is not a NP-HARD problem anymore!)
    %9.3 wi - ub*sigma <= 0
    sub = 0.2; %w upper bound
    for i = 1:nAssets
        A(auxc,i) = 1;
        A(auxc,i + nAssets) = -sub;
        b(auxc) = 0;
        auxc = auxc + 1;
    end
%{
 %WITH CARD AND FLOOR-CEILING CONSTRAINTS
    nAssets = length(P.initialAssetList);
    cardSectorSet = length(P.sectorSet); 
    nVars = 2*nAssets;
    %{
    nIneqLinConstraints = (r6)(2 + 2 + 2*length(P.sectorSet)) + (r9.1) (r9.2&r9.3)2*nAssets;
    %}
    nLinConstr = 2 + nAssets;%(2*nAssets); %+(length(P.sectorSet));
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    A=zeros(nLinConstr,nVars);
    rhs=zeros(1, nLinConstr);
    lhs=zeros(1, nLinConstr);
    auxc = 1; %auxiliar index that loops through constraints 
    auxw = 1; %auxiliar index that loops through weights vars
    
    %sum(w) = 1
    A(auxc,1:nAssets) = 1;  
    rhs(auxc) = 1;
    lhs(auxc) = 1;
    auxc = auxc + 1;
    %{
    %SECTOR CONSTRAINTS
    %r6.2: 0.05 <= sum(xi) <= 0.20, where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
        lhs(auxc) = 0;
        rhs(auxc) = 1;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = 1;
            end
        end
        auxc = auxc + 1;
    end
    %}
    % sigma constraints
    %9.1 20 <= sum(sigmai) <= 40
    for i = nAssets+1:nVars %loop through sigmas 
        A(auxc,i) = 1; 
    end
    lhs(auxc) = 10; %minimum number of assets
    rhs(auxc) = 10;%maximum number of assets
    auxc = auxc + 1; 
    
    %boundary*Sigma constraints (without asset limit (r9.1)! This is not a NP-HARD problem anymore!)
    %9.2-sub 0 <= -wi + ub*sigma <= ub-lb
    for i = 1:nAssets
        A(auxc,i) = -1;
        A(auxc,i + nAssets) = 1;
        lhs(auxc) = 0;
        rhs(auxc) = 1;
        auxc = auxc + 1;
    end
    %{
    %9.3 -(1-0.001) <= wi - ub*sigma <= 0
    %remake-> 0 <= -wi + lb*sigmai <= lb
    for i = 1:nAssets
        A(auxc,i) = -1;
        A(auxc,i + nAssets) = sub;
        lhs(auxc) = 0;
        rhs(auxc) = sub;
        auxc = auxc + 1;
    end
   %}
%}
end
    
