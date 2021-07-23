function [A, b] = getLinCon(P)
    fprintf('Go to getLinCon')
    nAssets = length(P.initialAssetList);
    cardSectorSet = length(P.sectorSet); 
    %{
    nIneqLinConstraints = (r5)(2*nAssets) + (r6)(2*length(P.sectorSet)) 
        + (r7)10 + (r8)2;
    %}
    nLinConstr = 12 + (2*nAssets) + (2*length(P.sectorSet));
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    A=zeros(nLinConstr,nAssets);
    b=zeros(1,nLinConstr);
    auxc = 1; %auxiliar index that loops through constraints 
    auxw = 1; %auxiliar index that loops through weights vars
    
    %r5.1: -Wi <= 0.05 - Wb
    for i = 1:nAssets
        A(auxc,i) = -1;
        b(auxc) = 0.05 - P.initialAssetList{1,i}.benchWeight;
        auxc = auxc + 1;
    end
    
    %r5.2: Wi <= 0.05 + Wb
    for i = 1:nAssets
        A(auxc,i) = 1;
        b(auxc) = 0.05 + P.initialAssetList{1,i}.benchWeight;
        auxc = auxc + 1;
    end
    %r6.1: -sum(Wi) <= 0.1 - sum(Wb), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = -1;
                b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
            end
        end
    auxc = auxc + 1;
    end
    %r6.2: sum(Wi) <= 0.1 + sum(Wb), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = 1;
                b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
            end
        end
    auxc = auxc + 1;
    end
    
    %r7.1: -sum(Wi) <= 0.1 - sum(Wb), where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                A(auxc,i) = -1;
                b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
            end
        end
        auxc = auxc + 1;
    end
    %r7.2: sum(Wi) <= 0.1 + sum(Wb), where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                A(auxc,i) = 1;
                b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
            end
        end
        auxc = auxc + 1;    
    end
    
    %r8
    b(auxc) = 0.1
    b(auxc + 1) = 0.1
    for i = 1:nAssets %loop through weights 
        A(auxc,i) = -P.initialAssetList{1,i}.beta; %r8.1: -sum(Wi*Betai) <= 0.1 - sum(Wb*Betai)
        b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
        A(auxc + 1,i) = P.initialAssetList{1,i}.beta; %r8.2: sum(Wi*Betai) <= 0.1 + sum(Wb*Betai)
        b(auxc + 1) = b(auxc + 1) + P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
    end
%}
%{
%WITH CARD AND FLOOR-CEILING CONSTRAINTS
nAssets = length(P.initialAssetList);
    cardSectorSet = length(P.sectorSet); 
    nVars = 2*nAssets;
    %{
    nIneqLinConstraints = (r5)(2*nAssets) + (r6)(2*length(P.sectorSet)) 
        + (r7)10 + (r8)2 + (r9.1) (r9.2&r9.3)2*nAssets;
    %}
    nLinConstr = 14 + (4*nAssets) + (2*length(P.sectorSet));
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    A=zeros(nLinConstr,nVars);
    b=zeros(1,nLinConstr);
    auxc = 1; %auxiliar index that loops through constraints 
    auxw = 1; %auxiliar index that loops through weights vars
    
    %r5.1: -Wi <= 0.05 - Wb
    for i = 1:nAssets
        A(auxc,i) = -1;
        b(auxc) = 0.05 - P.initialAssetList{1,i}.benchWeight;
        auxc = auxc + 1;
    end
    
    %r5.2: Wi <= 0.05 + Wb
    for i = 1:nAssets
        A(auxc,i) = 1;
        b(auxc) = 0.05 + P.initialAssetList{1,i}.benchWeight;
        auxc = auxc + 1;
    end
    %r6.1: -sum(Wi) <= 0.1 - sum(Wb), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = -1;
                b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
            end
        end
    auxc = auxc + 1;
    end
    %r6.2: sum(Wi) <= 0.1 + sum(Wb), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                A(auxc,i) = 1;
                b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
            end
        end
    auxc = auxc + 1;
    end
    
    %r7.1: -sum(Wi) <= 0.1 - sum(Wb), where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                A(auxc,i) = -1;
                b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
            end
        end
        auxc = auxc + 1;
    end
    %r7.2: sum(Wi) <= 0.1 + sum(Wb), where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:nAssets %loop over weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                A(auxc,i) = 1;
                b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
            end
        end
        auxc = auxc + 1;    
    end
    
    %r8
    b(auxc) = 0.1;
    b(auxc + 1) = 0.1;
    for i = 1:nAssets %loop through weights 
        A(auxc,i) = -P.initialAssetList{1,i}.beta; %r8.1: -sum(Wi*Betai) <= 0.1 - sum(Wb*Betai)
        b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
        A(auxc + 1,i) = P.initialAssetList{1,i}.beta; %r8.2: sum(Wi*Betai) <= 0.1 + sum(Wb*Betai)
        b(auxc + 1) = b(auxc + 1) + P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
    end
    auxc = auxc + 2; %2 constraints added
    % sigma constraints
    %9.1 50 <= sum(sigmai) <= 70
    for i = nAssets+1:nVars %loop through sigmas 
        A(auxc,i) = -1; %r9.1: -sum(Sigmai) <= -50
        A(auxc + 1,i) = 1; %r9.2: sum(Sigmai) <= 70)
    end
    b(auxc) = -50;
    b(auxc + 1) = 70;
    auxc = auxc + 2; %2 constraints added
    %boundary*Sigma constraints (without asset limit (r9.1)! This is not a NP-HARD problem anymore!)
    %9.2 -wi + lb*sigma <= 0
    slb = 0.001; %w lower bound
    for i = 1:nAssets
        A(auxc,i) = -1;
        A(auxc,i + nAssets) = slb;
        b(auxc) = 0;
        auxc = auxc + 1;
    end
    %9.3 wi - ub*sigma <= 0
    sub = 1; %w upper bound
    for i = 1:nAssets
        A(auxc,i) = 1;
        A(auxc,i + nAssets) = -sub;
        b(auxc) = 0;
        auxc = auxc + 1;
    end
%}
end