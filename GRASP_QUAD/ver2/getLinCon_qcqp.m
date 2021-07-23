function [A, b] = getLinCon_qcqp(Portfolios, S, t, numOfPastPeriods)
    %using di as x
    P = Portfolios{1,t};
    nVars = length(S);
    cardSectorSet = length(P.sectorSet); 
    %{
    nIneqLinConstraints = (r6)(2*length(P.sectorSet)) + (r7)10 + (r8)2;
    %}
    nLinConstr = 12  + (2*length(P.sectorSet)) + nVars + 1;
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    A=zeros(nLinConstr,2*nVars);
    b=zeros(1,nLinConstr);
    auxc = 1; %auxiliar index that loops through constraints 
    auxw = 1; %auxiliar index that loops through weights vars
    
    %TRASFORM r5.1 and r5.2 in lb and ub, respectively
    %-di <= -0.001 + wbi
    for i = 1:nVars
        A(auxc,i) = -1;
        b(auxc) = -0.001 + P.initialAssetList{1, S(i)}.benchWeight;
        auxc = auxc+1;
    end
    
    %r6.1: -sum(di) <= 0.1 - (sum(bi) if i is not contained in S), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over active weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                if(ismember(i,S)) %if i belongs to S 
                    %{
                        ismember(i, S) returns a logical value 
                        ismember(S, i) returns a logical vector 
                    %}
                    %fprintf('\n%d equals to %d', i, S(find(ismember(S,i))))
                    A(auxc,find(ismember(S,i))) = -1; %search for i index in S
                else
                    b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
                end 
            end
        end
    auxc = auxc + 1;
    end
    %r6.2: sum(di) <= 0.1 + (sum(bi) if i is not contained in S), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over active weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                if(ismember(i,S))%if i belongs to S
                    A(auxc,find(ismember(S,i))) = 1;
                else
                    b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
                end 
            end
        end
    auxc = auxc + 1;
    end
    
    %r7.1: -sum(di) <= 0.1, where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over active weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                if(ismember(i,S))%if i belongs to S
                    A(auxc,find(ismember(S,i))) = -1;
                else
                    b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
                end 
            end
        end
        auxc = auxc + 1;
    end
    %r7.2: sum(di) <= 0.1 + sum(Wb), where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over active weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                if(ismember(i,S))%if i belongs to S
                    A(auxc,find(ismember(S,i))) = 1;
                else
                    b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
                end 
            end
        end
        auxc = auxc + 1;    
    end
    %r8
    b(auxc) = 0.1;
    b(auxc + 1) = 0.1;
    for i = 1:length(P.initialAssetList) %loop through active weights 
        if(ismember(i,S))%if i belongs to S
            A(auxc,find(ismember(S,i))) = -P.initialAssetList{1,i}.beta; %r8.1: -sum(di*Betai) <= 0.1 
            A(auxc + 1,find(ismember(S,i))) = P.initialAssetList{1,i}.beta; %r8.2: sum(di*Betai) <= 0.1 
        else
            b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
            b(auxc + 1) = b(auxc + 1) + P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
        end
    end
    auxc = auxc + 2;
    %min return constraint: -r_port <= -r_bench 
    b(auxc) = -get_historical_MaxReturn(Portfolios,...
        t - numOfPastPeriods, t - 1);
    scaleFactor = 1e4;
    for i = 1:length(P.initialAssetList) %loop through weights 
        expectedReturn = scaleFactor*P.initialAssetList{1,i}.alphaScore;
        A(auxc,find(ismember(S,i))+length(S)) = -expectedReturn;
    end

%{
    %original formulation
    nVars = length(S);
    cardSectorSet = length(P.sectorSet); 
    %{
    nIneqLinConstraints = (r6)(2*length(P.sectorSet)) + (r7)10 + (r8)2;
    %}
    nLinConstr = 12  + (2*length(P.sectorSet));
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    %A=zeros(nLinConstr,nVars);
    %b=zeros(1,nLinConstr);
    auxc = 1; %auxiliar index that loops through constraints 
    auxw = 1; %auxiliar index that loops through weights vars
    %{
    %TRASFORM r5.1 and r5.2 in lb and ub, respectively
    %r5.1: -Wi <= 0.05 - Wb  (but wi>=0.001) because it is included in the portfolio (GRASP SOLUTION)
    inclusion_lb = 0.001;
    for i = 1:nVars
        A(auxc,i) = -1;
        candidatelb = 0.05 - P.initialAssetList{1,S(i)}.benchWeight;  
        if(candidatelb < inclusion_lb), candidatelb = inclusion_lb;  end %if candidatelb is less than inc_lb, then we use inc_lb!!
        b(auxc) = candidatelb;
        auxc = auxc + 1;
    end
    
    %r5.2: Wi <= 0.05 + Wb That's the upper bound
    for i = 1:nVars
        A(auxc,i) = 1;
        b(auxc) = 0.05 + P.initialAssetList{1,S(i)}.benchWeight;
        auxc = auxc + 1;
    end
    %}
    %should i add the other remaining benchmarks?
    %r6.1: -sum(Wi) <= 0.1 - sum(Wb), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over all N weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                if(ismember(i,S)) %if i belongs to S 
                    %{
                        ismember(i, S) returns a logical value 
                        ismember(S, i) returns a logical vector 
                    %}
                    %fprintf('\n%d equals to %d', i, S(find(ismember(S,i))))
                    A(auxc,find(ismember(S,i))) = -1; %search for i index in S
                end 
                b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
            end
        end
    auxc = auxc + 1;
    end
    %r6.2: sum(Wi) <= 0.1 + sum(Wb), where i belongs to sector j
    for j = 1:cardSectorSet %loop over sectors
    b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over weights
            isector = P.initialAssetList{1,i}.sector;
            if(strcmp(isector,P.sectorSet(j))) %i belongs to sector j?
                if(ismember(i,S))%if i belongs to S
                    A(auxc,find(ismember(S,i))) = 1;
                end 
                b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
            end
        end
    auxc = auxc + 1;
    end
    
    %r7.1: -sum(Wi) <= 0.1 - sum(Wb), where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                if(ismember(i,S))%if i belongs to S
                    A(auxc,find(ismember(S,i))) = -1;
                end 
                b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight;
            end
        end
        auxc = auxc + 1;
    end
    %r7.2: sum(Wi) <= 0.1 + sum(Wb), where i belongs to mcapQ k
    for k  = 1:5 %loop over mcapQ
        b(auxc) = 0.1;
        for i = 1:length(P.initialAssetList) %loop over weights
            imcapQ = P.initialAssetList{1,i}.mCapQ;
            if(imcapQ == k) %i belongs to mcapQ k?
                if(ismember(i,S))%if i belongs to S
                    A(auxc,find(ismember(S,i))) = 1;
                end 
                b(auxc) = b(auxc) + P.initialAssetList{1,i}.benchWeight;
            end
        end
        auxc = auxc + 1;    
    end
    %r8
    b(auxc) = 0.1;
    b(auxc + 1) = 0.1;
    for i = 1:length(P.initialAssetList) %loop through weights 
        if(ismember(i,S))%if i belongs to S
            A(auxc,find(ismember(S,i))) = -P.initialAssetList{1,i}.beta; %r8.1: -sum(Wi*Betai) <= 0.1 - sum(Wb*Betai)
            A(auxc + 1,find(ismember(S,i))) = P.initialAssetList{1,i}.beta; %r8.2: sum(Wi*Betai) <= 0.1 + sum(Wb*Betai)
        end
        b(auxc) = b(auxc) - P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
        b(auxc + 1) = b(auxc + 1) + P.initialAssetList{1,i}.benchWeight*P.initialAssetList{1,i}.beta;
    end
%}
end