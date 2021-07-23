function [A,b] = get_qplincon(P)
%using di as x
    nVars = length(P.initialAssetList);
    cardSectorSet = length(P.sectorSet); 
    %{
    nIneqLinConstraints = (r6)(2*length(P.sectorSet)) + (r7)10 + (r8)2;
    %}
    nLinConstr = 2*nVars;
    %Linear Inequality Constraints (Ax <= b)
    %i -> constraint index, j-> element index(j= 1 = x(1))
    A=zeros(nLinConstr,nVars);
    b=zeros(1,nLinConstr);
    %-wi <= 0
    auxc = 1; %auxiliar index that loops through constraints 
    for i = 1:nVars
        A(auxc,i+nVars) = -1;
        b(auxc) = 0;
        auxc = auxc+1;
    end
    %wi <= 1
    for i = 1:nVars
        A(auxc,i+nVars) = 1;
        b(auxc) = 1;
        auxc = auxc+1;
    end