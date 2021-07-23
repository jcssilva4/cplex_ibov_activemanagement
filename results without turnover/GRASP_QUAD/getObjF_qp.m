%ObjFunc
function f = getObjF_qcqp(P, lambda, S, COVIJ)
    %{
    assets = P.initialAssetList;
    %since cplex qcqp solver deals with a vectorized obj function
    %it is necessary to expand di^2 and di*dj terms
    f = zeros(1, length(S)); %initialzie linear coeff vec: 
    for i=1:length(S)  
        for j = 1:length(S)
            if(i==j)  %add -lambda*wi*alphai 
                f(i) = f(i) - (lambda*assets{1,S(i)}.alphaScore); 
            end
            %add -covij*wi*wbj - covij*wj*wbi 
                f(i) = f(i) - 2*COVIJ(i,j)*assets{1,S(j)}.benchWeight; %because we multiplied H by 2
                %f(j) = f(j) - COVIJ(i,j)*assets{1,S(i)}.benchWeight;
        end
    end
    %}
    
    %using di as a variable
    assets = P.initialAssetList;
    %since cplex qcqp solver deals with a vectorized obj function
    %it is necessary to expand di^2 and di*dj terms
    f = zeros(1, 2*length(S)); %initialzie linear coeff vec: 
    for i=1:length(S)  
        f(i) = -(lambda*assets{1,S(i)}.alphaScore); 
    end
end



