%ObjFunc
function f = getObjF_qcqp(P, lambda, S, COVIJ)
    assets = P.initialAssetList;
    %since cplex solvers deals with a vectorized obj function
    %it is necessary to expand di^2 and di*dj terms
    f = zeros(1, length(S)); %initialzie linear coeff vec: 
    for i=1:length(S)  
        for j = 1:length(S)
            if(i==j)  %add -lambda*wi*alphai 
                f(i) = f(i) - (lambda*assets{1,S(i)}.alphaScore); 
            end
            %add -covij*wi*wbj - covij*wj*wbi 
                f(i) = f(i) - COVIJ(i,j)*assets{1,S(j)}.benchWeight;
                f(j) = f(j) - COVIJ(i,j)*assets{1,S(i)}.benchWeight;
        end
    end
end



