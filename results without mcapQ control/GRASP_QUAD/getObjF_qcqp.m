%ObjFunc
function f = getObjF_qcqp(P, lambda, S, COVIJ, nVars)
    nVars = length(S);
    assets = P.initialAssetList;
    %since cplex solvers deals with a vectorized obj function
    %it is necessary to expand di^2 and di*dj terms
    f = zeros(1, nVars); %initialzie linear coeff vec: 
                            %   f = [-lambda*alpha
                            %   - 2*vari*wbi - 2*(sum(COVij*wbj(i!=j)))]  1xn  
    for i=1:length(S)  
        for j = 1:length(S)
            if(i~=j) 
            else %add -lambda*alphai - 2*vari*wbi to fi
                f(i) = f(i) - (lambda*assets{1,S(i)}.alphaScore); 
            end %add -2*(COVij*wbj) to fi
            f(i) = f(i) - (2*COVIJ(i, j)*assets{1,S(j)}.benchWeight);
        end
    end
end



