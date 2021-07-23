%ObjFunc
function [H,f] = getObjFunc(P, lambda)
    
    nAssets = length(P.initialAssetList);
    assets = P.initialAssetList;
    nVars = nAssets;
    %since cplex solvers deals with a vectorized obj function
    %it is necessary to expand di^2 and di*dj terms
    H = zeros(nVars);  %initialize coeff matrix:
                            %   H = [COVij]  nxn
    f = zeros(1, nVars); %initialzie linear coeff vec: 
                            %   f = [-lambda*alpha
                            %   - 2*vari*wbi - 2*(sum(COVij*wbj(i!=j)))]  1xn  
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t)) %add -2*(COVij*wbj) to fi and -2*(COVji*wbi) to fj 
            f(P.i_asset(t)) = f(P.i_asset(t)) - (2*P.covij(t)*assets{1,P.j_asset(t)}.benchWeight); 
            f(P.j_asset(t)) = f(P.j_asset(t)) - (2*P.covij(t)*assets{1,P.i_asset(t)}.benchWeight); 
            H(P.j_asset(t),P.i_asset(t)) = P.covij(t); %get lower triangular
        else %add -lambda*alphai - 2*vari*wbi to fi
            f(P.i_asset(t)) = f(P.i_asset(t)) - (lambda*assets{1,P.i_asset(t)}.alphaScore); 
            f(P.i_asset(t)) = f(P.i_asset(t)) - (2*P.covij(t)*assets{1,P.i_asset(t)}.benchWeight); 
        end
        H(P.i_asset(t),P.j_asset(t)) = P.covij(t); %get upper triangular and diagonal
    end
    %{
    for t=(nAssets+1):2*nAssets  %loop over sigmas and multiply them by a little value
        f = f + (1e-6)*x(t);
    end
    %}
    %f = f';
end



