%ObjFunc
function f = getObjF_classical(P, lambda, S, COVIJ)  
    assets = P.initialAssetList;
    f = zeros(1, length(S)); %initialzie linear coeff vec: 
    for i=1:length(S)  
        f(i) = -(lambda*assets{1,S(i)}.alphaScore); 
    end
end



