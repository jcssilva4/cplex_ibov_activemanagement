function [ropt_t,rbench_t] = get_ropt_t(P,Pbench)
    ropt_t = 0;
    rbench_t = 0;
    assets = P.initialAssetList;
    for i = 1:length(assets)
        ropt_t = ropt_t + assets{1,i}.weight*assets{1,i}.r;
    end
    assets = Pbench.initialAssetList;
    for i = 1:length(assets)
        rbench_t = rbench_t + assets{1,i}.weight*assets{1,i}.r;
    end
end