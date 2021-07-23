function [x0] = getx0(P) 
    nAssets = length(P.initialAssetList);
    x0 = zeros(1,2*nAssets);
    benchw = zeros(1,nAssets);
    %get benchw
    benchwsum = 0;
    for asset = 1:nAssets
        benchw(asset) = P.initialAssetList{1,asset}.benchWeight;
        benchwsum = benchwsum +benchw(asset);
    end
    %to generate our population, we just need to redistribute points of
    %assets that are zero (w < 0.001) randomly. 
    excess = 0; 
    for asset = 1:nAssets
        if(benchw(asset) < 0.001) 
            excess = excess + benchw(asset); %save asset i's contribution (to sum the rest of the w's up to 1)
            benchw(asset) = 0; %now asset i is not part of the benchmark anymore
        end
    end
    %start random distribution for each individual
    extraw = 0;
    currentExcess = 0; %current excess of each individual
    lastAssetIdx = 1;
    newsum = 0;
    perc = 0;
    currentExcess = excess;
    newsum = 0;
    for asset = 1:(nAssets-1)
        %check which assets are included in the benchmark and add extra
        %weight
        if(benchw(asset) > 0)
            perc = rand; %generate random number contained in (0,1)
            extraw = (perc)*(currentExcess); %get a percentage of currentExcess
            x0(asset) = benchw(asset) + extraw; %add this percentage to asset i weight
            currentExcess = currentExcess - extraw; %reduce currentExcess
            x0(nAssets + asset) = 1; %if wbench_i >= 0, then sigma_i = 1       
            lastAssetIdx = asset; %memorize the last asset idx included in the benchmark
        end
        newsum = newsum + x0(asset);
    end
    x0(lastAssetIdx) = x0(lastAssetIdx) + currentExcess; %the rest goes to the last asset
    newsum = newsum +  x0(lastAssetIdx);
    fprintf('w sum of x0 = %f\n',newsum);
end