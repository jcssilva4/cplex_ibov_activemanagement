%return an asset's numeric integer associated with a SEDOL
function idx = getIdx(SEDOL, Assets)
    idx = 0; %initialize idx variable
    for i = 1: length(Assets)
        if(strcmp(SEDOL, Assets{1,i}.sedol))
            idx = i; %index identified
        end
    end
end