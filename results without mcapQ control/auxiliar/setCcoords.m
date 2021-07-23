function [i, j, c] = setCcoords(y, m, currDay)
    %set "calendar coordinates"
    c = calendar(y,m);
    i = 1; %ith calendar row
    j = 1; %jth calendar column
    while (c(i,j) ~= currDay) 
        j = j + 1;
        if (j > 7)
            i = i + 1;
            j = 1;
        end
    end
end