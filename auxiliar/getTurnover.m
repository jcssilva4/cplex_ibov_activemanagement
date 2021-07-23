function [turnover] = getTurnover( currentW , prevW )

somadifW = 0;

%calculo do modulo da dif entre o peso atual e o anterior
if currentW > prevW
    difW = currentW - prevW;
else
    difW = prevW - currentW;
end
somadifW = somadifW + difW;

turnover = somadifW;
end