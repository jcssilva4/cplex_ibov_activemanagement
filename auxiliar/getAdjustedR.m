function [Rtxadj] = getAdjustedR ( Ropt, turnover)

Rtxadj = Ropt - 0.005*turnover;

end