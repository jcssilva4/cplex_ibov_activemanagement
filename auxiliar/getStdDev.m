function stdDev = getStdDev(w, P)
  %Calculate stdDev (NOT STD DEV, IT IS VAR!)
    stdDev = 0;
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
        end
        stdDev = stdDev + w(P.i_asset(t))*cov*w(P.j_asset(t));
    end
    %stdDev = sqrt(stdDev);
end