function TE = getTE(xreduced, P, S) %get tracking error
    %using di as x
    %fprintf('Calculating Tracking Error')
    assets = P.initialAssetList;
    cov = 0;
    TE = 0;
    %get original di vector(S%P500) solved by GRASP+CPLEX
    x = zeros(1,length(P.initialAssetList)); %initialize
    for i = 1:length(P.initialAssetList)
        if(ismember(i,S)) %if stock i is included in the Portfolio S
          %  fprintf('\nlength(xred): %d requested idx:%d val: %f',...
           %     length(xreduced),find(ismember(S,i)), xreduced(find(ismember(S,i))));
            %{
            if(length(xreduced(find(ismember(S,i))))>1)
                S 
                ismember(S,i)
            end
        %}
            x(i) = xreduced(find(ismember(S,i)));
        end
    end
    %Calculate TE
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
        end
        di = x(P.i_asset(t));
        dj = x(P.j_asset(t));
        TE = TE + di*cov*dj;
    end
    TE = sqrt(TE);
    %{
%fprintf('Calculating Tracking Error')
    assets = P.initialAssetList;
    cov = 0;
    TE = 0;
    xreduced;
    %get original vector(S%P500) solved by GRASP+CPLEX
    x = zeros(1,length(P.initialAssetList)); %initialize
    for i = 1:length(P.initialAssetList)
        if(ismember(i,S)) %if stock i is included in the Portfolio S
          %  fprintf('\nlength(xred): %d requested idx:%d val: %f',...
           %     length(xreduced),find(ismember(S,i)), xreduced(find(ismember(S,i))));
            %{
            if(length(xreduced(find(ismember(S,i))))>1)
                S 
                ismember(S,i)
            end
        %}
            x(i) = xreduced(find(ismember(S,i)));
        end
    end
    %Calculate TE
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
        end
        di = x(P.i_asset(t))-assets{1,P.i_asset(t)}.benchWeight;
        dj = x(P.j_asset(t))-assets{1,P.j_asset(t)}.benchWeight;
        TE = TE + di*cov*dj;
    end
    TE = sqrt(TE);
%}
end