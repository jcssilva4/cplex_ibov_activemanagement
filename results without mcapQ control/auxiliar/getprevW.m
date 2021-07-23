function [prevW] = getprevW ( lastR ,lastW)

somaAdj = 0;

lastAdjW = lastW *( 1 + lastR);

for i = 1:1:n
somaAdj = somaAdj +  lastW *( 1 + lastR);  
end

prevW = lastAdjW/somaAdjW;
end

%incomplete