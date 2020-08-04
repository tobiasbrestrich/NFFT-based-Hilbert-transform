function [weighting] = getWeighting(randomNodes)
%GETWEIGHTING Calculates weighted 1D-Areas around the random source nodes
%   randomNodes     1D coordinates of the random source nodes
%RETURN:
%   weighting       weightings of source nodes

    N = length(randomNodes);
    weighting = zeros(N,1);
    
    %calculate firs and last area
    weighting(1) = (randomNodes(2)-randomNodes(1));
    weighting(N) = (randomNodes(N) - randomNodes(N-1));
    
    %calculate areas
    for i = 2 : N-1
        weighting(i) = (randomNodes(i+1)-randomNodes(i-1))/2;
    end

end

