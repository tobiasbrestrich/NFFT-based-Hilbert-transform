function [w, y] = getVoronoiWeightingAndRefPoints(x, scaling, plot, rad)
%GETVORONOIWEIGHTING Calculates area of Voronoi cells with special
%attention to border cells and returns vertices within 'rad'.
%   xy              The seeds around which the voronoi cells shall be generated
%   scaling         Must be >= 1. 
%                   For withConvexHull=1, decides how close the convex hull is
%                   to the outter most boundary point.
%                   For withConvexHull=0, 
                
%   plot            If 'on' then the vornoi cells will be ploted (if 'off' not)
%   rad             Radius, in which the vertices are accepted
%RETURN:
%   w               Individual Voronoi poylong's surface area
%   y               Coordinates of vertices within 'rad'





% Outer cells are set to zero
[vertices,cell] = voronoin(x);
if strcmp(plot,'on')
    figure(101);
    hold on
    voronoi(x(:,1), x(:,2));
end
w = zeros(length(cell),1);
for i = 1:length(cell)
    v1 = vertices(cell{i},1) ; 
    v2 = vertices(cell{i},2) ;
    w(i) = polyarea(v1,v2);
end
% Ensure no NaN or too big areas are replaced (e.g. with median/zero) 
w_median = median(w,'omitnan');
repl_indices = [find(w > scaling.*w_median); find(isnan(w))];
w(repl_indices) = 0;

%get y_fastsum
nv = length(vertices(:,1));
v = NaN(nv, 2);
for i = 1 : nv
    if norm(vertices(i,:)) < rad && ~any(isnan(vertices(i,:))) && ~any(isinf(vertices(i,:))) %ensure NaN, Inf and out of circle points are excluded
        v(i,:) = vertices(i,:);
    end
end
ri = ~isnan(v(:,1)) | ~isnan(v(:,2));
y = v(ri,:);

end