function [g] = hiTrafo1D1NonUniform(f,x,y,dx)
%HITRAFO1D Hilbert transfrom without singularity (changed var. approach)
%   Adjusted the Hilbert transform integral, to not have a singularity
%   f   coefficient of the Hilbert transform
%   x   start nodes of the transform (1D)
%   y   target nodes of the transform (1D)
%   dx  weighting of f
%RETURN:
%   g   solution of the Hilbert transform


nx = length(x);
ny = length(y);
g = zeros(ny,1);
x = x+200;

%zero-padding
f = [zeros(nx,1);f;zeros(nx,1)];

% convolution
for i = 1 : ny
    for j = 1 : nx
       g(i) = g(i) + (f(nx+i+j) - f(nx+i-j))/x(j)*dx(j); 
    end
end
g = g./pi;

end

