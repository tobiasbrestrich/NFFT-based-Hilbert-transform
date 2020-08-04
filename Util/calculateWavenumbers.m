function [KX, KY, KXY] = calculateWavenumbers(nx, ny, dx, dy)
%CALCULATEWAVENUMBERS Computes wavenumbers (freq. domain) for an equispaced
%grid according data given in the space domain.
%   nx       number of points in x-direction
%   nx       number of points in y-direction
%   dx       spacing of points in x-direction
%   dx       spacing of points in y-direction
%RETURN:
%   KX       wavenumber for x coordinates
%   KY       wavenumber for y coordinates
%   KXY      normed wave numbers

[KX, KY] = ndgrid( ...
    [0:ceil(nx/2)-1, -floor(nx/2):-1], ...
    [0:ceil(ny/2)-1, -floor(ny/2):-1]);

KX = KX / floor(nx/2) * pi / dx;
KY = KY / floor(ny/2) * pi / dy;

KXY = sqrt(KX.^2 + KY.^2);
