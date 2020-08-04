function [HX, HY] = provideHilbertOperator(KX, KY)
%PROVIDEHILBERTOPERATOR Calculates the Hilbert transform operators.
%   KX          Wavenumbers in x-direction
%   KY          Wavenumbers in y-direction
%RETURN:
%   HX          the x-direction Hilbert transform operator in the frequency domian
%   HY          the y-direction Hilbert transform operator in the frequency domian

KXY = sqrt(KX.^2 + KY.^2);

KXY(KXY == 0) = 1; %Inf;

HX = 1i * KX ./ KXY;
HY = 1i * KY ./ KXY;
