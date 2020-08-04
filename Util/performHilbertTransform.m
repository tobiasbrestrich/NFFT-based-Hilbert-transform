function [dX, dY] = performHilbertTransform(dZ, HX, HY)
%PERFORMHILBERTTRANSFORM Calculates the 3D-Hilbert transform in the 
%frequency domain
%   dZ          the component to be transformed
%   HX          the x-direction Hilbert transform operator in the frequency domian
%   HY          the y-direction Hilbert transform operator in the frequency domian
%RETURN:
%   dX          the x-component of the Hilbert transform solution
%   dY          the y-component of the Hilbert transform solution


spec_dZ = fft2(dZ);

dX = real(ifft2(HX .* spec_dZ));
dY = real(ifft2(HY .* spec_dZ));
