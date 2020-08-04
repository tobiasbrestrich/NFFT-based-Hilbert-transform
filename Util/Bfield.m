function B = Bfield(m, xp, xq)
%BFIELD Calculates the magnetic induction of a magnetized sphere
%   m           magnetic moment of sphere
%   xp          cooridnates of magnetic field
%   xq          coordinates of center of sphere
%RETURN:
%   B           magnetic induction at xq. Contains X,Y and Z components

mu = pi * 4e-7;
m = m(:);
x = xp(:) - xq(:);
R = norm(x);

B = mu * (3 * dot(x, m) * x - m * R^2) / (4 * pi * R^5);
