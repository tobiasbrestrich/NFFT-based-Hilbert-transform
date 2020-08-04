function g = attrHorCylinder(r, r0, a, drho)
%ATTRHORCYLINDER Calculates attraction of a infinite horizontal cylinder
%   r          Profile point, where attraction is to be calculated
%   r0         Center of the cylinder
%   a          Radius of the cylinder 
%   drho       Difference of density of the cylinder to the surrounding
%              rock
%RETURN:
%   g          Attraction in [mgal] at r  

    assert(length(r0)==3,'r0 must be 3-vector');
    assert(length(a)==1,'a must be scalar');
    assert(length(drho)==1,'drho must be scalar');
    
    gravKonst = 6.673*10^(-11);    
    si2mgal = 1e5;
    
    dr=norm(r-r0);
    m=drho*a^2*pi;
    
    g=-2*gravKonst*m*(r-r0)/dr^2;
    
    g=g*si2mgal;
end    