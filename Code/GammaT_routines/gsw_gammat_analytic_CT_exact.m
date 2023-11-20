function [gammat,zref,pref,sigref] = gsw_gammat_analytic_CT_exact(sr,ct)
%
% fun_gammat_analytic: Compute thermodynamic neutral density based on an
% analytical expression of Lorenz reference density
%
% INPUT:
%   
%   sr          : reference composition salinity (g/kg) 
%   ct          : Conservative Temperature (deg C)
%
% OUTPUT: 
%   zref        : Reference position
%   pref        : Reference pressure
%   sigref      : Reference density 
%   gammat      : Thermodynamic neutral density 
%
% DEPENDENCIES: 
%   gsw_rho_CT_exact(sr,ct,p): equation of state for seawater as a function
%   of reference composition salinity (sr), conservative temperature (ct)
%   and pressure (p, in dbars) 
%
% AUTHOR: Remi Tailleux, University of Reading, 8 July 2020
%==========================================================================
% Set values of coefficients
% --------------------------
a = 4.56016575;
b = -1.24898501;
c =  0.00439778209;
d =  1030.99373;
e =  8.32218903; 

% Set polynomial correction

%      Linear model Poly8:
%      f(x) = p1*x^8 + p2*x^7 + p3*x^6 + p4*x^5 + 
%                     p5*x^4 + p6*x^3 + p7*x^2 + p8*x + p9
%        where x is normalized by mean 1440 and std 1470
%      Coefficients (with 95% confidence bounds):
p1 =   0.0007824;  % (0.0007792, 0.0007855)
p2 =   -0.008056;  % (-0.008082, -0.008031)
p3 =     0.03216;  % (0.03209, 0.03223)
p4 =    -0.06387;  % (-0.06393, -0.06381)
p5 =     0.06807;  % (0.06799, 0.06816)
p6 =    -0.03696;  % (-0.03706, -0.03687)
p7 =    -0.08414;  % (-0.08419, -0.0841)
p8 =       6.677;  % (6.677, 6.677)
p9 =       6.431;  % (6.431, 6.431)

% Set value of gravity 
% --------------------
grav = 9.81; 

% Define the different analytical functions
% -----------------------------------------
% drhordz = @(z) a.*(z+e).^b + c;
rhor = @(z) a.*(z+e).^(b+1)./(b+1) + c.*z + d;
pr = @(z) grav.*(a.*(z+e).^(b+2)./((b+1)*(b+2)) + c.*z.^2/2 + ...
    d.*z - a.*e.^(b+2)./((b+1)*(b+2))  )./1d4; 

% Polynomial correction
%----------------------
f = @(x) p9 + x.*( p8 + x.*( p7 + x.*( p6 + x.*( p5 ...
    + x.*( p4 + x.*( p3 + x.*( p2 + x.*p1)))))));

% Compute the reference positions
% -------------------------------
zmin = 0.; zmax = 6000.; 
zref = sr; ztop = sr; zbot=sr; 
%buoyancy = 0.*sr; 

ztop(:) = zmin;
zbot(:) = zmax;
zref(:) = 2000.;
%zref_new = zref; 

% Valid points
%-------------

for n=1:30
    %disp(n)
    % Compute buoyancy 
    buoyancy = rhor(zref) - gsw_rho_CT_exact(sr,ct,pr(zref)); 
    % Compute sign of buoyancy
    ss = sign(buoyancy);
    % Redefine zref depending on sign of buoyancy
    zref_new = 0.25.*(1.+ss).*(ztop+zref) + 0.25*(1.-ss).*(zbot+zref);
    ztop = 0.5*(1+ss).*ztop + 0.5*(1-ss).*zref;
    zbot = 0.5*(1-ss).*zbot + 0.5*(1+ss).*zref;
    zref = zref_new;
    %buo_max = max(abs(buoyancy(:)));
    %dz_max = max(zbot(:)-ztop(:)); 
    %disp([buo_max,dz_max])
end

% Compute analytic gammat
% -----------------------
pmean = 1440;
pstd = 1470; 

pref = pr(zref); 
x = (pref-pmean)./pstd;
sigref = gsw_rho_CT_exact(sr,ct,pref) - 1000.;
gammat = sigref - f(x); 

end

