function [gammat,zref,pref,sigref] = gsw_gammat_analytic_theta_eos80(s,th)
%
% fun_gammat_analytic: Compute thermodynamic neutral density based on an
% analytical expression of Lorenz reference density
%
% INPUT:
%   
%   s          : practical salinity
%   th          : potential temperature (deg C)
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
% Coefficients for numerator of equation of state for density
%------------------------------------------------------------
a0 =             9.9984085444849347d+02 ;
a1 =             7.3471625860981584d+00 ;
a2 =            -5.3211231792841769d-02 ;
a3 =             3.6492439109814549d-04 ;
b0 =             2.5880571023991390d+00 ;
b1 =            -6.7168282786692355d-03 ;
b2 =             1.9203202055760151d-03 ;
c0 =             1.1798263740430364d-02 ;
c1 =             9.8920219266399117d-08 ;
c2 =             4.6996642771754730d-06 ;
c3 =             2.5862187075154352d-08 ;
c4 =             3.2921414007960662d-12 ;

% Coefficients for denominator of equation of state for density
%--------------------------------------------------------------
d0 =             1.0000000000000000d+00 ;
d1 =             7.2815210113327091d-03 ;
d2 =            -4.4787265461983921d-05 ;
d3 =             3.3851002965802430d-07 ;
d4 =             1.3651202389758572d-10 ;
e0 =             1.7632126669040377d-03 ;
e1 =             8.8066583251206474d-06 ;
e2 =             1.8832689434804897d-10 ;
e3 =             5.7463776745432097d-06 ;
e4 =             1.4716275472242334d-09 ;
f0 =             6.7103246285651894d-06 ;
f1 =             2.4461698007024582d-17 ;
f2 =             9.1534417604289062d-18 ;

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

% Define equation of state for density 
%-------------------------------------
anum = @(s,sqrts,th,th2,p,pth) a0 + th.*( a1 + th.*( a2 + th*a3)) + ...
        s.*( b0 + th*b1 + s*b2) + ...
        p.*( c0 + th2*c1 + s*c2 - p.*(c3 + th2*c4));
aden = @(s,sqrts,th,th2,p,pth) ...
        d0 + th.*( d1 + th.*( d2 + th.*( d3 + th*d4))) + ...
        s.*( e0 - th.*( e1 + th2*e2 ) + sqrts.*( e3 + th2*e4)) + ...
        p.*( f0 - pth.*(th2*f1 + p*f2)); 
rho = @(s,sqrts,th,th2,p,pth) anum(s,sqrts,th,th2,p,pth)./ ...
    aden(s,sqrts,th,th2,p,pth); 

% Compute the reference positions
% -------------------------------
zmin = 0.; zmax = 6000.; 
zref = s; ztop = s; zbot=s; 
%buoyancy = 0.*sr; 

ztop(:) = zmin;
zbot(:) = zmax;
zref(:) = 2000.;
%zref_new = zref; 

% Valid points
%-------------

th2 = th.*th; sqrts = sqrt(s); 

for n=1:30
    %disp(n)

    % Compute density at reference pressure pr(zref)
    p = pr(zref); pth = p.*th; 
    rho_ref = rho(s,sqrts,th,th2,p,pth); 

    % Compute buoyancy 
    buoyancy = rhor(zref) - rho_ref; 

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
sigref = rho(s,sqrts,th,th2,pref,pref.*th) - 1000.;
gammat = sigref - f(x); 

end

