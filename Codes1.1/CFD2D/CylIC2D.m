function [rho,rhou,rhov,Ener] = CylIC2D(x, y, time);
  
%  function [rho,rhou,rhov,Ener] = CylIC2D(x, y, time)
%  Purpose: Impose uniform plane flow 

% Example is Mach ** 0.4 ** flow in wind tunnel
gamma = 1.4;

% Inflow conditions -- uniform inflow
rhoin = 1.4; uin = 0.4; vin = 0.0; pin = 1.0;
Ein = pin/(gamma-1.0) + 0.5*rhoin*(uin^2+vin^2);

% pack modified conserved variables

rho  = rhoin*ones(size(x));
rhou = rhoin*(1/.41)^2*6*(y+.2).*(0.41 - (y+.2));
rhov = 0;
Ener = Ein + 0.5*(rhou.^2 + rhov.^2)./rho;
return
