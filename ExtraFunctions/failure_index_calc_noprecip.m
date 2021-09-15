function [FI] = failure_index_calc_noprecip(DEM, pw, ps, So)
% failure_index_calc takes six inputs:
% DEM --> digital elevation model as a topotoolbox GRIDobj
% precip --> Mean annual precipetation grid clipped to the same extent as
% the DEM in mm yr^-1 (function will convert to m yr^-1)
% T --> soil transmissivity in m^2 yr^-1 (a good value is 1e3)
% pw --> water density kg m^-3 (a good value is 1e3)
% ps --> sediment density kg m^-3
% So --> threshold slope (degrees)
% b --> cell size (m)
% The function generates one output, FI, which the failure index (basically
% the inverse of the factory of safety for a cohesionless material)

So = tand(So);

DEM = fillsinks(DEM);

slope = gradient8(DEM,'deg');
slope.Z(slope.Z > 70) = 45; % this is just put in to take of weird values from Katies Grid



FI = (tand(slope.Z) ./ So) .* (1./( 1 - (pw/ps)));

minThresh = 0.5;                                
maxThresh = 1.5;

FI = (FI./minThresh) - 1;
FI(FI < 0) = 0;
FI = FI./maxThresh;
FI(FI > 1) = 1;  
end
