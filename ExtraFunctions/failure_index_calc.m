function [FI] = failure_index_calc(DEM, precip, T, pw, ps, So, b)
% failure_index_calc takes six inputs:
% DEM --> digital elevation model as a topotoolbox GRIDobj
% precip --> Mean annual precipitation grid clipped to the same extent as
% the DEM in mm yr^-1 (function will convert to m yr^-1)
% T --> soil transmissivity in m^2 yr^-1 (a good value is 1e3)
% pw --> water density kg m^-3 (a good value is 1e3)
% ps --> sediment density kg m^-3
% So --> threshold slope (degrees)
% b --> cell size (m)
% The function generates one output, FI, which the failure index (basically
% the inverse of the factor of safety for a cohesionless material)

So = tand(So);

DEM = fillsinks(DEM);

Flow = FLOWobj(DEM);
DrainA = flowacc(Flow).*(Flow.cellsize^2);

slope = gradient8(DEM,'deg');
slope.Z(slope.Z > 70) = 45; % this is just put in to take off weird values from Katies Grid

precip_meters = precip;
precip_meters.Z = precip_meters.Z./1000;

W = (precip_meters.Z .* DrainA.Z) ./ (b * T .*sind(slope.Z));
W(W > 1) = 1;

FI = (tand(slope.Z) ./ So) .* (1./( 1 - W.*(pw/ps)));

minThresh = 0.5;                                %Can all this be folded into the failure_index_calc function?
maxThresh = 1.5;

FI = (FI./minThresh) - 1;
FI(FI < 0) = 0;
FI = FI./maxThresh;
FI(FI > 1) = 1;  
end
