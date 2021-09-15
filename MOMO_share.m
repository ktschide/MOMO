% close all
% clear all
addpath(genpath('../topotoolbox-master'));
addpath(genpath('.'));

%in general:    XXX capitals for grids and matrices
%               Xxx for vectors
%               xxx for variables
rng('shuffle');
filename = 'MOMO_testrun';
max_ls_size = 7;            %power of max landslide size for random sampling (ex. 7 = 10^7)
model_time = 100000;        %model time [years]

%% INPUT BACKGROUND EROSION - Model works in cm in z direction, m in x & y, grams, years, atoms
backgE = 0.001;                                      %[cm/yr]
totalE = 0.02;                                       %[cm/yr]
lslide_E = totalE - backgE;                          %[cm/yr]

%% LOAD DATA

%RASTER RESOLUTION
cellsize_m = 30;                                    %[m]
cellsize_m2 = (cellsize_m)^2;                       %[m^2]

%COSMOGENIC PRODUCTION GRIDS
load('Production_BK_hindi.mat');             %Contains neutron, muon production grids as well as atmospheric pressure. Generated with Calc_Production function
PN = Production_BK_hindi.Pn;                 %Neutron production grid [at/g/yr]
PMU = Production_BK_hindi.Pmu;               %Muon production grid [at/g/yr]
PPRES = Production_BK_hindi.pres;            %Atmospheric pressure grid [hPa]

%PREDEFINED Leffs GRIDS and make sure production and precipitation grids are within limits
load Leffs_KS.mat
load consts_LSD
pp = Leffs_KS.pp;                           %pp -> pressure [hPa]
ee = 0.1 * Leffs_KS.ee;                     %ee -> erosion rate [mm converted to cm]
Leff10 = Leffs_KS.Leff10;                   %Leff10 -> e-folding length 10Be muons [g/cm^2]
minee = min(min(ee));
backgE(backgE < minee) =  minee;
maxee = max(max(ee));
backgE(backgE > maxee) = maxee;
minp = min(min(pp));
PPRES(PPRES < minp) =  minp;
maxp = max(max(pp));
PPRES(PPRES > maxp) = maxp;
XGRID = PPRES;
XGRID(~isnan(XGRID)) = backgE;
LEFF = interp2(pp,ee,Leff10,PPRES,XGRID);   %Make LEFF grid (see Balco, 2017) [g/cm^2]

%DEM for stability analysis [m.asl]
DEM = load('hindi_DEMc.mat');%                         DEM = load('Trisuli_resortDEM.mat');
DEM = DEM.DEMc;        %                               DEM = DEM.Trisulu_resortDEM;
DEMr = resample(DEM, cellsize_m);   %                  DEMr = resample(DEM, cellsize_m);

%PRECIPITATION for stability analysis [mm/yr]
PRECIP = load('hindi_precip.mat');
PRECIP = PRECIP.DEMprecip;
PRECIPr = resample(PRECIP,cellsize_m);

%% resize grids to same extent and size
% Recast PN, PMU and PPRES as GRIDobj
PN2 = DEMr;
PN2.Z = PN;
PN2.size = size(PN);
PN = PN2;

PMU2 = DEMr;
PMU2.Z = PMU;
PMU2.size = size(PMU);
PMU = PMU2;

PPRES2 = DEMr;
PPRES2.Z = PPRES;
PPRES2.size = size(PPRES);
PPRES = PPRES2;

clear('PN2', 'PMU2', 'PPRES2');

% make sure grids are the same size
[xMin(1,1),xMax(1,1),yMin(1,1),yMax(1,1)] = findCorners(PN);
[xMin(2,1),xMax(2,1),yMin(2,1),yMax(2,1)] = findCorners(PMU);
[xMin(3,1),xMax(3,1),yMin(3,1),yMax(3,1)] = findCorners(PPRES);
[xMin(4,1),xMax(4,1),yMin(4,1),yMax(4,1)] = findCorners(DEMr);
[xMin(5,1),xMax(5,1),yMin(5,1),yMax(5,1)] = findCorners(PRECIPr);

% resizing grids to largest overlapping area
xMinP = max(xMin);
xMaxP = min(xMax);
yMinP = max(yMin);
yMaxP = min(yMax);

PN = gridReSize(PN,xMinP,xMaxP,yMinP,yMaxP);
PMU = gridReSize(PMU,xMinP,xMaxP,yMinP,yMaxP);
PPRES = gridReSize(PPRES,xMinP,xMaxP,yMinP,yMaxP);
DEMr = gridReSize(DEMr,xMinP,xMaxP,yMinP,yMaxP);
PRECIPr = gridReSize(PRECIPr,xMinP,xMaxP,yMinP,yMaxP);

% clear unneeded values
clear('xMinP','xMaxP','yMinP','yMaxP','xMin','xMax','yMin','yMax');

% define refmat
DEMr.refmat = PN.refmat;
PRECIPr.refmat = PN.refmat;

PN = PN.Z;
PMU = PMU.Z;
PPRES = PPRES.Z;

%% MODEL PARAMETERS
sa = nansum(nansum(PPRES./PPRES))*cellsize_m2;    %Drainage area [m^2]

%COSMOGENIC NUCLIDES
att=160;                                          %attenuantion of CRN production curve [g/cm^2]
rho = 2.7;                                        %density of substrate/bedrock [g/cm^3]
lambda = (4.62*(10^-7));                          %CRN decay rate [1/yr]

%LANDSLIDES
amax = 1e7;                                       %max landslide area [m^2]
amin = cellsize_m2;                               %minimum landslide area [m^2]

%Larsen's Volume-Area scaling
gamV = 1.36;                                      %Also note that Larsen uses the convension gamma for the exponent (here gamV)
alphV = 10^-0.59; %
% for reference Himalaya soil gamV = 1.25, alphV = 10^-0.44; mixed gamV =
% 1.36, alphV = -0.59; BR gamV = 1.34, alphV = 10^-0.49

% Frequency-Area scaling
k = 39817;                                      %exp(4.68) Roback fig 3 %K = 39817; %fig 3 Roback et al 2017 (Log10p = (-2.47±0.11)Log10A + (4.68±0.54))
beta = 1.5;
p = @(a) k.*a.^(-beta);

%MODEL TIME STEPS
                                %total time period the model runs [years]
dt = 1;                                           %time spacing [years]
Tarray = 0:dt:model_time;                         %time vector [years]

[row,col] = size(PN);

%% OUTPUT VECTORS
Lsmassout = zeros(1,length(Tarray));  %total mass of sediment exported from domain at each time step from landslide erosion [g]
Mout = zeros(1,length(Tarray));
Aout = zeros(1,length(Tarray));

%% build up initial concentration array - steady state erosion concentration for E = background erosion
SCAPE_PN = PN ./ (lambda + (rho*backgE/att));
SCAPE_PMU = PMU ./ (lambda + ((rho*backgE)./LEFF));


%% set up Failure Index
tf = 1000;                                      %soil transmissivity [m^2/yr]
pw = 1000;                                      %water density [kg/m3]
ps = 2700;                                      %sediment density [kg/m3]
So = 35;                                        %threshold slope [degrees]

[FI] = failure_index_calc(DEMr, PRECIPr, tf, pw, ps, So, cellsize_m);
%[FI] = failure_index_calc_noprecip(DEMr, pw, ps, So);      %if not including precipitation data

report = model_time/100;

%% to record all erosion
ALL_SLIDES = zeros(size(PN));
bigslides = zeros(1,length(Tarray));

%%subbasins
smallbasins = load('trib_basins.mat');
smallbasins = smallbasins.trib_basins;
smallbasins = smallbasins.Z;

%nested catchments
trunkbasins = load('trunkbasins_hindi.mat');
trunkbasins = trunkbasins.trunkbasins;


%% TIME LOOP FOR MODELING OF EACH TIME STEP
for ttime=1:length(Tarray)
    
    DEPTH_SCAPE = zeros(size(PN)) + (backgE * dt);   %erosion depth during current time step [cm]
    DEPTH_SCAPE(isnan(PN)) = nan;
    
    landslidecount=0;
    time=ttime*dt;
    
    if ttime>report
        percent_done = 100*(1-(model_time-ttime)/model_time)       %tracking progress
        report = report+(model_time/100);
    end
    
    % update surface concentrations for decay and production
    SCAPE_PN = (SCAPE_PN + (PN * dt)) .* exp(-lambda*dt);                  %nuclide loss to radioactive decay and gain by production during dt - neutrons [at/g]
    SCAPE_PMU = (SCAPE_PMU + (PMU * dt)) .* exp(-lambda*dt);               %nuclide loss to radioactive decay and gain by production during dt - muons [at/g]
    
    %% Landslide decision making
    [Sizes] = random_LSareas(p,3,max_ls_size,9000);                     %generate list of landslide areas based on powerlaw expression [m^2]
    accepted = rand(size(FI)) < FI;                                     %Randomly choose cells more likely to slide
    
    total_Ls = 0;
    ls_year = 0;
    
    while ls_year < lslide_E
        
        landslidecount=landslidecount+1;
        la = Sizes(landslidecount);             %choose landslide area [m^2]
        
        %% Do rejection sampling
        %pick landslide location
        k = 0;
        while k < 1
            xloc=(round(rand*(col-1)))+1;
            yloc=(round(rand*(row-1)))+1;
            
            if ~isnan(PN(yloc,xloc)) && (accepted(yloc,xloc)) == 1     %check that landslide is within the grid and meets FI threshold
                disturbed_cells=la/(cellsize_m2);
                sqcells=round(sqrt(disturbed_cells));
                ydir = sqcells-1;
                xdir = sqcells-1;
                if (ydir+yloc)>row
                    ydir=row-yloc;
                end
                
                if (xdir+xloc)>col
                    xdir=col-xloc;
                end
                 FItest = FI(yloc:yloc+ydir,xloc:xloc+xdir);
                 if nanmean(nanmean(FItest)) > 0.5
                     k = k +1;
                 end
            end
        end
        
        Elevations(ttime,landslidecount) = DEMr.Z(yloc,xloc);
        
        [Locyarray, Locxarray, la] = fitlandslide(la, cellsize_m2, PN, xloc,yloc);                        %make sure landslide fits bounds of basin outline
        
        if la > 10^5
            bigslides(ttime) = la;
        end
        
        %% Use Larsen scaling to get depth
        volume = (alphV*la^gamV);                                               %[m^3]
        depth=(volume/la)*100;                                                  %landslide depth [cm]
        
        lstemp = Lsmassout(1,ttime)+(rho*depth*la*10000);
        ls_year = lstemp/(rho*sa*10000*dt);
        
        if ls_year > lslide_E && rem(ttime,2) == 0
            % the nothing
        else
            Lsmassout(1,ttime)=Lsmassout(1,ttime)+(rho*depth*la*10000);                   %track mass removed from landsliding [g]
            
            ls_year = (Lsmassout(1,ttime)/(rho*sa*10000*dt));                             %measure atoms and mass eroded by landslides in this year [cm]
            
            DEPTH_SCAPE(Locyarray,Locxarray) = DEPTH_SCAPE(Locyarray,Locxarray) + depth;   %add landsliding depth to background erosion [cm]
            
            Locxarray=[];
            Locyarray=[];
            
        end
    end
    
    MASS = DEPTH_SCAPE.* rho .* cellsize_m2 .* 10000;                                                   %total eroded mass during time step [g]
    AOUT_PN = cellsize_m2 .* 10000 .* SCAPE_PN .* att .* (1 - exp(-DEPTH_SCAPE .* rho ./ att));         %calculate neutrons removed during time step [atoms]
    AOUT_PMU = cellsize_m2 .* 10000 .* SCAPE_PMU .* LEFF .* (1 - exp(-DEPTH_SCAPE .* rho ./ LEFF));     %calculate muons removed during time step [atoms]
    
    %%WHOLE BASIN
    Mout(1,ttime) = nansum(nansum(MASS));                                             %mass of sediments exported by background erosion and landsliding [g]
    Aout(1,ttime) = (nansum(nansum(AOUT_PN)) + (nansum(nansum(AOUT_PMU))));           %total nuclide atoms exported (neutrons + muons) [atoms]
    
    %%subbasins
    fn = fieldnames(trunkbasins);
    for basin=1:numel(fn)
        trunkstats.depth(basin,ttime) = nansum(nansum(DEPTH_SCAPE(trunkbasins.(fn{basin}) == 1)));
        trunkstats.Aout(basin,ttime) = nansum(nansum(AOUT_PN(trunkbasins.(fn{basin}) == 1))) + nansum(nansum(AOUT_PMU(trunkbasins.(fn{basin}) == 1)));
        trunkstats.Mout(basin,ttime) = nansum(nansum(MASS(trunkbasins.(fn{basin}) == 1)));  
    end
    
    for sbasin=1:max(max(smallbasins))
        sbasin_stats.depth(sbasin,ttime) = nansum(nansum(DEPTH_SCAPE(smallbasins == sbasin)));
        sbasin_stats.Aout(sbasin,ttime) = nansum(nansum(AOUT_PN(smallbasins == sbasin))) + nansum(nansum(AOUT_PMU(smallbasins == sbasin)));
        sbasin_stats.Mout(sbasin,ttime) = nansum(nansum(MASS(smallbasins == sbasin)));  
    end
    
    SCAPE_PN = SCAPE_PN .* exp(-DEPTH_SCAPE .* rho ./ att);                      %updated surface concentration after removal of background and landsliding erosion
    SCAPE_PMU = SCAPE_PMU .* exp(-DEPTH_SCAPE .* rho ./ LEFF);                   %updated surface concentration after removal of background and landsliding erosion
    
    ALL_SLIDES = ALL_SLIDES + DEPTH_SCAPE;
end

tmass= sum(Mout);                           %total mass eroded [g]
tatoms = sum(Aout);                         %total atoms eroded [atoms]
Etotal=(tmass/(rho*sa*10000*time))*10;      %reports total long term erosion rate [mm]
C = Aout./Mout(1,:);                        %average concentration over time [atoms/g]

M = (movmean(Aout,500))./(movmean(Mout,500));  %moving window average

%%save workspace

save(filename, 'trunkstats','sbasin_stats','ALL_SLIDES', 'Etotal', 'Tarray', 'Mout', 'C', 'M', 'Aout', 'sa','DEMr','FI','SCAPE_PN','SCAPE_PMU','max_ls_size','backgE','totalE')
