%%% Script to make a spatial plot of Foehn Wind Occurance in Dronning Maud
%%% Land using MARv2 data. This script with calculate the occurance of
%%% Foehn events based on given Foehn Conditions (from Banwell et al.,
%%% 2021) for each of the gridpoints of MAR, and will finally export the
%%% Foehn events map as tiff. 

%%% 26-May-2023, (c) anirudha.mahagaonkar@npolar.no 
clc; clear; close;

%% Data Location and details
% dataFolder = 'E:\MARv2\';
dataFolder = '\\nett.npolar\personlig\TromsoA2M\anirudha.mahagaonkar\DATA_REPOSITORY\MARv2';
fileWind = 'WIND.2014-2022.ANj.nc'; fpath_wind = fullfile(dataFolder, fileWind);
fileTemp = 'TEMP.2014-2022.ANj.nc';fpath_temp = fullfile(dataFolder, fileTemp);
fileRHum = 'HUM.2014-2022.ANj.nc';  fpath_RHum = fullfile(dataFolder, fileRHum);
fileLatLon = 'GEO.2014-2022.ANj.nc'; fpath_latlon = fullfile(dataFolder, fileLatLon);

% Save location
rawDataOutput = 'D:\NPI_Work_Research\MARv2\data\processed\';
outputFolder = 'D:\NPI_Work_Research\MARv2\outputs\Foehn\';

% Set Thresholds
WSTH = 3.5; % WindSpeed Threshold
TEMPTH = 1; % Temperature Threshold
RHTH = -5; % Relative Humidity Threshold
WindFrom = 0; % Wind Direction From
WindTo = 360; % wind Direction Till

% See information
% ncdisp(fpath_RHum); 

%% Read Information
% Read Time
disp('Reading Time from the NC file ...');
rawtime = ncread(fpath_temp, 'TIME'); % Days since 01-03-2010
time = rawtime + datenum(2010,03,01); 
ddatetime = datetime(time, 'ConvertFrom', 'datenum');

% Read X and Y from file
disp('Reading Lat/Lon from the NC file ...');
lat = ncread(fpath_latlon, 'LAT');
lon = ncread(fpath_latlon, 'LON');

% Read Wind Speed
disp('Reading Wind Speeds from the NC file ...');
% wSpeed = ncread(fpath_wind, 'UV');
% wSpeed = reshape(wSpeed(:,:,1,:), [500, 190, size(wSpeed, 4)]);
% save([rawDataOutput, 'wSpeed.mat'], 'wSpeed');
load([rawDataOutput, 'wSpeed.mat']);

% Read U and V for wind direction
disp('Reading UU Data from the NC file ...');
% UU = ncread(fpath_wind, 'UU');
% UU = reshape(UU(:,:,1,:), [500, 190, size(UU, 4)]);
% save([rawDataOutput, 'UU.mat'], 'UU');
load([rawDataOutput, 'UU.mat']);

disp('Reading VV Data from the NC file ...');
% VV = ncread(fpath_wind, 'VV');
% VV = reshape(VV(:,:,1,:), [500, 190, size(VV, 4)]);
% save([rawDataOutput, 'VV.mat'], 'VV');
load([rawDataOutput, 'VV.mat']);

% Read Temperature
disp('Reading Temperature from the NC file ...');
% rawtemp = ncread(fpath_temp, 'ST2');
% temp = reshape(rawtemp(:,:,1,:), [500, 190, size(rawtemp, 4)]);
% save([rawDataOutput, 'temp.mat'], 'temp');
load([rawDataOutput, 'temp.mat']);

% Read Relative Humidity
disp('Reading Relative Humidity from the NC file ...');
% relhum = ncread(fpath_RHum, 'RH');
% relhum = reshape(relhum(:,:,1,:), [500, 190, size(relhum, 4)]);
% save([rawDataOutput, 'relhum.mat'],  'relhum');
load([rawDataOutput, 'relhum.mat']);

% Calculate Wind Direction
dirRad = atan2(vv, uu);
wdir = mod(180 + (180/pi .* dirRad), 360);

%% Loop
% Loop to extract information from each required leeward and windward
% points. 
disp('Estimating Foehn activity point by point...');

% Pre-define a variable to accumulate days with Foehn Activity
foehnDays = zeros(size(lat));

% Pre-define a variable to make a daily foehn array.
dailyFoehnRecords = zeros([size(relhum)]);

for yyear = 2014:2021
    % disp(['Meltyear: ', num2str(yyear), '-', num2str(yyear+1)])
    startdate = datetime(yyear, 01, 01); % disp(['Start Date: ', datestr(startdate)]);
    enddate = datetime(yyear+1, 01, 01); % disp(['End Date: ', datestr(enddate-1)]);

    % Pre-define a variable to accumulate days with Foehn Activity
    % This variable only for Foehn in the selected time
    foehnDaysSelectedTime = zeros(size(lat));

    % Select Time
    DJFTime = ddatetime(ddatetime >= startdate & ddatetime < enddate);

    % Number of days in last year
    if yyear == 2014
        NumDays = 0;
    else
        NumDays = NumDays + daysact(datetime(yyear-1, 01, 01), datetime(yyear, 01, 01));
    end

    % Foehn Assessment pixel by pixel
    for rowi = 1:size(lat, 1)
        for colj = 1:size(lat, 2)
            % Display the point that is currently being processed
            disp(['Processing for point (row/col) :', num2str(rowi), '/', num2str(colj), ':', num2str(yyear)]);
            
            % That the point is now chosen, extract data for that point and
            % process year by year. 
            % Select data
            DJFptWind = wSpeed(rowi, colj, ddatetime >= startdate & ddatetime < enddate);
            DJFptTemp = temp(rowi, colj, ddatetime >= startdate & ddatetime < enddate);
            DJFptRH = relhum(rowi, colj, ddatetime >= startdate & ddatetime < enddate);
            DJFptWDir = wdir(rowi, colj, ddatetime >= startdate & ddatetime < enddate);
            
            % Assessment of Foehn Events
            % Start with T1
            tstep = 1;
    
            % Variable to count the number of events in a given season
            foehnEventCount = 0;
            
            % Loop through every timestep
            for i = 1:length(DJFptRH)
                
                % Check that the tstep value doesn't exceed the length of array
                if tstep == length(DJFptRH)
%                     disp('Ending the Loop');
                    break
                end
                
                % Check for Foehn conditions between tstep & tstep+1 
                if (DJFptTemp(tstep+1)-DJFptTemp(tstep) >= TEMPTH &&...
                        DJFptWind(tstep+1)-DJFptWind(tstep) >= WSTH &&...
                        DJFptRH(tstep+1)-DJFptRH(tstep) <= RHTH &&...
                        isAngBetween(DJFptWDir(tstep), WindFrom, WindTo))
        
                    % Since a Foehn event is found add one to event count
                    foehnEventCount = foehnEventCount + 1;

                    % If the condition is true => Foehn Activity
                    % current tstep => pre foehn condition.
                    % tstep+1 => first timestep with foehn condition
                    preFStart = tstep;
        
                    % From the pre condition, use the next tsteps to assess how
                    % long the Foehn Condition exists
                    foehnDayscount = 0;
        
                    searchstep = preFStart;
                    while (DJFptTemp(searchstep+1)-DJFptTemp(preFStart) >= TEMPTH &&...
                            DJFptWind(searchstep+1)-DJFptWind(preFStart) >= WSTH &&...
                            DJFptRH(searchstep+1)-DJFptRH(preFStart) <= RHTH &&...
                            isAngBetween(DJFptWDir(tstep), WindFrom, WindTo))
        
                        % Add one to the day count each time the condition
                        % is met.
                        foehnDayscount = foehnDayscount + 1;
%                         disp(['Foehn Day at TStep :', num2str(searchstep+1)]);
                        
                        % Assign '1' to the tstep (searchstep) in
                        % dailyFoehnRecords
                        dailyFoehnRecords(rowi,colj, NumDays+searchstep) = 1;
        
                        % Move to next tstep
                        searchstep = searchstep+1;
                        
                        % Break the loop if searchstep goes beyond the
                        % length of the selected timeperiod array
                        if searchstep == length(DJFptRH)
                            break
                        end
                    end

                    % Push the loop to jump the identified days and find new
                    % Foehn event
                    tstep = searchstep;

                    % Add this FoehnDaysCount value to predefined array.
                    foehnDays(rowi, colj) = foehnDays(rowi, colj) + foehnDayscount;
                    foehnDaysSelectedTime(rowi, colj) = foehnDaysSelectedTime(rowi, colj) + foehnDayscount;
                else
                    tstep = tstep + 1;
                end
            end
        end
    end
    % Save FoehnDays variable - for each loop step (e.g., Meltyear)
    saveName = (['FoehnDays_', datestr(startdate, 'mmmyyyy'), '-', datestr(enddate-1, 'mmmyyyy'), '.mat']);
%     save([outputFolder, saveName], "foehnDaysSelectedTime");
end

% Save FoehnDays variable - for each loop step (e.g., Meltyear)
saveName2 = (['MeanFoehnDays_', datestr(startdate, 'mmm'), '-', datestr(enddate-1, 'mmm'), '.mat']); 
% save([outputFolder, saveName2], "foehnDays");
% load('D:\NPI_Work_Research\MARv2\data\processed\dailyFoehnRecords_3pt5.mat');

%% Save the dailyFoehnRecords variable as .mat file. 
% This can be used for calculating number of Foehn Days for a chosen
% period. 
% save([rawDataOutput, 'dailyFoehnRecords_5ms.mat'], "dailyFoehnRecords", '-v7.3');

%% Convert to TIFF
% Define variables
cellsize = 5000;    % pixel size in m or distance between 2 grid points m

% Covert lat/lon to XY 3031
[X3031, Y3031] = geod2utm(lon, lat, 'spolar');

% Create refmat arrays
xPS = ceil(min(X3031(:))) : cellsize : floor(max(X3031(:))); %define projected grid
yPS = ceil(min(Y3031(:))) : cellsize : floor(max(Y3031(:)));

% Create MeshGrid
[X, Y] = meshgrid(xPS, yPS);

% Make RefMat
R = make_refmat(X, Y, cellsize, 'center');
R = double(R);

% Flip Var data
var1 = flip(foehnDays');

% Export as GTIFF
geotiffwrite(fullfile(outputFolder, [saveName2(1:end-4), '012014-122021_5.tif']), var1, R, 'CoordRefSysCode', 3031);

%% PCOLOR PLOT

shpFile1 = 'D:\NPI_Work_Research\DML Digitizing\DML_Coastline_wgs84.shp';
shpFile2 = 'D:\NPI_Work_Research\DML Digitizing\DML_Nunataks_wgs84.shp';

[dmlLon, dmlLat] = extractshpxy(shpFile1);
[NunLon, NunLat] = extractshpxy(shpFile2);

figure('visible','on');
set(gcf,'position',[200, 200, 1500, 550]);

% Set up projection
m_proj('stereographic','lat', -90,'lon', 0);

% Plot
p3 = m_pcolor(lon, lat, sum(dailyFoehnRecords(:,:,1:2921), 3, 'omitnan')); shading interp; 
p3.FaceAlpha = 0.70;

m_line(dmlLat, dmlLon, 'Color', [0 0 0], 'LineWidth', 1);
m_line(NunLat, NunLon, 'Color', [0 0 0], 'LineWidth', 1);

% Plot Extents
xlim([-0.1 0.22]); ylim([0.25 0.37]);
set(gca,'XTick',[], 'YTick',[]);

m_ruler([.35 .55], 0.1, 3,'fontsize',12, 'ticklen',.01);

numColors = 256;
positions = linspace(0, 1, numColors);
redC = [0.7, 0, 0]; lightRedC = [1, 0.5, 0.5]; whiteC = [1, 1, 1];
red2wht = interp1([0, 0.5, 1], [redC; lightRedC; whiteC], positions);
colormap(flipud(red2wht)); colorbar;

caxis([100 500])

%% Plotting Wind Rose for selected points

% ptx = 97; pty = 109;  name = 'RL';
% ptx = 417; pty = 106; name = 'RBE';
ptx = 118; pty = 128; name = 'Ekstrom';

figure
title(name);

ptSpeed = reshape(wSpeed(ptx, pty, :), 3165, []);
ptDir = reshape(0 - wdir(ptx, pty, :), 3165, []);
ptFoehn = reshape(dailyFoehnRecords(ptx, pty, :), 3165, []);

WRSpeed = ptSpeed(ptFoehn > 0);
WRDir = ptDir(ptFoehn > 0); 

% Define Plotting options for Wind Rose
Options = {'anglenorth', 0, 'angleeast', 90, ...
    'LegendType', 2, ...
    'labels', {'N (0°)', 'E (90°)', 'S (180°)', 'W (270°)'}, ...
    'freqlabelangle', 45, ...
    'MaxFrequency', 50};

subplot(121)
Options1 = [Options, {'axes', gca}];
WindRose(WRDir, WRSpeed, Options1);

subplot(122)
Options2 = [Options, {'axes', gca}];
WindRose(ptDir, ptSpeed, Options2);

% Add windrose plots on projected map
j = [97, 118, 417]; k = [109, 128, 106];
m_windrose(lon(97,118),lat(97,118),wdir(97,118),wspeed(97,118)*0.1, 'size', 1);
colormap(m_colmap('jet')); caxis([0 24]);
