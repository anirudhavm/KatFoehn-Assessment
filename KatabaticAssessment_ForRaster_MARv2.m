%%% Script to make a spatial plot of Katabatic Wind Occurance in Dronning Maud
%%% Land using MARv2 data. This script with calculate the occurance of
%%% Katabatic winds based on the wind direction and surface topography or
%%% the aspect, and will finally export the Katabatic events map as tiff. 

%%% 30-May-2023, (c) anirudha.mahagaonkar@npolar.no 
clc; clear; close;

%% Data Location and details
dataFolder = '\\nett.npolar\personlig\TromsoA2M\anirudha.mahagaonkar\DATA_REPOSITORY\MARv2';
fileWind = 'WIND.2014-2022.ANj.nc'; fpath_wind = fullfile(dataFolder, fileWind);
fileLatLon = 'GEO.2014-2022.ANj.nc'; fpath_latlon = fullfile(dataFolder, fileLatLon);

% Save location
rawDataOutput = 'D:\NPI_Work_Research\MARv2\data\processed\';
outputFolder = 'D:\NPI_Work_Research\MARv2\outputs\Katabatics\';

% Define Thresholds
aspectRange = 15;

% See information
ncdisp(fpath_wind); 
atmlay = ncread(fpath_wind, 'ATMLAY');

% ncdisp('S:\DATA_REPOSITORY\MARv2\ICE.2014-2022.ANj.nc');

%% Read Information to workspace
% Read Time
disp('Reading Time from the NC file ...');
rawtime = ncread(fpath_wind, 'TIME'); % Days since 01-03-2010
time = rawtime + datenum(2010,03,01); 
ddatetime = datetime(time, 'ConvertFrom', 'datenum');

% Read X and Y from file
disp('Reading Lat/Lon from the NC file ...');
lat = ncread(fpath_latlon, 'LAT');
lon = ncread(fpath_latlon, 'LON');

% Read Wind Speed
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

% Read Aspect information
disp('Reading surface aspect information...');
load([rawDataOutput, 'smoothAsp.mat']); aspect = smoothAsp;
% aspect = flip(imread('D:\NPI_Work_Research\REMA_DEM\Aspect10km_NearNeigh_MAREx.tif')',2);

%% For Saving to TIFF file at the end
% Covert lat/lon to XY 3031
[X3031, Y3031] = geod2utm(lon, lat, 'spolar');

% Create refmat arrays
cellsize = 5000;    % pixel size in m or distance between 2 grid points m
xPS = ceil(min(X3031(:))) : cellsize : floor(max(X3031(:))); %define projected grid
yPS = ceil(min(Y3031(:))) : cellsize : floor(max(Y3031(:)));

% Create MeshGrid
[X, Y] = meshgrid(xPS, yPS);

% Make RefMat
R = make_refmat(X, Y, cellsize, 'center');
R = double(R);

%% Pixel by pixel estimation
disp('Estimating Katabatic activity pixel by pixel...');

% Initiate variable to store katabatic occurances
meanKatabatics = zeros(size(aspect));
dailyKatRecords = zeros(size(uu));

% Loop
for yyear = 2014:2021
    % Define period of interest
    startdate = datetime(yyear, 01, 01); disp(['Start Date: ', datestr(startdate)]);
    enddate = datetime(yyear+1, 01, 01); disp(['End Date: ', datestr(enddate-1)]);
    
    % Select Time
    DJFTime = ddatetime(ddatetime >= startdate & ddatetime < enddate);
    
    % Initiate variable to store katabatic occurances
    selectedKatabatics = zeros(size(aspect));

    % Number of days in last year
    if yyear == 2014
        NumDays = 0;
    else
        NumDays = NumDays + daysact(datetime(yyear-1, 01, 01), datetime(yyear, 01, 01));
    end

    % Choose/Loop through pixels
    for rowi = 1:size(lat, 1)
        for colj = 1:size(lat, 2)
            
            % Display the point that is currently being processed
            disp(['Processing for point (row/col) :', num2str(rowi), '/', num2str(colj)]);

            % Select Aspect of the pixel/point
            pixelAspect = aspect(rowi,colj);

            % If aspect is -9999 = ignore the step
            if pixelAspect == -9999 
%                 disp('Skipping point as aspect is -9999');
                % Do nothing
            else
                % Select U and V Data
                DJFUU = reshape(uu(rowi,colj, ddatetime >= startdate & ddatetime < enddate), size(DJFTime));
                DJFVV = reshape(vv(rowi,colj, ddatetime >= startdate & ddatetime < enddate), size(DJFTime));
                
                % Calculate Wind Direction
                DJFdirRad = atan2(DJFVV, DJFUU);
                DJFwdir = mod(180 + (180/pi .* DJFdirRad), 360);
                
                % Calculate Wind Speed
                DJFWSpeed = sqrt(DJFUU.^2 + DJFVV.^2);

                % Define High and Low Range
                DirLow = mod(pixelAspect - (aspectRange+15), 360);
                DirHigh = mod(pixelAspect + aspectRange, 360);

                for dayNumber = 1:length(DJFwdir)

                    % Check condition - W Dir in line with surface aspect
                    if isAngBetween(DJFwdir(dayNumber), DirLow, DirHigh) && ...
                       DJFWSpeed(dayNumber) >= 0

                        % Add 1 as Katabatic Activity
                        meanKatabatics(rowi, colj) = meanKatabatics(rowi, colj) +1;
                        selectedKatabatics(rowi, colj) = selectedKatabatics(rowi, colj) +1;

                        % Mark respective day as Katabatic Day
                        dailyKatRecords(rowi, colj, NumDays+dayNumber) = 1;
                    end
                end
            end
        end
    end

    % Save the yearly meltyear (Selected period) plots 
%     var1 = flip(selectedKatabatics');
%     saveName1 = (['StrongKatDays_',  num2str(aspectRange), '_', datestr(startdate, 'mmmyyyy'), '-', datestr(enddate-1, 'mmmyyyy'), '.tif']); 
%     geotiffwrite(fullfile(outputFolder, saveName1), var1, R, 'CoordRefSysCode', 3031);
end

% Save mat variable
% save([rawDataOutput, '/dailyStrongKatDays_15Coriolis.mat'], "dailyKatRecords", "-v7.3")

%% Save as TIFF
% Flip Var data
var2 = flip(sum(dailyKatRecords, 3)');

% Export as GTIFF
saveName2 = (['KatDaysCoriolis_', num2str(aspectRange),'_', datestr(startdate, 'mmm'), '-', datestr(enddate-1, 'mmm'), '.tif']); 
geotiffwrite(fullfile(outputFolder, saveName2), var2, R, 'CoordRefSysCode', 3031);
