%%% Script to indentify the temporal distribution of Katabatic winds
%%% for a given point. 

%%% 08-Feb-2024, (c) anirudha.mahagaonkar@npolar.no 
clc; clear; close;

%% Data Location and details
dataFolder = '\\nett.npolar\personlig\TromsoA2M\anirudha.mahagaonkar\DATA_REPOSITORY\MARv2\';
fileWind = 'WIND.2014-2022.ANj.nc'; fpath_wind = fullfile(dataFolder, fileWind);
fileLatLon = 'GEO.2014-2022.ANj.nc'; fpath_latlon = fullfile(dataFolder, fileLatLon);

% Save location
rawDataOutput = 'D:\NPI_Work_Research\MARv2\data\processed\';
outputFolder = 'D:\NPI_Work_Research\MARv2\outputs\Katabatics\';

% Define Thresholds
aspectRange = 15;

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

KatMonthlyAll = [];

% Loop to extract information from each required leeward and windward
% points. 
for caseNo = 6%1:5
    % Assign indices of points
    if caseNo == 1; area_key = 'RBW'; areaName = 'Roi Baudouin West';
        disp('Processing for RB West');
        rowi = 409; colj = 84;
%         figure;
    elseif caseNo == 2; area_key = 'RBC'; areaName = 'Roi Baudouin Center';
        disp('Processing for RB Center');
        rowi = 381; colj = 88;
%         figure;
    elseif caseNo == 3; area_key = 'RBE'; areaName = 'Roi Baudouin East';
        disp('Processing for RB East');
        rowi = 361; colj = 89;
%         figure;
    elseif caseNo == 4; area_key = 'FB'; areaName = 'Fimbulisen';
        disp('Processing for NV');
        rowi = 204; colj = 135;
%         figure;
    elseif caseNo == 5; area_key = 'NV'; areaName = 'Nivlisen';
        disp('Processing for NV');
        rowi = 274; colj = 129;
%         figure;
    elseif caseNo == 6; area_key = 'JO'; areaName = 'Jotneisen';
        disp('Processing for JO');
        rowi = 340; colj = 101;
%         figure;
    end

    % Initiate variable to store katabatic occurances
    katDaysAcc = [];

    % Loop
    for yyear = 2014:2021
        % Define period of interest
        startdate = datetime(yyear, 01, 01); disp(['Start Date: ', datestr(startdate)]);
        enddate = datetime(yyear+1, 01, 01); disp(['End Date: ', datestr(enddate-1)]);

        % Select Time
        DJFTime = ddatetime(ddatetime >= startdate & ddatetime < enddate);
        
        % Select U and V Data
        DJFUU = reshape(uu(rowi,colj, ddatetime >= startdate & ddatetime < enddate), size(DJFTime));
        DJFVV = reshape(vv(rowi,colj, ddatetime >= startdate & ddatetime < enddate), size(DJFTime));

        % Calculate Wind Direction
        DJFdirRad = atan2(DJFVV, DJFUU);
        DJFwdir = mod(180 + (180/pi .* DJFdirRad), 360);

        % Calculate Wind Speed
        DJFWSpeed = sqrt(DJFUU.^2 + DJFVV.^2);

        % Initiate variable to store katabatic occurances
        KatDays = false(size(DJFTime));

        % Select Aspect of the pixel/point
        pixelAspect = aspect(rowi,colj);

        % If aspect is -9999 = ignore the step
        if pixelAspect == -9999
            %                 disp('Skipping point as aspect is -9999');
            % Do nothing
        else
            % Define High and Low Range
            DirLow = mod(pixelAspect - (aspectRange+15), 360);
            DirHigh = mod(pixelAspect + aspectRange, 360);

            for dayNumber = 1:length(DJFwdir)

                % Check condition - W Dir in line with surface aspect
                if isAngBetween(DJFwdir(dayNumber), DirLow, DirHigh) && ...
                        DJFWSpeed(dayNumber) >= 0

                    % Mark respective day as Katabatic Day
                    KatDays(dayNumber) = 1;
                end
            end
        end

        % Populate yearly KatDays into KatDaysAcc variable
        % If leap year, eliminate value of Feb29. Feb29 = Index Pos 60.
        if length(KatDays) == 366
            KatDays(60) = [];
        end
        katDaysAcc = [katDaysAcc, KatDays];

    end
        
    % Calculate MeanKatDaysAcc
    mean_katDaysAcc = mean(katDaysAcc, 2, 'omitnan');
    
    % Make Timetable for mothly group summary
    KatAccTT = timetable([datetime(2014,01,01):datetime(2014,12,31)]', mean_katDaysAcc);

    % Group monthly
    katAccMonthly = retime(KatAccTT, "monthly", "sum");

    % Plot Figure
    figure
    bar(1:12, katAccMonthly.mean_katDaysAcc); 
%     ylim([0 15])
    title(areaName)
    grid on; box on;

    % To calculate average of all points, populate Mean KatDaysAcc into one
    % variable 
    KatMonthlyAll = [KatMonthlyAll, katAccMonthly.mean_katDaysAcc];
end

figure
bar(1:12, mean(KatMonthlyAll, 2, 'omitnan'), 0.9, 'FaceColor', 'k' ,'FaceAlpha', 1);
box on; grid on;
xticklabels({'J', 'F', 'M', 'A', 'M', 'J',...
             'J', 'A', 'S', 'O', 'N', 'D'})
xlabel('Month'); ylabel('Days with Katabatic Occurence');

%% Plot Foehn and Kat together
% run the full script of foehn wind assessment for points before running
% this following code. 

figure; hold on;
yyaxis right;
bar(1:12, mean(FoehnMonthlyAll, 2, 'omitnan'), 0.9, 'FaceColor', 'k' ,'FaceAlpha', 1);
ylim([0 10]);
ylabel('Days with Foehn Occurences');

yyaxis left;
bar(1:12, mean(KatMonthlyAll, 2, 'omitnan'), 0.9, 'FaceColor', 'k', 'FaceAlpha', 0.5);
grid on; box on;
xticks(1:12);
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',...
             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
xlabel('Month'); ylabel('Days with Katabatic Occurences');
legend('Katabatic', 'Foehn');




