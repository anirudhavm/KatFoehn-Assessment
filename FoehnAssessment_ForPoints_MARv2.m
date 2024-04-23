%%% Script to make temporal plots of Wind Speed, Temperature and Relative
%%% Humidity for a chosen *particular point* from MARv2 gridpoints. This
%%% script is for assessing/identifying the Foehn wind activity using the
%%% Foehn condition defined by Banwell et al., 2021

%%% This script is also used for estimating the warming caused by adjacent
%%% ice rises by measuring the windward-to-leeside temp difference. This is
%%% shown through a box plot. 

%%% 22-May-2023, (c) anirudha.mahagaonkar@npolar.no 
clc; clear; close;

% Data Location and details
dataFolder = '\\nett.npolar\personlig\TromsoA2M\anirudha.mahagaonkar\DATA_REPOSITORY\MARv2\';
% dataFolder = 'E:\MARv2\';
fileWind = 'WIND.2014-2022.ANj.nc'; fpath_wind = fullfile(dataFolder, fileWind);
fileTemp = 'TEMP.2014-2022.ANj.nc';fpath_temp = fullfile(dataFolder, fileTemp);
fileRHum = 'HUM.2014-2022.ANj.nc';  fpath_RHum = fullfile(dataFolder, fileRHum);
fileLatLon = 'GEO.2014-2022.ANj.nc'; fpath_latlon = fullfile(dataFolder, fileLatLon);

% Save location
rawDataOutput = 'D:\NPI_Work_Research\MARv2\data\processed\';
outputFolder = 'D:\NPI_Work_Research\MARv2\outputs\';

% See information
% ncdisp(fpath_RHum); 

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
load([rawDataOutput, 'wSpeed.mat'])

% Read Temperature
disp('Reading Temperature from the NC file ...');
% rawtemp = ncread(fpath_temp, 'ST2');
% temp = reshape(rawtemp(:,:,1,:), [500, 190, size(rawtemp, 4)]);
% save([rawDataOutput, 'temp.mat'], 'temp');
load([rawDataOutput, 'temp.mat'])

% Read Relative Humidity
disp('Reading Relative Humidity from the NC file ...');
% relhum = ncread(fpath_RHum, 'RH');
% relhum = reshape(relhum(:,:,1,:), [500, 190, size(relhum, 4)]);
% save([rawDataOutput, 'relhum.mat'],  'relhum');
load([rawDataOutput, 'relhum.mat'])

% Read Specific Humidity
disp('Reading Specific Humidity from the NC file ...');
% spchum = ncread(fpath_RHum, 'QQ');
% spchum = reshape(spchum(:,:,1,:), [500, 190, size(spchum, 4)]);
% save([rawDataOutput, 'spchum.mat'],  'spchum');
load([rawDataOutput, 'spchum.mat'])

% Set Thresholds
WSTH = 3.5; % WindSpeed Threshold
TEMPTH = 1; % Temperature Threshold
RHTH = -5; % Relative Humidity Threshold

% Store Vars
maxTemp = zeros(8,5);
meanTemp = zeros(8,5);

% For BoxPlot
allFoehnDays = [];
allTempDiffStats = [];
allTempDiff = [];
allptLeeTemp = [];
allSummerTag = [];
grp = [];

% Loop to extract information from each required leeward and windward
% points. 
disp('Processing one meltyear at a time...');
for caseNo = 1:6
    % Assign indices of points
    if caseNo == 1; area_key = 'SO'; areaName = 'Riiser Larsen';
        disp('Processing for Søråsen Point');
        Llon_pt_row = 97; Llat_pt_col = 109;
        Wlon_pt_row = 110; Wlat_pt_col = 120;
        figure;
    elseif caseNo == 2; area_key = 'HR'; areaName = ' Ekstrom';
        disp('Processing for Halvfarryggen Point');
        Llon_pt_row = 118; Llat_pt_col = 128;
        Wlon_pt_row = 130; Wlat_pt_col = 132;
        figure;
    elseif caseNo == 3; area_key = 'NV'; areaName = 'Nivlisen';
        disp('Processing for Nivlisen point');
        Llon_pt_row = 261; Llat_pt_col = 138;
        Wlon_pt_row = 276; Wlat_pt_col = 138;
        figure;
    elseif caseNo == 4; area_key = 'MN'; areaName = 'Muninisen';
        disp('Processing for Muninisen point');
        Llon_pt_row = 325; Llat_pt_col = 116;
        Wlon_pt_row = 337; Wlat_pt_col = 114;
        figure;
    elseif caseNo == 5; area_key = 'RL'; areaName = 'Roi Baudouin';
        disp('Processing for Riiser Larsen Halvøya Point');
        Llon_pt_row = 413; Llat_pt_col = 94; % New correct point
        Wlon_pt_row = 433; Wlat_pt_col = 88;
        figure;
    elseif caseNo == 6; area_ket = 'REF'; areaName = 'Reference';
        disp('Processing for Reference Point');
        Llon_pt_row = 161; Llat_pt_col = 146;
        Wlon_pt_row = 176; Wlat_pt_col = 146;
    end
    
    T = tiledlayout(8,1,"TileSpacing","tight","Padding","tight");
    title(T, ['Foehn Activity on ', areaName]);
    countn = 1;
    
    ptTempDiff = [];
    ptLeeTemp = [];
    ptSummerTag = [];
    ptFoehnDays = [];

    for yyear = 2014:2022-1
        disp(['Meltyear: ', num2str(yyear), '-', num2str(yyear+1)])
        startdate = datetime(yyear, 01, 01); % disp(['Start Date: ', datestr(startdate)]);
        enddate = datetime(yyear+1, 01, 01); % disp(['End Date: ', datestr(enddate-1)]);
        
        DJFTime = ddatetime(ddatetime >= startdate & ddatetime < enddate);
        DJFWSpeed = wSpeed(:,:, ddatetime >= startdate & ddatetime < enddate);
        DJFTemp = temp(:,:, ddatetime >= startdate & ddatetime < enddate);
        DJFrh = relhum(:,:, ddatetime >= startdate & ddatetime < enddate);
        DJFsh = spchum(:,:, ddatetime >= startdate & ddatetime < enddate);

        % Check lat lon of selected points. 
%         disp('Manually check these Lat/Lon points to be correct...')
%         disp([num2str(lat(Llon_pt_row, Llat_pt_col)), ':', num2str(lon(Llon_pt_row, Llat_pt_col))]);
        
        % Extract values for point locations
%         disp('Extracting Wind, Temp and RH values for given point locations...')
        LptWind = reshape(DJFWSpeed(Llon_pt_row, Llat_pt_col, :), [size(DJFTemp, 3),1]);
        LptTemp = reshape(DJFTemp(Llon_pt_row, Llat_pt_col, :), [size(DJFTemp, 3),1]);
        LptRH = reshape(DJFrh(Llon_pt_row, Llat_pt_col, :), [size(DJFTemp, 3),1]);
        LptSH = reshape(DJFsh(Llon_pt_row, Llat_pt_col, :), [size(DJFTemp, 3),1]);

        WptWind = reshape(DJFWSpeed(Wlon_pt_row, Wlat_pt_col, :), [size(DJFTemp, 3),1]);
        WptTemp = reshape(DJFTemp(Wlon_pt_row, Wlat_pt_col, :), [size(DJFTemp, 3),1]);
        WptRH = reshape(DJFrh(Wlon_pt_row, Wlat_pt_col, :), [size(DJFTemp, 3),1]);
        WptSH = reshape(DJFsh(Wlon_pt_row, Wlat_pt_col, :), [size(DJFTemp, 3),1]);

        % Assessment of Foehn Events
        % Variables to find Days with Foehn Activity in a given season
        foehnDays = false(size(LptRH));
        FDsummerTag = false(size(LptRH)); % Summer - True; Winter - False.
        FDTempDiff = nan(size(LptTemp));
        FDTempDiff_wrt_wside = nan(size(LptTemp));
        FDSHDiff_wrt_wside = nan(size(LptTemp));
        FDWSDiff = nan(size(LptTemp));
        FDRHDiff = nan(size(LptTemp));

        % Start with T1
        tstep = 1;

        % Variable to count the number of events in a given season
        foehnEventCount = 0;
        
        % Loop through every timestep
        for i = 1:length(LptRH)
            % Check TStep
%             disp(num2str(tstep));
            
            % Check that the tstep value doesn't exceed the length of array
            if tstep == length(LptRH)
%                 disp('Ending the Loop');
                break
            end
            
            % Check for Foehn conditions between tstep & tstep+1 
            if (LptTemp(tstep+1)-LptTemp(tstep) >= TEMPTH &&...
                    LptWind(tstep+1)-LptWind(tstep) >= WSTH &&...
                    LptRH(tstep+1)-LptRH(tstep) <= RHTH)
    
                % If the condition is true => Foehn Activity
                % current tstep => pre foehn condition.
                % tstep+1 => first timestep with foehn condition
                preFStart = tstep;
    
                % From the pre condition, use the next tsteps to assess how
                % long the Foehn Condition exists
%                 disp('<strong>Found a Foehn Event.</strong>'); 
%                 disp('Checking the duration of the Foehn event...')
                
                % Add one to event count
                foehnEventCount = foehnEventCount + 1;
                foehnDayscount = 0;
    
                searchstep = preFStart;
                while (LptTemp(searchstep+1)-LptTemp(preFStart) >= TEMPTH &&...
                        LptWind(searchstep+1)-LptWind(preFStart) >= WSTH &&...
                        LptRH(searchstep+1)-LptRH(preFStart) <= RHTH)
                    foehnDays(searchstep+1) = true;

                    % Summer Tag
                    if month(DJFTime(searchstep+1)) == 1 || month(DJFTime(searchstep+1)) == 12 % || month(DJFTime(searchstep+1)) == 12
                        FDsummerTag(searchstep+1) = true;
%                         disp(DJFTime(searchstep+1));
                    end

                    % Calculate differences 
                    FDTempDiff(searchstep+1) = LptTemp(searchstep+1)-LptTemp(preFStart);
                    FDWSDiff(searchstep+1) = LptWind(searchstep+1)-LptWind(preFStart);
                    FDRHDiff(searchstep+1) = LptRH(searchstep+1)-LptRH(preFStart);

                    % Calculate temperature difference for the Foehn Day
                    % with respect to the temperature on windward side. 
                    FDTempDiff_wrt_wside(searchstep+1) = LptTemp(searchstep+1)-WptTemp(searchstep+1);
                    FDSHDiff_wrt_wside(searchstep+1) = LptSH(searchstep+1)-WptSH(searchstep+1);
    
                    % Add one to the day count each time the condition
                    % is met.
                    foehnDayscount = foehnDayscount + 1;
%                     disp(['Foehn Day at TStep :', num2str(searchstep+1)]);
    
                    % Move to next tstep
                    searchstep = searchstep+1;

                    % Break the loop if searchstep goes beyond the
                    % length of the selected timeperiod array
                    if searchstep == length(LptRH)
                        break
                    end
                end
    
                % Display information
%                 disp(['Information about Foehn Event # ', num2str(foehnEventCount)]);
%                 disp(['Pre-Foehn Day: ', datestr(DJFTime(tstep), 'dd-mmm-yyyy')]);
%                 disp(['Foehn Activity Started on: ', datestr(DJFTime(tstep+1), 'dd-mmm-yyyy')]);
%                 disp(['Foehn Activity Ended on: ', datestr(DJFTime(searchstep), 'dd-mmm-yyyy')]);
%                 disp(['Duration of Foehn Event: ', num2str(foehnDayscount), ' Day(s).']);
%                 disp(' ');
    
                % Push the loop to jump the identified days and find now
                % Foehn event
                tstep = searchstep;
            else
                tstep = tstep + 1;
            end
        end
        
        % Store foehnDays
        ptFoehnDays = [ptFoehnDays; foehnDays];

        % Store Max Temp Diff Values
        maxTemp(countn, caseNo) = max(FDTempDiff_wrt_wside);
        meanTemp(countn, caseNo) = mean(FDTempDiff_wrt_wside, 'omitnan');
        
        nonNanID = ~isnan(FDTempDiff_wrt_wside);
        nonNanTempDiff = FDTempDiff_wrt_wside(nonNanID);
        nonNanTempDiff(nonNanTempDiff < 0) = [];
        
        nonNanFDSummerTag = FDsummerTag(nonNanTempDiff>=0);
        ptSummerTag = [ptSummerTag; nonNanFDSummerTag];
        
        ptLeeTemp = [ptLeeTemp; LptTemp];
        ptTempDiff = [ptTempDiff; nonNanTempDiff];
        allTempDiff = [allTempDiff; nonNanTempDiff];
        grp = [grp, repelem(caseNo, length(nonNanTempDiff))];
        allSummerTag = [allSummerTag; nonNanFDSummerTag];


        countn = countn + 1;

        % Find all days when temperature on Leeward side was more than
        % temperature on Wind side. By plotting this, we will get to know
        % how well Foehn is captured.
        TempDiff = LptTemp-WptTemp;
        % False Positive: Its a Foehn Day and temp diff wrt wpt is negative - 1
        falsePositiveCond = foehnDays == 1 & TempDiff < 0;
        falsePos = TempDiff; 
        falsePos(falsePositiveCond == 0) = 0;

        % False Negative: Its a non-Foehn Day and temp diff wrt wpt is positive - 1
        falseNegativeCond = foehnDays == 0 & TempDiff > 0;
        falseNeg = TempDiff;
        falseNeg(falseNegativeCond == 0) = 0;
        
        % Plot Figure
        % Initiate a figure with required subplots
        nexttile(T); hold on;
        
%         t = tiledlayout(T,2,1,"TileSpacing","none","Padding","tight");
%         t.Layout.Tile = countn;
%         t.Layout.TileSpan = [1,1];
%         nexttile(t); hold on;

        % Foehn Days
        bar(DJFTime, foehnDays.*-20, 1,'FaceColor',[.75 .75 .75], 'FaceAlpha', 0.5);
        bar(DJFTime, foehnDays.*20, 1,'FaceColor',[.75 .75 .75], 'FaceAlpha', 0.5);
        
        % True Positive - Temp Diff % Blue
        bar(DJFTime, FDTempDiff_wrt_wside*1, 1,'FaceColor', [0 0 0.5], 'FaceAlpha', 0.8);
        % False Positive % Red
        bar(DJFTime, falsePos, 1,'FaceColor',[0.5 0 0], 'FaceAlpha', 0.8);
%         % False Negative % Red
%         bar(DJFTime, falseNeg, 1.5,'FaceColor',[0.5 0 0], 'FaceAlpha', 0.8);
        ylim([-5 +7])
        grid on; box off;

        % Text to print on Figure
        txt = [datestr(startdate, "mmm-YYYY"), ' to ', datestr(enddate-1, "mmm-YYYY")];
        text(0.05, 0.9, txt, "Units", "normalized", 'FontSize', 10)

%         nexttile(t);
%         yyaxis right;
%         plot(DJFTime, LptTemp, '-r');
%         ylim([-40 0]); grid on; box on;
    end

    allFoehnDays = [allFoehnDays, ptFoehnDays];
    allptLeeTemp = [allptLeeTemp, ptLeeTemp];
end

% legend('','', 'True Positive', 'False Positive', 'Temp')
% lgd = legend('Orientation',"horizontal");
% lgd.Layout.Tile = 'south';

%% WARMING POTENTIAL BOX PLOT USING WINDSIDE - LEESIDE TEMP

figure;
% s = scatter(grp, allTempDiff, 50, allTempDiff, 'filled');
s1 = scatter(grp, allTempDiff, (allSummerTag+1).*50, allSummerTag, 'filled', MarkerFaceAlpha='0.2', MarkerEdgeAlpha='0.2');
hold on;
cond = allSummerTag == 1;
s2 = scatter(grp(cond), allTempDiff(cond), (allSummerTag(cond)+1).*50, allSummerTag(cond), 'filled', MarkerFaceAlpha='0.7', MarkerEdgeAlpha='0.7');
xlim([0 6]); ylim([0 17]); clim([-0.5 1.5])
colormap("jet");

hold on;
boxplot(allTempDiff, grp, 'Whisker', inf, 'Colors', 'k');
grid on; box on;
xticks([1,2,3,4,5]);
xlabel([]);

% Calculate Number of days in Summer 
areaCode = 5;
disp(['Total Foehn Days: ', num2str(sum(grp==areaCode))]);
disp(['Summer Days: ', num2str(sum(grp==areaCode & allSummerTag'==1))]);
disp(['Summer Precentage: ', num2str(sum(grp==areaCode & allSummerTag'==1)/sum(grp==areaCode)*100, 3)]);

for areaCode=1:5
    % areaCode = 2;
    % Check Median and Stddev values, max and min values
    disp([num2str(median(allTempDiff(grp == areaCode)), 3), ' ± ', ...
        num2str(std(allTempDiff(grp == areaCode)), 3), ' (', ...
        num2str(min(allTempDiff(grp == areaCode)), 3), '-', ...
        num2str(max(allTempDiff(grp == areaCode)), 3), ')']);
end

%% BOXPLOT FOR TEMPERATURE DURING FOEHN DAYS

% Variables needed for this section
% ddatetime: Datetime array (2922 x 1)
% allFoehnDays: Binary Foehn Days Array (2922 X 6; one column for each point of interest)
% allptLeeTemp: Lee Side Temperature (2922 X 6)

% Identify Summer - December and January days
ddatetime = ddatetime(1:2922);
DJInd = month(ddatetime)==1 | month(ddatetime) == 12 | month(ddatetime) == 11;
DJLeeTemp = allptLeeTemp(DJInd, :);
DJFoehn = allFoehnDays(DJInd, :);

% Calculate mean summer temperature
mean_DJLeeTemp = mean(DJLeeTemp);

% Non-Summer
NDJInd = month(ddatetime)~=1 & month(ddatetime) ~= 12 & month(ddatetime) ~= 11;
NDJLeeTemp = allptLeeTemp(NDJInd,:);
NDJFoehn = allFoehnDays(NDJInd, :);

% Calculate mean NDJ temperature
mean_NDJLeeTemp = mean(NDJLeeTemp);

ttemp = [DJLeeTemp; NDJLeeTemp];
tfoehn = [DJFoehn; NDJFoehn];
tgrp = [repelem(0, length(DJLeeTemp)), repelem(1, length(NDJLeeTemp))]';
tcolor1 = [repelem("#ff7d04", length(DJLeeTemp)), repelem("#0285ff", length(NDJLeeTemp))]';
tcolor2 = [repelem(1, length(DJLeeTemp)), repelem(0, length(NDJLeeTemp))]';
tareaNames = ["Riiser Larsen", "Ekstrom", "Nivlisen", "Muninisen", "Roi Baudouin", "Reference"];

figure; colormap("jet");
tiledlayout(1,5,"TileSpacing","compact","Padding","compact");
for jj = 1:5
    nexttile; hold on;
    title(tareaNames(jj))

    % Plot all temps
    scatter(tgrp+1, ttemp(:,jj), 100, 'kx', MarkerFaceAlpha='0.05', MarkerEdgeAlpha='0.05');
    yline(mean_DJLeeTemp(jj), '--', 'Color', "k", 'LineWidth', 1);
    yline(mean_NDJLeeTemp(jj), '.', 'Color', "k", 'LineWidth', 1);
%     scatter(1, mean_DJLeeTemp(jj), 50, 'k.');
%     boxplot(ax2, ttemp(:,jj), tgrp, 'Whisker', inf, 'Colors', 'k'); hold on;

    % Plot temp on foehn days
    % 1 - Calculate
    foehnTemp = ttemp(:,jj);
    foehnID = tfoehn(:,jj);
    foehngrp = tgrp;
    
    foehnTemp = foehnTemp(foehnID == 1);
    foehngrp = foehngrp(foehnID == 1);
    foehnColor = tcolor2(foehnID == 1);

    % 2 - Plot
    scatter(foehngrp+1.30, foehnTemp, 100, foehnColor, 'filled', MarkerFaceAlpha='0.2', MarkerEdgeAlpha='0.2');
%     boxplot([ttemp(:,jj); foehnTemp], [(foehngrp+0.30); tgrp], 'Whisker', inf, 'Colors', 'k')

    clim([-0.5 1.5]);
    xlim([0.25 3.25]); ylim([-50 0])
    xticklabels([]);
    box on; grid on;

    if jj ~=1
        yticklabels([]);
    end
end

% Save plot 
exportgraphics(gcf,['D:\NPI_Work_Research\MARv2\figures\', 'Boxplot_DJNDJTemp_regoinWise.png'],'Resolution',300)


%% TEMPERATURE DURING MONTHLY FOEHN DAYS - FIGURE IN SUPPLEMENT. 

% Variables needed for this section
% ddatetime: Datetime array (2922 x 1)
% allFoehnDays: Binary Foehn Days Array (2922 X 6; one column for each point of interest)
% allptLeeTemp: Lee Side Temperature (2922 X 6)
% mean_DJLeeTemp: Mean Summer Temperature 
% tareaNames: Names of the areas of interest. 

% ddatetime till 31-12-2021 to make it consistent
ddatetime = ddatetime(1:2922);

% Monthly processing
figure; hold on;
tiledlayout(3,2,"TileSpacing","compact","Padding","compact");
for jj = 1:5
    jjTemp = allptLeeTemp(:, jj);
    jjFoehn = allFoehnDays(:, jj);
    
    % Initiate a figure
    nexttile; hold on;
    title(tareaNames(jj));

    % Mean Summer Temperature - Yline
    yline(mean_DJLeeTemp(jj), 'k--', LineWidth=1);

    % Plot Monthly data
    for months = 1:12
        monthInd = month(ddatetime)==months;
        monthLeeTemp = jjTemp(monthInd);
        monthFoehn = jjFoehn(monthInd);
        monthFoehnTemp = monthLeeTemp(monthFoehn == 1);
        
        scatter(repelem(months, length(monthLeeTemp)), monthLeeTemp, 50, 'kx', MarkerFaceAlpha='0.05', MarkerEdgeAlpha='0.05');
        scatter(repelem(months, length(monthFoehnTemp)), monthFoehnTemp, 150, 'MarkerFaceColor', [1,0.5,0.01], MarkerFaceAlpha='0.1', MarkerEdgeAlpha='0.1');
    end
    grid on; box on;
    ylim([-40 0]);xlim([0 13])
    
    % X Axes Labels
    xticks(1:12);
    xticklabels(["J", "F", "M", "A", "M", "J", ...
        "J", "A", "S", "O", "N", "D"]);
end

legend(["Mean Summer Temperature", "Daily Temperature", "Temperature During Foehn Days"])

% Save plot 
exportgraphics(gcf,['D:\NPI_Work_Research\MARv2\figures\', 'Monthly_FoehnDayTemperatres.png'],'Resolution',300)
