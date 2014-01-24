function CSbatch
%CSBATCH find the CS from SPD of all .isd files in a folder

addpath('CS_CLA_postBerlinCorrMelanopsin_02Oct2012');

%% Have user select folder to load
startDir = fullfile([filesep,filesep,'ROOT'],'projects',...
    'NASA-Review','2014-01-22_spectrometerReadings');
fileDir = uigetdir(startDir,'Select folder to analyze');
dirName = regexprep(fileDir,'^.*[^\w-](\w*)$','$1');
calFile = fullfile(startDir,'cal.mat');

%% Find .isd files in folder
dirListing = dir(fullfile(fileDir,'*.isd'));

%% Import and process all files
nFiles = numel(dirListing);
% Preallocate variables
filePath = cell(nFiles,1);
wavelength = cell(nFiles,1);
irradiance = cell(nFiles,1);
spd = cell(nFiles,1);
Lux = zeros(nFiles,1);
x = zeros(nFiles,1);
y = zeros(nFiles,1);
CLA = zeros(nFiles,1);
CS = zeros(nFiles,1);
sampleTime = zeros(nFiles,1);
% Begin loop
for i1 = 1:nFiles
    % Import the current file
    filePath{i1} = fullfile(fileDir,dirListing(i1).name);
    [wavelength{i1},irradiance{i1}] = importISD(filePath{i1});
    % Calibrate SPD
    irradiance{i1} = calibrateSPD(wavelength{i1},irradiance{i1},calFile);
    % Merge SPD
    spd{i1} = [wavelength{i1}(:),irradiance{i1}(:)];
    % Perform analysis
    [Lux(i1),x(i1),y(i1)] = Lxy23Sep05(spd{i1});
    CLA(i1) = CLA_postBerlinCorrMelanopsin_02Oct2012(spd{i1});
    CS(i1) = CSCalc_postBerlin_12Aug2011(CLA(i1));
    % Extract sample time
    sampleTime(i1) = datenum(dirListing(i1).date);
end

%% Display results
figure;
hold on
[AX,~,~] = plotyy(sampleTime,CS,sampleTime,[Lux,CLA],'plot','semilogy');
datetick2;
title(dirName);
xlabel('sample time');
set(get(AX(1),'Ylabel'),'String','Circadian Stimulus (CS)') 
set(get(AX(2),'Ylabel'),'String','Illuminance (Lux)')
legend('CS','Lux','CLA','Location','SouthOutside','Orientation','Horizontal');

end

