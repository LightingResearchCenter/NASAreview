function results = CSbatch
%CSBATCH find the CS from SPD of all .isd files in a folder

addpath('CS_CLA_postBerlinCorrMelanopsin_02Oct2012');

%% Have user select folder to load
startDir = fullfile([filesep,filesep,'ROOT'],'projects',...
    'NASA-Review');
fileDir = uigetdir(startDir,'Select folder to analyze');
dirName = regexprep(fileDir,'^.*[^\w-](\w*)$','$1');
calFile = fullfile(startDir,'cal.mat');

%% Find .isd files in folder
dirListing = dir(fullfile(fileDir,'*.isd'));

%% Import and process all files
nFiles = numel(dirListing);
% Preallocate variables
results = dataset;
results.condition = cell(nFiles,1);
filePath = cell(nFiles,1);
wavelength = cell(nFiles,1);
irradiance = cell(nFiles,1);
spd = cell(nFiles,1);
results.Lux = zeros(nFiles,1);
results.x = zeros(nFiles,1);
results.y = zeros(nFiles,1);
results.CCT = zeros(nFiles,1);
results.CLA = zeros(nFiles,1);
results.CS = zeros(nFiles,1);
% Begin loop
for i1 = 1:nFiles
    % Set results.condition name
    results.condition{i1} = dirListing(i1).name(1:end-4);
    % Import the current file
    filePath{i1} = fullfile(fileDir,dirListing(i1).name);
    [wavelength{i1},irradiance{i1}] = importISD(filePath{i1});
    % Calibrate SPD
    irradiance{i1} = calibrateSPD(wavelength{i1},irradiance{i1},calFile);
    % Merge SPD
    spd{i1} = [wavelength{i1}(:),irradiance{i1}(:)];
    % Perform analysis
    [results.Lux(i1),results.x(i1),results.y(i1)] = Lxy23Sep05(spd{i1});
    results.CCT(i1) = CCT23Sep05(spd{i1});
    results.CLA(i1) = CLA_postBerlinCorrMelanopsin_02Oct2012(spd{i1});
    results.CS(i1) = CSCalc_postBerlin_12Aug2011(results.CLA(i1));
end

end

