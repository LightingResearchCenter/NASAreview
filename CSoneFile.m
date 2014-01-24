function [L,x,y,CLA,CS] = CSoneFile
%CSONEFILE find the CS from one SPD file

addpath('CS_CLA_postBerlinCorrMelanopsin_02Oct2012');

%% Have user select file to load
defaultDir = fullfile([filesep,filesep,'ROOT'],'projects',...
    'NASA-Review','2014-01-22_spectrometerReadings');
[fileName,fileDir,filterSpec] = uigetfile({'*.isd';'*.txt'},'Select SPD to analyze',defaultDir);
filePath = fullfile(fileDir,fileName);

%% Import the selected file
if filterSpec == 1
    [wavelength,irradiance] = importISD(filePath);
    % Calibrate SPD
    calFile = fullfile(defaultDir,'cal.mat');
    irradiance = calibrateSPD(wavelength,irradiance,calFile);
elseif filterSpec == 2
    [wavelength,irradiance] = importTXT(filePath);
else
    error('Unknown file type.');
end

%% Merge SPD
spd = [wavelength(:),irradiance(:)];

%% Perform analysis
[L,x,y] = Lxy23Sep05(spd);
CLA = CLA_postBerlinCorrMelanopsin_02Oct2012(spd);
CS = CSCalc_postBerlin_12Aug2011(CLA);

end

