function CAS120calibration
%CAS120CALIBRATION Generate calibration factors for CAS120 spectrometer
%   Uses SPD from reference lamp 42

%% Set file paths
refFile = '\\ROOT\projects\NASA-Review\Lamp42_HalfMeter_17Nov2011.txt';
cas120File = '\\ROOT\projects\NASA-Review\2014-01-22_spectrometerReadings\Lamp42_50cm.ISD';
calFile = '\\ROOT\projects\NASA-Review\2014-01-22_spectrometerReadings\cal.mat';

%% Import SPD from files
[refWav,refIrr] = importTXT(refFile);
[casWav,casIrr] = importISD(cas120File);

%% Use spline interpolation to resample reference to measured wavelengths
newRefIrr = spline(refWav,refIrr,casWav);

%% Generate calibration factor
calIrr = newRefIrr./casIrr;
calWav = casWav;

%% Save calibration factor file
save(calFile,'calWav','calIrr');


end

