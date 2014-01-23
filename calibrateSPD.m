function irradiance = calibrateSPD(wavelength,irradianceIn,calFile)
%CALIBRATESPD Apply calibration from file to SPD

load(calFile,'calWav','calIrr')

if wavelength == calWav
    irradiance = irradianceIn.*calIrr;
else
    newCalIrr = spline(calWav,calIrr,wavelength);
    irradiance = irradianceIn.*newCalIrr;
end


end

