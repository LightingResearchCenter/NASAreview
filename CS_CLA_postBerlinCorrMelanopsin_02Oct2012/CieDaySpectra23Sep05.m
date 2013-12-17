function spdDay = CieDaySpectra23Sep05(Tc,wave)
% Daylight Spectra
% Calculate the relative CIE Standard Daylight spectrum of given color temperature, Tc.
% Returns two vectors: wavelength in nanometers and the corresponding spectral values
% Function arguement:
% Tc - correlated color temperature in Kelvin

%wave = 350:5:1000;
if (Tc <= 25000)
	%load('Ciedaysn.mat','wavelength','S0','S1','S2');
    Table = load('CIEDaySN.txt');
    wavelength = Table(:,1);
    S0 = Table(:,2);
    S1 = Table(:,3);
    S2 = Table(:,4);
	if (Tc <= 7000)
		xd = -4.6070e9 / Tc.^3 + 2.9678e6 / Tc.^2 + 0.09911e3 / Tc + 0.244063;
	else
		xd = -2.0064e9 / Tc.^3 + 1.9018e6 / Tc.^2 + 0.24748e3 / Tc + 0.237040;
	end
	yd = -3.000*xd*xd + 2.870*xd - 0.275;
	M1 = (-1.3515 - 1.7703*xd + 5.9114*yd) / (0.0241 + 0.2562*xd - 0.7341*yd);
	M2 = (0.0300 - 31.4424*xd + 30.0717*yd) / (0.0241 + 0.2562*xd - 0.7341*yd);
	spdDay = S0 + M1*S1 + M2*S2;
    
    % Interpolate to desired wavelength intervals
    spdDay = interp1(wavelength,spdDay,wave,'pchip',0.0);
    spdDay(isnan(spdDay)) = 0.0;
else
	R = -1;
	return
end

