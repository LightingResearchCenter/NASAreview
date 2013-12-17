function spdBlackBody = blackBodySpectra23Sep05(Tc,wave)
% Black Body Spectra
% Calculates the Planktian black body spectrum of given color temperature, Tc.
% Function arguements are:
%	Tc - color temperature in Kelvin
%	wave - column vector specifiying the wavelength values at which the spd is evaluated

% 2002 CODATA recommended values
h = 6.6260693e-34;
c = 299792458;
k = 1.3806505e-23;

c1 = 2*pi*h*c^2;
c2 = h*c/k;

spdBlackBody = c1 * (1e-9*wave).^-5 ./ (exp(c2./(Tc.* 1e-9*wave)) - 1);
