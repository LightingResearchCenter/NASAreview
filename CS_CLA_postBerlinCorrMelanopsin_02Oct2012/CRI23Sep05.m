function [R,DC,Tc] = CRI23Sep05(spd,varargin)
% CRI Color Rendering Indices (1-8)
% Calculates the Color Rendering Indices (CRI) for the first 8 CIE
% test color samples according to CIE 13.3 (1995).
% Interpolates xbar, ybar, and zbar values to same increment as spd data. 
% The following files are needed:
%	CIE31_1.mat, CIEDaySn.mat, CCT_1.m, isoTempLines.mat, CIEtcsa.mat
% Function arguments are:
%	spd - spectral power distribution vector,
%	startw - starting wavelength of spd in nanometers,
%	endw - ending wavelength of spd in nanometers,
%	incrementw - increment of wavelength data in nanometers.
%  OR
%   spd - spectral power distribution 2-column matrix
%         column 1 = wavelength in nm (can be any arbitrary spacing)
%         column 2 = intensity values
%   put in place holders for startw, endw and incrementw, e.g.  [R,DC,Tc] = CRI_1(spd_matrix,0,0,0)
% Function returns:
%   R = vector of the color rendering indicies for each of the eight test color samples
%   DC = scalar tolorance value having to do with the appropriateness of applying CRI to colored light sources (see CIE documentation)
%   Tc = The correlated color temperature of the light source in Kelvin units

if length(varargin)==0
    [rows columns] = size(spd);
    if columns > 2
        error('Not column oriented data. Try transposing spd');
    end
    wavelength_spd = spd(:,1);
	spd = spd(:,2);
else
    startw = varargin{1}
    endw = varargin{2}
    incrementw = varargin{3}
    wavelength_spd = (startw:incrementw:endw)';
    [rows columns] = size(spd);
    if columns > 1
        error('Detected multiple columns of data. Try transposing spd');
    end
end

% Calculate Correlated Color Temperature, Tc.
Tc = CCT23Sep05([wavelength_spd,spd]);

%load('CIE31_1', 'wavelength','xbar','ybar','zbar');
Table = load('CIE31by1.txt');
wavelength = Table(:,1);
xbar = Table(:,2);
ybar = Table(:,3);
zbar = Table(:,4);

xbar = interp1(wavelength,xbar,wavelength_spd);
xbar(isnan(xbar)) = 0.0;
ybar = interp1(wavelength,ybar,wavelength_spd);
ybar(isnan(ybar)) = 0.0;
zbar = interp1(wavelength,zbar,wavelength_spd);
zbar(isnan(zbar)) = 0.0;

% Calculate Reference Source Spectrum, spdref.
if (Tc < 5000)
    spdref = blackBodySpectra23Sep05(Tc,wavelength_spd);
else
	if (Tc <= 25000)
		spdref = CieDaySpectra23Sep05(Tc,wavelength_spd);      
	else
		error('CCT above 25,000 K');
	end
end

TCS = load('Tcs14_23Sep09.txt'); % 'Tcs14.txt'
%Interpolate TCS values from 5 nm to spd nm increments
TCS_1 = zeros(length(wavelength_spd),14);
wavelength_5 = TCS(:,1);
TCS = TCS(:,2:end); % Remove wavelength column
TCS = TCS/1000;
for i = 1:14
	TCS_1(:,i) = interp1(wavelength_5,TCS(:,i),wavelength_spd,'linear',0);
	%TCS_1(isnan(TCS_1(:,i)),i) = 0.0; % remove NaN from vector, but it
	%should't have any since a zero value was specified for extrapolated
	%values in above interp1 function.
end

% Calculate u, v chromaticity coordinates of samples under test illuminant, uk, vk and
% reference illuminant, ur, vr.
uki = zeros(1,14);
vki = zeros(1,14);
uri = zeros(1,14);
vri = zeros(1,14);
X = trapz(wavelength_spd,spd .* xbar);
Y = trapz(wavelength_spd,spd .* ybar);
Z = trapz(wavelength_spd,spd .* zbar);
Yknormal = 100 / Y;
Yk = Y*Yknormal;
uk = 4*X/(X+15*Y+3*Z);
vk = 6*Y/(X+15*Y+3*Z);
X = trapz(wavelength_spd,spdref .* xbar);
Y = trapz(wavelength_spd,spdref .* ybar);
Z = trapz(wavelength_spd,spdref .* zbar);
Yrnormal = 100 / Y;
Yr = Y*Yrnormal;
ur = 4*X/(X+15*Y+3*Z);
vr = 6*Y/(X+15*Y+3*Z);
for i = 1:14
	X = trapz(wavelength_spd,spd .* TCS_1(:,i) .* xbar);
	Y = trapz(wavelength_spd,spd .* TCS_1(:,i) .* ybar);
	Z = trapz(wavelength_spd,spd .* TCS_1(:,i) .* zbar);
	Yki(i) = Y*Yknormal;
	uki(i) = 4*X/(X+15*Y+3*Z);
	vki(i) = 6*Y/(X+15*Y+3*Z);
	X = trapz(wavelength_spd,spdref .* TCS_1(:,i) .* xbar);
	Y = trapz(wavelength_spd,spdref .* TCS_1(:,i) .* ybar);
	Z = trapz(wavelength_spd,spdref .* TCS_1(:,i) .* zbar);
	Yri(i) = Y*Yrnormal;
	uri(i) = 4*X/(X+15*Y+3*Z);
	vri(i) = 6*Y/(X+15*Y+3*Z);
end
% Check tolorence for reference illuminant
DC = sqrt((uk-ur).^2 + (vk-vr).^2);

% Apply adaptive (perceived) color shift.
ck = (4 - uk - 10*vk) / vk;
dk = (1.708*vk + 0.404 - 1.481*uk) / vk;
cr = (4 - ur - 10*vr) / vr;
dr = (1.708*vr + 0.404 - 1.481*ur) / vr;

for i = 1:14
	cki = (4 - uki(i) - 10*vki(i)) / vki(i);
	dki = (1.708*vki(i) + 0.404 - 1.481*uki(i)) / vki(i);
	ukip(i) = (10.872 + 0.404*cr/ck*cki - 4*dr/dk*dki) / (16.518 + 1.481*cr/ck*cki - dr/dk*dki);
	vkip(i) = 5.520 / (16.518 + 1.481*cr/ck*cki - dr/dk*dki);
end

%  Transformation into 1964 Uniform space coordinates.
for i = 1:14
	Wstarr(i) = 25*Yri(i).^.333333 - 17;
	Ustarr(i) = 13*Wstarr(i)*(uri(i) - ur);
	Vstarr(i) = 13*Wstarr(i)*(vri(i) - vr);
	
	Wstark(i) = 25*Yki(i).^.333333 - 17;
	Ustark(i) = 13*Wstark(i)*(ukip(i) - ur); % after applying the adaptive color shift, u'k = ur
	Vstark(i) = 13*Wstark(i)*(vkip(i) - vr); % after applying the adaptive color shift, v'k = vr
end

% Determination of resultant color shift, delta E.
deltaE = zeros(1,14);
R = zeros(1,14);
for i = 1:14
	deltaE(i) = sqrt((Ustarr(i) - Ustark(i)).^2 + (Vstarr(i) - Vstark(i)).^2 + (Wstarr(i) - Wstark(i)).^2);
	R(i) = 100 - 4.6*deltaE(i);
end
%Ra = sum(R(1:8))/8;
