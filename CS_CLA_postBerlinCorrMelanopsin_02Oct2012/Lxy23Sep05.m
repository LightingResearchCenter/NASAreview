function [L,x,y] = Lxy23Sep05(spd,varargin)
% Photometry calculations
% Calculates luminous flux and (x,y) chromaticity coordinates
% Function arguements:
%	spd - spectral power distribution vector (when in units of Watts/nm th output, L, is in lumens)
%			(spd vector must be a column vector, size(n by 1))
%	startw - starting wavelength of spd in nanometers,
%	endw - ending wavelength of spd in nanometers,
%	incrementw - increment of wavelength data in nanometers.
% OR FOR NON-REGULAR WAVELENGTH INCREMENTS
% use a spd file with 2 columns vectors [wavelength value] in any increment, regular or not.
% When used with a 2-column spd arguement, startw, endw and increment arguements ar not used, but
% place holder values, such as zeros, must be specified. e.g., Tc = CCT_1(spd,0,0,0)

if isempty(varargin)
    [~, columns] = size(spd);
    if columns > 2
        error('Not column oriented data. Try transposing spd');
    end
    wavelength_spd = spd(:,1);
	spd = spd(:,2);
else
    startw = varargin{1};
    endw = varargin{2};
    incrementw = varargin{3};
    wavelength_spd = (startw:incrementw:endw)';
    [~, columns] = size(spd);
    if columns > 1
        error('Detected multiple columns of data. Try transposing spd');
    end
end

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

X = trapz(wavelength_spd,spd .* xbar);
Y = trapz(wavelength_spd,spd .* ybar);
Z = trapz(wavelength_spd,spd .* zbar);

L = 683*trapz(wavelength_spd',spd.*ybar);
x = X/(X+Y+Z);
y = Y/(X+Y+Z);
