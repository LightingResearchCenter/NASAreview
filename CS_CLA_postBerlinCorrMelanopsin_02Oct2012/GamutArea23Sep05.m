function Areauvprime = GamutArea23Sep05(spd,varargin)

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

%********************* CRI Test Color Samples Chromaticity Calc ****************
TCS = load('Tcs14.txt');
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
xki = zeros(1,8);
yki = zeros(1,8);
uki = zeros(1,8);
vki = zeros(1,8);
ukiprime = zeros(1,8);
vkiprime = zeros(1,8);
for i = 1:8
    X = sum(spd .* TCS_1(:,i) .* xbar);
    Y = sum(spd .* TCS_1(:,i) .* ybar);
    Z = sum(spd .* TCS_1(:,i) .* zbar);
    xki(i) = X/(X+Y+Z);
    yki(i) = Y/(X+Y+Z);
    uki(i) = 4*X/(X+15*Y+3*Z);
    vki(i) = 6*Y/(X+15*Y+3*Z);
    ukiprime(i) = uki(i);
    vkiprime(i) = vki(i)*1.5;
end
%figure(1)
%hold on
%xki = [xki xki(1)];
%yki = [yki yki(1)];
%plot(xki,yki,'r*-')
%Areaxy = polyarea(xki,yki);
%hold off

%figure(2)
%hold on
%uki = [uki uki(1)];
%vki = [vki vki(1)];
%plot(uki,vki,'r*-')
%Areauv = polyarea(uki,vki);
%hold off

%figure(3)
%hold on
%ukiprime = [ukiprime ukiprime(1)];
%vkiprime = [vkiprime vkiprime(1)];
%plot(ukiprime,vkiprime,'ks-')
Areauvprime = polyarea(ukiprime,vkiprime);
%hold off
