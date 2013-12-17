function FSI = fsi(spd,varargin)
% Full Spectrum Index
% Calculates the Full Spectrum Index according to the procedure documented
% in NLPIP Lighitng Answers: Specifying Color, 2004
% fsi = FSI(spd,startw,endw,incrementw
% spd is a single column or row vector and startw, endw and incrementw specify
% the starting, ending and increment of wavelength in nm.

% fsi = FSI(spd)
% spd is a two colunm matrix with wavelength values in column 1 and spd values in
% column 2. Weavelength values are assumed to be in units of nm.

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

% Interpolate to wavelength interval of 1 from 380nm to 730nm
numWave = 351;
t=(380:1:730)';
spd=interp1(wavelength_spd,spd,t,'spline');
spd(isnan(spd)) = 0.0;  
spd = spd/sum(spd); % Normalize the relative spd so that the total power equal 1

%Equal engery accumulative spd
EEspd=(1/numWave:1/numWave:1)';

%Calculate FSI
for j=1:numWave
    cum = cumsum(spd); % The built-in MatLab function cumsum replaces the commented-out for-loop below
    leastwo = (cum-EEspd).^2;
    %for i=1:numWave
    %    cum=sum(spd(1:i));% ./sum(spd); removed this term because spd is already normalized to a total power of 1.0
    %    leastwo(i)=(cum-EEspd(i)).^2;
    %end
    sumleastwo(j)=sum(leastwo);
    spd=circshift(spd,1);
end
FSI=mean(sumleastwo); 

%Normalize the FSI in such a way that equal energy source has a FSI of 100
%and the Low CRI HPS lamps has a FSI~0
%FSCI=100-5.1*FSI;
