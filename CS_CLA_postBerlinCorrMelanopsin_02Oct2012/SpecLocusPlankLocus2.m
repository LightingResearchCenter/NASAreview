%function [spectrumLocus,PlankLocus,isoTempLines] = SpecLocusPlankLocusFun()

%**** Put Temperature of Isotemperature Lines Here (K) ***********
Tiso = [2000 3000 4000 6000 8000];
% How long do you want the isotemperature lines?
CCTlimit = 0.03; % Distance in u',v' (CIE 1976) space

%data rounded to 4 decimal places 380 to 750 nm
%load('CIE31', 'wavelength','xbar','ybar','zbar');
% data to full precision 360 to 800 nm by 5 nm
%load xyz31.txt
%wavelength = (xyz31(:,1))';
%xbar = (xyz31(:,2))';
%ybar = (xyz31(:,3))';
%zbar = (xyz31(:,4))';
%Data to full precision 360 to 830 nm by 1 nm
%load('CIE31_1', 'wavelength','xbar','ybar','zbar');
Table = load('CIE31by1.txt');
wavelength = Table(:,1)';
xbar = Table(:,2)';
ybar = Table(:,3)';
zbar = Table(:,4)';
n = length(xbar);
% Generate spectrum locus for u,v space   
for i = 1:n
   spd = zeros(1,n);
   spd(i) = 1;
   X = sum(spd .* xbar);
	Y = sum(spd .* ybar);
	Z = sum(spd .* zbar);
   
   x(i) = X/(X+Y+Z);
   y(i) = Y/(X+Y+Z);
	u(i) = 4*X/(X+15*Y+3*Z);
   v(i) = 6*Y/(X+15*Y+3*Z);
   uprime(i) = u(i);
   vprime(i) = 1.5*v(i);
end
x(length(xbar)+1) = x(1);
y(length(xbar)+1) = y(1);
u(length(xbar)+1) = u(1);
v(length(xbar)+1) = v(1);
uprime(length(xbar)+1) = uprime(1);
vprime(length(xbar)+1) = vprime(1);

figure(1)
plot(x,y)
hold on
figure(2)
plot(u,v)
hold on
figure(3)
plot(uprime,vprime)
hold on

%spectrumLocus = [uprime' vprime'];
spectrumLocus = [x' y'];

%Calculate Plankian (blackbody) locus
for i = 1:(20000-1000)/50 
   Tc(i) = i*50+1000;
	wave = wavelength;
   c1 = 3.7418e-16;
	c2 = 1.4388e-2;
	spdref = c1 * (1e-9*wave).^-5 ./ (exp(c2./(Tc(i).* 1e-9*wave)) - 1);

	X = sum(spdref .* xbar);
	Y = sum(spdref .* ybar);
   Z = sum(spdref .* zbar);
   xb(i) = X/(X+Y+Z);
   yb(i) = Y/(X+Y+Z);
	ub(i) = 4*X/(X+15*Y+3*Z);
   vb(i) = 6*Y/(X+15*Y+3*Z);
   uprimeb(i) = ub(i);
   vprimeb(i) = vb(i)*1.5;
end
PlankLocus = [xb' yb'];

figure(1)
plot(xb,yb,'k-')
axis equal
title('CIE 1931 Chromaticity Diagram with Planktain Locus')
xlabel('x')
ylabel('y')
figure(2)
plot(ub,vb,'k--')
axis equal
title('CIE 1960 Chromaticity Diagram with Planktain Locus')
xlabel('u')
ylabel('v')
figure(3)
plot(uprimeb,vprimeb,'k--')
axis equal
title('CIE 1976 Chromaticity Diagram with Planktain Locus')
xlabel('u`')
ylabel('v`')

%*********** Calculation of Isotemperature Lines ******************************

ubar = (2/3)*xbar;
vbar = ybar;
wbar = -0.5*xbar + (3/2)*ybar + 0.5*zbar;

dwave = 1;

for i = 1:length(Tiso)
   spdref = c1 * (1e-9*wavelength).^-5 ./ (exp(c2./(Tiso(i).* 1e-9*wavelength)) - 1);
   spdref = spdref/max(spdref);
   wave = wavelength*1e-9;
   
   U = sum(spdref.*ubar);
   V = sum(spdref.*vbar);
   W = sum(spdref.*wbar);
   R = U+V+W;
   u(i) = U/R;
   v(i) = V/R;
      % take derivative of Plankian locus
   Uprime = c1*c2*(Tiso(i))^-2*sum(wave.^-6.*ubar.*exp(c2./(wave.*Tiso(i))).*(exp(c2./(wave.*(Tiso(i))))-1).^-2)*dwave;
   Vprime = sum(c1*c2*Tiso(i)^-2*wave.^-6.*vbar.*exp(c2./(wave.*Tiso(i))).*(exp(c2./(wave.*(Tiso(i))))-1).^-2)*dwave;
   Wprime = sum(c1*c2*Tiso(i)^-2*wave.^-6.*wbar.*exp(c2./(wave.*Tiso(i))).*(exp(c2./(wave.*(Tiso(i))))-1).^-2)*dwave;
   Rprime = Uprime+Vprime+Wprime;
   
   sl = (Vprime*R-V*Rprime)/(Uprime*R-U*Rprime);
   m = -1/sl;
   
   uLine = u(i);
   vLine = v(i);
   index = 0;
   clear uLineV;
   clear vLineV;
   while (sqrt((uLine-u(i))^2+(9/4)*(vLine-v(i))^2) < CCTlimit)
	   index = index + 1;
	   uLine = uLine + 0.0001;
	   vLine = m*uLine + (v(i)-m*u(i));
	   uLineV(index) = uLine;
	   vLineV(index) = vLine;
   end
   uLine = u(i);
   vLine = v(i);
   uLineV = fliplr(uLineV);
   vLineV = fliplr(vLineV);
   while (sqrt((uLine-u(i))^2+(9/4)*(vLine-v(i))^2) < CCTlimit)
	   index = index + 1;
	   uLine = uLine - 0.0001;
	   vLine = m*uLine + (v(i)-m*u(i));
	   uLineV(index) = uLine;
	   vLineV(index) = vLine;
   end
   isoTempLines = [uLineV',vLineV'];
   
  
   figure(2)
   hold on
   plot(uLineV,vLineV,'r-')
   hold off
   figure(1)
   hold on
   xLine = 27/4*uLineV./(9/2*uLineV - 18*vLineV + 9);
   yLine =  9/2*vLineV./(9/2*uLineV - 18*vLineV + 9);
   plot(xLine,yLine,'r-')
   hold off
   figure(3)
   hold on
   uPrimeLine = uLineV;
   vPrimeLine = vLineV*1.5;
   plot(uPrimeLine, vPrimeLine,'r-')
   hold off
	   
end

px = [.445,.452,.326,.279,.146];
py = [.419,.415,.349,.273,.0583];
labels = {'warm white','dim warm white','cool white','phase shift','blue'};
puprime = 4*px./(-2*px+12*py+3);
pvprime = 9*py./(-2*px+12*py+3);
pu = puprime;
pv = 2/3*pvprime;

figure(1)
hold on
h1 = plot(px,py,'bd');
set(h1,'MarkerFaceColor',[0,0,1]);
hold off
for i1 = 1:length(labels)
    text(px(i1)+.02,py(i1),labels(i1))
end
figure(2)
hold on
plot(pu,pv,'bs')
hold off
figure(3)
hold on
plot(puprime,pvprime,'bs')
hold off

