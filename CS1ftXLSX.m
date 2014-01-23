function [ L,x,y,CLA,CS ] = CS1ftXLSX
%CS1FTXLSX Summary of this function goes here
%   Detailed explanation goes here

addpath('CS_CLA_postBerlinCorrMelanopsin_02Oct2012');

spd = importXLSX;

samples = numel(spd);

L = zeros(1,samples);
x = zeros(1,samples);
y = zeros(1,samples);
CLA = zeros(1,samples);
CS = zeros(1,samples);

for i1 = 1:samples
    [L(i1),x(i1),y(i1)] = Lxy23Sep05(spd{i1});
    CLA(i1) = CLA_postBerlinCorrMelanopsin_02Oct2012(spd{i1});
    CS(i1) = CSCalc_postBerlin_12Aug2011(CLA(i1));
end

end

