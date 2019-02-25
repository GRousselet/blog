function [ ColorArray ] = ColorCoder( n )
% Generates a nx3 numerical array where ColorArray(i,:) are the RGB values
% that the user choses
%   version: <1.00> from 30/11/2015
% 
%       Manuel Lera Ramírez: manulera14@gmail.com
%
%       Developed during Master Thesis Internship in Marcos González-Gaitán's Lab
%       Departments of Biochemistry and Molecular Biology, Sciences II, 
%       30 Quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
%--------------------------------------------------------------------------
ColorArray=nan(n,3);
for i=1:size(ColorArray,1)
    ColorArray(i,:)=uisetcolor;    
end
end

