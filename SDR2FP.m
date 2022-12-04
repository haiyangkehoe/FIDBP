%Modified by Haiyang Kehoe from https://github.com/djpugh/MTplot/blob/master/MTconvert/SDR2FP.m
%University of Arizona
%Department of Geosciences
%4 July 2022
%Modified  4 July 2022

% Convert strike (°), dip (°) and rake  (°) to the slip and normal for a plane (i.e. the normals to the fault and auxiliary planes)
% Right hand rule following x=east, y=north, z=up

function [slip, normal] = SDR2FP(strike,dip,rake)

slip   = [sind(strike).*cosd(rake)-cosd(strike).*cosd(dip).*sind(rake);
          cosd(strike).*cosd(rake)+sind(strike).*cosd(dip).*sind(rake);
          sind(dip).*sind(rake)];
normal = [cosd(strike).*sind(dip);
         -sind(strike).*sind(dip);
          cosd(dip)];