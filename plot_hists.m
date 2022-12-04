%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%16 September 2022
%Modified 02 October 2022

clear;
close all;

%load('variables_global_recalc.mat')
%run_data = run_data_recalc;

load('variables_global.mat')

clean_cost = cost(~isinf(cost) & ~isnan(cost));
% h1 = histogram(clean_cost);
% h1.EdgeColor = [0 0 0];
% h1.FaceColor = [0 0 1];
% h1.BinWidth = (max(clean_cost)-min(clean_cost))/10;
cost_cutoff = median(clean_cost);

i = 5;
itr = 1;
figure(itr)
h1 = histogram(run_data(:,i));
h1.EdgeColor = [0 0 0];
h1.FaceColor = [0.7 0.7 0.7];
h1.BinWidth = (max(run_data(~isinf(run_data(:,i)),i)) - min(run_data(:,i)))/10;
hold on;
h2 = histogram(run_data(cost<cost_cutoff,i));
h2.EdgeColor = [0 0 0];
h2.FaceColor = [1 0 0];
h2.BinWidth = (max(run_data(~isinf(run_data(:,i)),i)) - min(run_data(:,i)))/10;
itr = itr+1;