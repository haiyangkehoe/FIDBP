%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%13 June 2022
%Modified 10 November 2022

clear;
close all;

%Flags
plot_static = 1; %=0 to plot IDBP movie, =1 to plot static IDBP result
save_param  = 1; %=1 to save static IDBP parameters

%User inputs
run_data = dlmread('run_data.txt');
dirname = 'Aftershock_0';
idx = dlmread('best_idxs.txt','\t',[0 1 0 1]);

%User input
tdirs_best = run_data(idx,2);
conts_best = run_data(idx,3);
rupv_acuts_best = run_data(idx,4);
fprintf('Reading files with %3.2fs timestep, %3.2f IDBP contour,and %3.2f IDBP amplitude cutoff.\n',tdirs_best,conts_best,rupv_acuts_best)
path = sprintf('%2.2fs/IDBP_%s_xcor_by_%s_syn_xcor_%3.2f/',tdirs_best,dirname,dirname,conts_best);
ms_path = sprintf('%2.2fs/%s_xcor/',tdirs_best,dirname);

%GCMT Solution
gcmt = dlmread('gcmt.txt');
strike1 = gcmt(1);
strike2 = gcmt(2);
dip1    = gcmt(3);
dip2    = gcmt(4);
rake1   = gcmt(5);
rake2   = gcmt(6);
c_lat   = gcmt(7);
c_lon   = gcmt(8);
c_dep   = gcmt(9);
c_dur   = gcmt(10);

%Read in data
% hyp         = dlmread(strcat(path,'hyp.txt'));
% cont_pad_in = dlmread(strcat(path,'cont.txt'));
% cont_t      = cont_pad_in(:,1);
% cont_locs   = cont_pad_in(:,2:4);

%Find back-projection dimensions
% tot_pts = length(dlmread(strcat(path,'3D_0001.txt')));
% t_max   = length(dir([path '3D_*.txt']));

%Read back-projection data
% data = zeros(tot_pts*t_max,5);
% for t = 1:t_max
%     data(tot_pts*(t-1)+1:tot_pts*t,:) = [t*ones(tot_pts,1) dlmread(strcat(path,sprintf('3D_%04d.txt',t)),'')];
% end

%Load data
load(strcat(path,'variables.mat'),'ms_hyp')
hyp = ms_hyp;
load(strcat(path,'variables.mat'),'cont_pad_in')
cont_t      = cont_pad_in(:,1);
cont_locs   = cont_pad_in(:,2:4);

%Load back-projection dimensions
load(strcat(path,'variables.mat'),'ms_tot_pts')
tot_pts = ms_tot_pts;
load(strcat(path,'variables.mat'),'ms_t_max')
t_max = ms_t_max;

%Load back-projection data
load(strcat(path,'variables.mat'),'idbp_data')
data = idbp_data;

%Fault plane determination values
lond = round(mean(diff(unique(data(:,2)))),4);
latd = round(mean(diff(unique(data(:,3)))),4);
depd = round(mean(diff(unique(data(:,4)))),4);
min_rup_len_val = max([deg2km(lond) deg2km(latd) depd]) + 0.0001;
misfit_max = min_rup_len_val; %km
misfit2_max = min_rup_len_val; %km
misfit_ratio_max = 0.5;
vmisfit_max = min_rup_len_val; %km
vmisfit2_max = min_rup_len_val; %km
vmisfit_ratio_max = 0.5;

%Plot
if plot_static == 1
    parameters = fidbp_plot_static(hyp,data,rupv_acuts_best,cont_pad_in,strike1,strike2,dip1,dip2,rake1,rake2,c_lat,c_lon,c_dep,c_dur,misfit_max,misfit2_max,misfit_ratio_max,vmisfit_max,vmisfit2_max,vmisfit_ratio_max);
    R = 6371-parameters.mean_loc(3);
    rupv = [deg2km(parameters.mean_loc(4)-parameters.mean_loc(1),R);
            deg2km(parameters.mean_loc(5)-parameters.mean_loc(2),R);
                  -parameters.mean_loc(6)+parameters.mean_loc(3)];
    horiz = [0;0;-1];
    rupv_ang_from_horiz = acosd(dot(rupv,horiz)/(norm(rupv)*norm(horiz))) - 90;
    R = 6371-((parameters.mean_loc(3)+parameters.s_loc(3))/2);
    s_loc_dist = sqrt(deg2km(parameters.mean_loc(1)-parameters.s_loc(1),R).^2 + ...
                      deg2km(parameters.mean_loc(2)-parameters.s_loc(2),R).^2 + ...
                            (parameters.mean_loc(3)-parameters.s_loc(3)).^2);
    fprintf(sprintf('\nRupture properties:\nvr/vs = %3.2f\nvr = %3.2f\naz = %3.2f\nrup_dur = %3.2f\nrup_len = %3.2f\npref_strike = %i\npref_dip = %i\npref_rake = %i\nother_strike = %i\nother_dip = %i\nother_rake = %i\nRupture angle from horizontal = %3.2f\n',...
                    parameters.mean_vel/parameters.vs,...
                    parameters.mean_vel,...
                    parameters.az,...
                    parameters.rup_dur,...
                    parameters.rup_len,...
                    parameters.pref_strike,...
                    parameters.pref_dip,...
                    parameters.pref_rake,...
                    parameters.other_strike,...
                    parameters.other_dip,...
                    parameters.other_rake,...
                    rupv_ang_from_horiz))
    if save_param == 1
        file = fopen('parameters.txt', 'w');
        fprintf(file, '%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %3.2f %3.2f %3.2f %3.2f %3.2f %i %i %i %i %i %i %3.2f %3.2f %3.2f %3.2f\n',...
                      [parameters.mean_loc...
                       parameters.mean_vel...
                       parameters.mean_vel/parameters.vs...
                       parameters.az...
                       parameters.rup_dur...
                       parameters.rup_len...
                       parameters.pref_strike...
                       parameters.pref_dip...
                       parameters.pref_rake...
                       parameters.other_strike,...
                       parameters.other_dip,...
                       parameters.other_rake,...
                       rupv_ang_from_horiz,...
                       parameters.vangmisfit,...
                       parameters.vangmisfit2,...
                       parameters.vangmisfit/parameters.vangmisfit2]);
        fclose(file);
        file = fopen('parameters_slab2.txt', 'w');
        fprintf(file, '%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n',...
                       [parameters.s_loc,...
                       s_loc_dist,...
                       parameters.s_dip,...
                       parameters.s_str,...
                       parameters.s_thk]);
        fclose(file);
        idbp_radiators = [(parameters.data_idx_plane(:,1)-hyp(1)).*hyp(5) parameters.data_idx_plane(:,2:5)];
        dlmwrite(strcat(ms_path,'idbp_radiators.txt'),idbp_radiators,'delimiter','\t','precision',7);
        ms_tot_pts = size(data,1)/t_max;
        [~,closest_idx] = min(sum(abs(data(1:ms_tot_pts,2:4)-mean(idbp_radiators(:,2:4))),2));
        idbp_radiators_cent = data(closest_idx,2:4);
        dlmwrite(strcat(ms_path,'idbp_radiators_cent.txt'),idbp_radiators_cent,'delimiter','\t','precision',7);
    end
elseif plot_static == 0
    fidbp_plot(hyp,cont_t,cont_locs,data,rupv_acuts_best)
end