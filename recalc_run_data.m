%Haiyang Kehoe8
%University of Arizona
%Department of Geosciences
%12 September 2022
%Modified 20 October 2022

%recalc_run_data.m recalculates run_data.txt generated from a run through
%fidbp_auto.m. It does not update variables.mat nor variables_global.mat,
%nor the individual run_xxx_xxx.txt files.

clear;
close all;

tic

%Flag
write_txt = 0;

load('variables_global.mat')

%Remove GCMT constraints -HK 09/27/2022
c_lat = 0;
c_lon = 0;
c_dep = 0;
c_dur = 0;

%run_data_recalc = zeros(size(run_data));
for l = 1:size(run_data,1)
    tdirs_best = run_data(l,2);
    conts_best = run_data(l,3);
    rupv_acuts_best = run_data(l,4);
    %fprintf('Reading files with %3.2fs timestep, %3.2f IDBP contour,and %3.2f IDBP amplitude cutoff at %f seconds.\n',tdirs_best,conts_best,rupv_acuts_best,toc)
    ms_path=sprintf('%2.2fs/%s_xcor/',tdirs_best,dirname);
    as_path=sprintf('%2.2fs/%s_syn_xcor/',tdirs_best,dirname);
    as_pad_path=sprintf('%2.2fs/Pad_%s_%3.2f/',tdirs_best,as_path(7:end-1),conts_best);
    idbp_path=sprintf('%2.2fs/IDBP_%s_by_%s_%3.2f/',tdirs_best,ms_path(7:end-1),as_path(7:end-1),conts_best);
    load(strcat(idbp_path,'variables.mat'),'ms_hyp','idbp_data','cont_pad_in')
    [misfit,misfit2,vmisfit,vmisfit2,v_loc_misfit,v_loctime_misfit,grdvol_misfit,grdvoltime_misfit,cont_misfit,conttime_misfit,~,mean_vel,vs,az,rup_dur,rup_len,num_idbp,rupv_ratio,~,~,~,~,~,~,~,~,~] = ...
    calc_misfit_rupture_prop(ms_hyp,idbp_data,rupv_acuts_best,cont_pad_in,strike1,strike2,dip1,dip2,rake1,rake2,c_lat,c_lon,c_dep,c_dur);
    run_data_recalc(l,:) = [l,run_data(l,2),run_data(l,3),run_data(l,4),run_data(l,5),mean_vel/vs,az,rup_dur,rup_len,num_idbp,rupv_ratio,misfit,misfit2,vmisfit,vmisfit2,v_loc_misfit,v_loctime_misfit,grdvol_misfit,grdvoltime_misfit,cont_misfit,conttime_misfit,misfit/misfit2,vmisfit/vmisfit2];
end

%Calculate cost and save run_data textfile
[best_cost,best_idx,cost] = cost_func(idbp_data,run_data_recalc,'run_data_recalc.txt');
fprintf('Done recalculating run_data at %f seconds.\n',toc)
save('variables_global_recalc')

%Write associated text files
if write_txt==1
    write_txt_func('run_data.txt',dirname,best_idx,0,1)
end