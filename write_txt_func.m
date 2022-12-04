%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%11 September 2022
%Modified 27 September 2022

%write_txt.m writes back-projection and IDBP textfiles from variables.mat files

function write_txt_func(run_data_file,dirname,idx,as_pad_out_best,save_idbp_best)

tic

%Flags
%as_pad_out_best = 0;
%save_idbp_best = 1;

%User inputs
run_data = dlmread(run_data_file);
% dirname = 'Aftershock_00';
% idx = 241;

tdirs_best = run_data(idx,2);
conts_best = run_data(idx,3);
fprintf('Reading files with %3.2fs timestep and %3.2f IDBP contour.\n',tdirs_best,conts_best)
ms_path=sprintf('%2.2fs/%s_xcor/',tdirs_best,dirname);
as_path=sprintf('%2.2fs/%s_syn_xcor/',tdirs_best,dirname);
as_pad_path=sprintf('%2.2fs/Pad_%s_%3.2f/',tdirs_best,as_path(7:end-1),conts_best);
idbp_path=sprintf('%2.2fs/IDBP_%s_by_%s_%3.2f/',tdirs_best,ms_path(7:end-1),as_path(7:end-1),conts_best);
load(strcat(idbp_path,'variables.mat'))

%Write textfiles for the lowest misfit runs
if as_pad_out_best==1
    if exist(as_pad_path, 'dir')~=7
        mkdir(as_pad_path)
    elseif isempty(dir(strcat(as_pad_path,'*.txt')))==0 || ...
           isempty(dir(strcat(as_pad_path,'*.pdf')))==0
        fprintf(sprintf('Directory ''%s'' contains protected files.\n',as_pad_path))
        fprintf('Exiting...\n')
        return
    end
    %Write optimized reference back-projection result to text files
    fprintf('Writing optimized reference back-projection to text files at %f seconds.\n',toc)
    dlmwrite(strcat(as_pad_path,'hyp.txt'),as_pad_hyp,'delimiter','\t')
    as_pad_amps_sum = zeros(as_pad_tot_pts,as_pad_t_max);
    for t = 1:as_pad_t_max
        dlmwrite(strcat(as_pad_path,sprintf('3D_%04d.txt',t)), as_data_pad(as_pad_tot_pts*(t-1)+1:as_pad_tot_pts*t,2:5),'delimiter',' ','precision','%.4f')
        as_pad_amps_sum(:,t) = as_data_pad(as_pad_tot_pts*(t-1)+1:as_pad_tot_pts*t,5);
    end
    as_pad_data_sum = cat(2,as_data_pad(1:as_pad_tot_pts,2:4),sum(as_pad_amps_sum,2)/max(sum(as_pad_amps_sum,2)));
    dlmwrite(strcat(as_pad_path,'sum.txt'), as_pad_data_sum,'delimiter',' ','precision','%.4f')
    fprintf('Optimized reference back-projection completed at %f seconds.\n',toc)
end
if save_idbp_best==1
    if exist(idbp_path, 'dir')~=7
        mkdir(idbp_path)
    elseif isempty(dir(strcat(idbp_path,'*.txt')))==0 || ...
           isempty(dir(strcat(idbp_path,'*.pdf')))==0
        fprintf(sprintf('Directory ''%s'' contains protected files.\n',idbp_path))
        fprintf('Exiting...\n')
        return
    end
    fprintf('Writing inversion result to text files at %f seconds.\n',toc)
    dlmwrite(strcat(ms_path,'hyp.txt'),ms_hyp,'delimiter','\t')
    dlmwrite(strcat(as_path,'hyp.txt'),as_hyp,'delimiter','\t')
    dlmwrite(strcat(idbp_path,'hyp.txt'),ms_hyp,'delimiter','\t')

    %Write static plot files
    ms_amps_sum = zeros(ms_tot_pts,ms_t_max);
    for t = 1:ms_t_max
        ms_amps_sum(:,t) = ms_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5);
    end
    ms_data_sum = cat(2,ms_data(1:ms_tot_pts,2:4),sum(ms_amps_sum,2)/max(sum(ms_amps_sum,2)));
    dlmwrite(strcat(ms_path,'sum.txt'), ms_data_sum,'delimiter',' ','precision','%.4f')

    as_amps_sum = zeros(as_tot_pts,as_t_max);
    for t = 1:as_t_max
        as_amps_sum(:,t) = as_data(as_tot_pts*(t-1)+1:as_tot_pts*t,5);
    end
    as_data_sum = cat(2,as_data(1:as_tot_pts,2:4),sum(as_amps_sum,2)/max(sum(as_amps_sum,2)));
    dlmwrite(strcat(as_path,'sum.txt'), as_data_sum,'delimiter',' ','precision','%.4f')

    idbp_amps_sum = zeros(ms_tot_pts,ms_t_max);
    for t = 1:ms_t_max
        dlmwrite(strcat(idbp_path,sprintf('3D_%04d.txt',t)), idbp_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,2:5),'delimiter',' ','precision','%.4f')
        idbp_amps_sum(:,t) = idbp_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5);
    end
    idbp_data_sum = cat(2,idbp_data(1:ms_tot_pts,2:4),sum(idbp_amps_sum,2)/max(sum(idbp_amps_sum,2)));
    dlmwrite(strcat(idbp_path,'sum.txt'), idbp_data_sum,'delimiter',' ','precision','%.4f')
    dlmwrite(strcat(idbp_path,'cont.txt'), [cont_t cont_locs],'delimiter',' ','precision','%.4f')
end