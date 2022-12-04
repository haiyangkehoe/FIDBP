%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%26 August 2021
%Modified 15 November 2022
 
tic

%fidbp is FAST image deconvolution back-projection

clear;
close all;

%Flags
check_in = 0; %=1 to check inputs
as_pad_out = 0; %=1 to write padded aftershock files
normamp = 0; %=1 to select contour after normalizing each time step, =0 to select contour across entire normalized back-projection
plot_idbp = 0; %=1 to plot IDBP movie in MATLAB
save_idbp = 1; %=1 to save IDBP textfiles
save_var = 1; %=1 to save MATLAB variables
as_pad_out_best = 0; %=1 to write lowest misfit padded aftershock files
save_idbp_best = 0; %=1 to save lowest misfit IDBP textfiles

%User inputs
dirname = 'Aftershock_0';
tdirs = [0.2 0.25];
conts = sort(0.05:0.05:0.95,'descend');
rupv_acuts = 0.05:0.05:0.95;

pad_loc_lat_1 = 1; %Grid points to pad S
pad_loc_lat_2 = 1; %Grid points to pad N
pad_loc_lon_1 = 1; %Grid points to pad W
pad_loc_lon_2 = 1; %Grid points to pad E
pad_loc_dep_1 = 1; %Grid points to pad up
pad_loc_dep_2 = 1; %Grid points to pad down
pad_t_1 = 0; %Timesteps to pad in negative time
pad_t_2 = 0; %Timesteps to pad in positive time
amp_cutoff = 0.001; %Select minimum amplitude to plot

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

itr=1;
for l = 1:length(tdirs)
    breakinner=0;
    ms_path=sprintf('%2.2fs/%s_xcor/',tdirs(l),dirname);
    as_path=sprintf('%2.2fs/%s_syn_xcor/',tdirs(l),dirname);
    ms_ts=tdirs(l);
    as_ts=tdirs(l);
    ms_ot = round((15-2)/ms_ts)+1; %Mainshock origin timestep
    as_ot = round((20-2)/as_ts)+1; %Reference event origin timestep
            
    %Find back-projection dimensions
    ms_tot_pts = length(dlmread(strcat(ms_path,'3D_0001.txt')));
    ms_t_max = length(dir([ms_path '3D_*.txt']));

    as_tot_pts = length(dlmread(strcat(as_path,'3D_0001.txt')));
    as_t_max = length(dir([as_path '3D_*.txt']));
    
    for ll = 1:length(conts)
        if breakinner==1
            continue
        end
        as_pad_path=sprintf('%2.2fs/Pad_%s_%3.2f/',tdirs(l),as_path(7:end-1),conts(ll));
        idbp_path=sprintf('%2.2fs/IDBP_%s_by_%s_%3.2f/',tdirs(l),ms_path(7:end-1),as_path(7:end-1),conts(ll));
        cont = conts(ll);

        %Create write directories
        if save_idbp==1 || save_var==1
            if exist(idbp_path, 'dir')~=7
                mkdir(idbp_path)
            elseif isempty(dir(strcat(idbp_path,'*.txt')))==0 || ...
                   isempty(dir(strcat(idbp_path,'*.pdf')))==0
                fprintf(sprintf('Directory ''%s'' contains protected files.\n',idbp_path))
                fprintf('Exiting...\n')
                return
            else
                rmdir(idbp_path,'s')
                mkdir(idbp_path)
            end
        end
        if as_pad_out==1
            if exist(as_pad_path, 'dir')~=7
                mkdir(as_pad_path)
            elseif isempty(dir(strcat(as_pad_path,'*.txt')))==0 || ...
                   isempty(dir(strcat(as_pad_path,'*.pdf')))==0
                fprintf(sprintf('Directory ''%s'' contains protected files.\n',as_pad_path))
                fprintf('Exiting...\n')
                return
            else
                rmdir(as_pad_path,'s')
                mkdir(as_pad_path)
            end
        end

        %Read in mainshock back-projection
        ms_data = zeros(ms_tot_pts*ms_t_max,5);
        for t = 1:ms_t_max
            ms_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,:) = [t*ones(ms_tot_pts,1) dlmread(strcat(ms_path,sprintf('3D_%04d.txt',t)),'')];
            %Normalize each time step
            %ms_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5) = ms_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5)/max(ms_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5));
        end

        %Read in reference event back-projection
        as_data = zeros(as_tot_pts*as_t_max,5);
        for t = 1:as_t_max
            as_data(as_tot_pts*(t-1)+1:as_tot_pts*t,:) = [t*ones(as_tot_pts,1) dlmread(strcat(as_path,sprintf('3D_%04d.txt',t)),'')];
            %Normalize each time step
            %as_data(as_tot_pts*(t-1)+1:as_tot_pts*t,5) = as_data(as_tot_pts*(t-1)+1:as_tot_pts*t,5)/max(as_data(as_tot_pts*(t-1)+1:as_tot_pts*t,5));
        end

        %Check grid spacing and timesteps for mainshock and reference back-projection are the same
        ms_lond = round(mean(diff(unique(ms_data(:,2)))),4);
        ms_latd = round(mean(diff(unique(ms_data(:,3)))),4);
        ms_depd = round(mean(diff(unique(ms_data(:,4)))),4);
        as_lond = round(mean(diff(unique(as_data(:,2)))),4);
        as_latd = round(mean(diff(unique(as_data(:,3)))),4);
        as_depd = round(mean(diff(unique(as_data(:,4)))),4);
        if ms_lond==as_lond && ...
           ms_latd==as_latd && ...
           ms_depd==as_depd && ...
           ms_ts==as_ts
            fprintf('Grid spacing and timesteps between mainshock and reference event are consistent.\n')
        elseif ms_ts~=as_ts
            fprintf('Timesteps between mainshock and reference event are NOT consistent.\n')
            fprintf('Exiting...\n')
            return
        else
            fprintf('Grid spacing between mainshock and reference event are NOT consistent.\n')
            fprintf('Exiting...\n')
            return
        end

        %Estimate hypocenter/origin times from back-projections
        [~, ms_hyp_idx_tmp] = max(ms_data(ms_tot_pts*(ms_ot-1)+1:ms_tot_pts*ms_ot,5));
        ms_hyp_idx = ms_hyp_idx_tmp + ms_tot_pts*(ms_ot-1);
        ms_hyp = [ms_data(ms_hyp_idx,1:4) ms_ts];

        [~, as_hyp_idx_tmp] = max(as_data(as_tot_pts*(as_ot-1)+1:as_tot_pts*as_ot,5));
        as_hyp_idx = as_hyp_idx_tmp + as_tot_pts*(as_ot-1);
        as_hyp = [as_data(as_hyp_idx,1:4) as_ts];

        % ms_hyp_idxs = find(ms_data(:,5)>=ms_hyp_cont & ms_data(:,1)==ms_ot);
        % ms_hyp_tmp = mean(ms_data(ms_hyp_idxs,1:4));
        % difference = [1e10 1e10 1e10 1e10];
        % for i=1:length(ms_hyp_idxs)
        %     if any(abs(ms_data(ms_hyp_idxs(i),1:4)-ms_hyp_tmp) < difference)
        %         ms_hyp_idx=ms_hyp_idxs(i);
        %         difference=abs(ms_data(ms_hyp_idxs(i),1:4)-ms_hyp_tmp);
        %     end
        % end
        % ms_hyp = [ms_data(ms_hyp_idx,1:4) ms_ts];
        % 
        % as_hyp_idxs = find(as_data(:,5)>=as_hyp_cont & as_data(:,1)==as_ot);
        % as_hyp_tmp = mean(as_data(as_hyp_idxs,1:4));
        % difference = [1e10 1e10 1e10 1e10];
        % for i=1:length(as_hyp_idxs)
        %     if any(abs(as_data(as_hyp_idxs(i),1:4)-as_hyp_tmp) < difference)
        %         as_hyp_idx=as_hyp_idxs(i);
        %         difference=abs(as_data(as_hyp_idxs(i),1:4)-as_hyp_tmp);
        %     end
        % end
        % as_hyp = [as_data(as_hyp_idx,1:4) as_ts];

        %Convert data to be relative of mainshock hypocenter/time
        ms_data_rel = ms_data;
        ms_data_rel(:,1:4) = ms_data(:,1:4) - ms_hyp(1:4);
        as_data_rel = as_data;
        as_data_rel(:,1:4) = as_data(:,1:4) - as_hyp(1:4);

        %Contour selection
        if normamp==0
            cont_idx = find(ms_data_rel(:,5)>=cont);
        elseif normamp==1
            for t = 1:ms_t_max
                ms_data_rel_norm(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5) = ms_data_rel(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5)/max(ms_data_rel(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5));
            end
            cont_idx  = find(ms_data_rel_norm(:,5)>=cont);
        elseif normamp==222 %Multiply mainshock back-projection data by a Guassian scaled and centered at each time step, then sum/normalize across all time steps
            sigma = ms_tot_pts*ms_t_max/8;
            gauss = zeros(ms_t_max,ms_tot_pts*ms_t_max);
            gauss_ms_data_rel = zeros(ms_t_max,ms_tot_pts*ms_t_max);
            for t = 1:ms_t_max
                [gauss_amp, gauss_idx] = max(ms_data_rel(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5));
                gauss(t,:) = gauss_amp*exp(-(((1:ms_tot_pts*ms_t_max)-(gauss_idx+ms_tot_pts*(t-1))).^2)/(2*sigma^2));
                gauss_ms_data_rel(t,:) = gauss(t,:)'.*ms_data_rel(:,5);
            end
            ms_data_rel_gauss(:,5) = sum(gauss_ms_data_rel,1)/max(sum(gauss_ms_data_rel,1));
            cont_idx  = find(ms_data_rel_gauss(:,5)>=cont);
        %     f4 = figure(4);
        %     plot(1:ms_tot_pts*ms_t_max,ms_data_rel(:,5))
        %     hold on;
        %     plot(1:ms_tot_pts*ms_t_max,ms_data_rel_gauss(:,5))
        %     plot(1:ms_tot_pts*ms_t_max,ms_data_rel(:,5)-ms_data_rel_gauss(:,5)-0.5)
        %     scatter(cont_idx,ms_data_rel(cont_idx,5),50)
        %     legend('BP Amplitudes','Gaussed up BP Amplitudes','Difference (-0.5)','Points above contour')
        %     fprintf('Waiting for Figure 4 to be closed\n')
        %     uiwait(f4)
        elseif normamp==333
            cont_tbnd = [(1:ms_t_max)' zeros(ms_t_max,1)];
            cont_tbnd(ms_ot-5/ms_ts:ms_ot+15/ms_ts,2) = 1;
        %     idx = find(ms_data_rel(1:ms_tot_pts,2) >= -0.1-1e-10 & ms_data_rel(1:ms_tot_pts,2) <= 0.1+1e-10 & ...
        %                ms_data_rel(1:ms_tot_pts,3) >= -0.2-1e-10 & ms_data_rel(1:ms_tot_pts,3) <= 0.2+1e-10 & ...
        %                ms_data_rel(1:ms_tot_pts,4) >= -5.0-1e-10 & ms_data_rel(1:ms_tot_pts,4) <= 5.0+1e-10);
        %     cont_locbnd = [(1:ms_tot_pts)' zeros(ms_tot_pts,1)];
        %     cont_locbnd(idx,2) = 1;
            for t = 1:ms_t_max
                ms_data_rel_bnd(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5) = cont_locbnd(:,2).*(cont_tbnd(t,2)*ms_data_rel(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5));
            end
            cont_idx  = find(ms_data_rel_bnd(:,5));
        %     f4 = figure(4);
        %     plot(1:length(ms_data),ms_data_rel(:,5))
        %     hold on;
        %     plot(1:length(ms_data),ms_data_rel_bnd(:,5))
        %     fprintf('Waiting for Figure 4 to be closed\n')
        %     uiwait(f4)
        end

        cont_idx_t_pad_tmp = zeros(length(cont_idx)*(pad_t_2+pad_t_1+1),1);
        cnt=1;
        for t = -pad_t_1:pad_t_2
            cont_idx_t_pad_tmp(length(cont_idx)*(cnt-1)+1:length(cont_idx)*cnt) = cont_idx + ms_tot_pts*t;
            cnt=cnt+1;
        end
        cont_idx_t_pad = sort(unique(cont_idx_t_pad_tmp(cont_idx_t_pad_tmp>=1 & cont_idx_t_pad_tmp<=size(ms_data_rel,1))));

        cont_pad_tmp = zeros(length(cont_idx_t_pad)*(pad_loc_lat_2+pad_loc_lat_1+1)*(pad_loc_lon_2+pad_loc_lon_1+1)*(pad_loc_dep_2+pad_loc_dep_1+1),4);
        cnt=1;
        for k = -pad_loc_dep_1*ms_depd:ms_depd:pad_loc_dep_2*ms_depd
            for j = -pad_loc_lat_1*ms_latd:ms_latd:pad_loc_lat_2*ms_latd
                for i = -pad_loc_lon_1*ms_lond:ms_lond:pad_loc_lon_2*ms_lond
                    cont_pad_tmp(length(cont_idx_t_pad)*(cnt-1)+1:length(cont_idx_t_pad)*cnt,1) = ms_data_rel(cont_idx_t_pad,1);
                    cont_pad_tmp(length(cont_idx_t_pad)*(cnt-1)+1:length(cont_idx_t_pad)*cnt,2:4) = ms_data_rel(cont_idx_t_pad,2:4) + [i j k];
                    cnt=cnt+1;
                end
            end
        end
        cont_pad = uniquetol(cont_pad_tmp,1e-4,'ByRows',true);
        cont_pad_in = cont_pad(cont_pad(:,2) >= round(min(ms_data_rel(:,2)),4)-1e-10 & cont_pad(:,2) <= round(max(ms_data_rel(:,2)),4)+1e-10 & ...
                               cont_pad(:,3) >= round(min(ms_data_rel(:,3)),4)-1e-10 & cont_pad(:,3) <= round(max(ms_data_rel(:,3)),4)+1e-10 & ...
                               cont_pad(:,4) >= round(min(ms_data_rel(:,4)),4)-1e-10 & cont_pad(:,4) <= round(max(ms_data_rel(:,4)),4)+1e-10,:);
        cont_t    = cont_pad_in(:,1);
        cont_locs = cont_pad_in(:,2:4);

        %Test that G-matrix can be stored in memory
        if size(cont_locs,1) > 70000
            fprintf('G-matrix is too large (%i x %i), continuing...\n',size(cont_locs,1),size(cont_locs,1))
            breakinner = 1;
            continue
%         elseif size(cont_locs,1) < 10000
%             fprintf('G-matrix is too small (%i x %i), continuing...\n',size(cont_locs,1),size(cont_locs,1))
%             continue
        end
        G = zeros(size(cont_locs,1)); %Must be type double (default) for lsqnonneg
        G_idx = zeros(size(cont_locs,1),'single');
        fprintf('G-matrix will be size (%i x %i)...\n',length(G),length(G))

        %Check input parameters
        if check_in==1
            cont_data = ms_data;
            cont_data(:,5) = 0;
            fidbp_plot(ms_hyp,cont_t,cont_locs,cont_data,amp_cutoff)
            return
        end

        %Determine smallest box around ?% contour at all time steps
        %ms_rel_ts   = unique(ms_data_rel(:,1))';
        ms_rel_ts   = unique(cont_t);
        ms_rel_lons = round(min(cont_locs(:,1)),4):ms_lond:round(max(cont_locs(:,1)),4);
        ms_rel_lats = round(min(cont_locs(:,2)),4):ms_latd:round(max(cont_locs(:,2)),4);
        ms_rel_deps = round(min(cont_locs(:,3)),4):ms_depd:round(max(cont_locs(:,3)),4);

        %Calculate minimum spatiotemporal dimensions of the reference event back-projection required for deconvolution
        %Since ?% contour is over all time steps, padding in time cannot be reduced
        as_rel_pad_ts   = round(ms_rel_ts(1)-ms_rel_ts(end),4):1:round(ms_rel_ts(end)-ms_rel_ts(1),4);
        as_rel_pad_lons = round(ms_rel_lons(1)-ms_rel_lons(end),4):as_lond:round(ms_rel_lons(end)-ms_rel_lons(1),4);
        as_rel_pad_lats = round(ms_rel_lats(1)-ms_rel_lats(end),4):as_latd:round(ms_rel_lats(end)-ms_rel_lats(1),4);
        as_rel_pad_deps = round(ms_rel_deps(1)-ms_rel_deps(end),4):as_depd:round(ms_rel_deps(end)-ms_rel_deps(1),4);

        %Determine reference back-projection spatiotemporal dimensions
        as_rel_ts   = unique(as_data_rel(:,1))';
        as_rel_lons = unique(as_data_rel(:,2))';
        as_rel_lats = unique(as_data_rel(:,3))';
        as_rel_deps = unique(as_data_rel(:,4))';

        %Check if reference back-projection event data matches calculated minimum spatiotemporal dimensions exactly
        if (length(as_rel_pad_ts)   == length(as_rel_ts)   && ...
            length(as_rel_pad_lons) == length(as_rel_lons) && ...
            length(as_rel_pad_lats) == length(as_rel_lats) && ...
            length(as_rel_pad_deps) == length(as_rel_deps)) && ...
           (all(abs(as_rel_pad_ts   - as_rel_ts)   <= 1e-10) && ...
            all(abs(as_rel_pad_lons - as_rel_lons) <= 1e-10) && ...
            all(abs(as_rel_pad_lats - as_rel_lats) <= 1e-10) && ...
            all(abs(as_rel_pad_deps - as_rel_deps) <= 1e-10))
            fprintf('Reference back-projection spatiotemporal dimensions are perfectly sized for deconvolution.\n')
            as_data_rel_pad = as_data_rel;
            as_pad_hyp = as_hyp;
            as_pad_tot_pts = as_tot_pts;
        else
            as_pad_tot_pts = length(as_rel_pad_lons)*length(as_rel_pad_lats)*length(as_rel_pad_deps);
            as_pad_t_max = length(as_rel_pad_ts);
            as_data_rel_pad_tmp = zeros(length(as_rel_pad_ts)*as_pad_tot_pts,5);
            fprintf(['Reference back-projection spatiotemporal dimensions are being optimized for the %f deconvolution contour.\n', ...
                    '(Original length: %e, Optimized length: %e)\n'],...
                    cont,length(as_data_rel),length(as_data_rel_pad_tmp))
            inbounds_idx = find(as_data_rel(:,1) >= as_rel_pad_ts(1)  -1e-10 & as_data_rel(:,1) <= as_rel_pad_ts(end)  +1e-10 & ...
                                as_data_rel(:,2) >= as_rel_pad_lons(1)-1e-10 & as_data_rel(:,2) <= as_rel_pad_lons(end)+1e-10 & ...
                                as_data_rel(:,3) >= as_rel_pad_lats(1)-1e-10 & as_data_rel(:,3) <= as_rel_pad_lats(end)+1e-10 & ...
                                as_data_rel(:,4) >= as_rel_pad_deps(1)-1e-10 & as_data_rel(:,4) <= as_rel_pad_deps(end)+1e-10);
            as_data_rel_pad_tmp(1:length(inbounds_idx),:) = as_data_rel(inbounds_idx,:);
            cnt=length(inbounds_idx)+1;
            for t = as_rel_pad_ts
                for k = as_rel_pad_deps
                    for j = as_rel_pad_lats
                        for i = as_rel_pad_lons
                            if (~ismembertol(t,as_rel_ts,1e-10)    || ...
                                ~ismembertol(k,as_rel_deps,1e-10)  || ...
                                ~ismembertol(j,as_rel_lats,1e-10)  || ...
                                ~ismembertol(i,as_rel_lons,1e-10))
                                as_data_rel_pad_tmp(cnt,:) = [t i j k 0.0];
                                cnt=cnt+1;
                            end
                        end
                    end
                end
            end
            as_data_rel_pad = sortrows(as_data_rel_pad_tmp,[1 4 3 2]);
            as_data_pad = as_data_rel_pad;
            as_data_pad(:,1:4) = as_data_rel_pad(:,1:4) + as_hyp(1:4);
            as_pad_hyp = [ceil(length(as_rel_pad_ts)/2) as_hyp(2:4) as_ts];
            if as_pad_out==1
                %Write optimized reference back-projection result to text files
                fprintf('Writing optimized reference back-projection to text files at %f seconds.\n',toc)
                dlmwrite(strcat(as_pad_path,'hyp.txt'),as_pad_hyp,'delimiter','\t','precision','%7.4f')
                as_pad_amps_sum = zeros(as_pad_tot_pts,as_pad_t_max);
                for t = 1:as_pad_t_max
                    dlmwrite(strcat(as_pad_path,sprintf('3D_%04d.txt',t)), as_data_pad(as_pad_tot_pts*(t-1)+1:as_pad_tot_pts*t,2:5),'delimiter',' ','precision','%.4f')
                    as_pad_amps_sum(:,t) = as_data_pad(as_pad_tot_pts*(t-1)+1:as_pad_tot_pts*t,5);
                end
                as_pad_data_sum = cat(2,as_data_pad(1:as_pad_tot_pts,2:4),sum(as_pad_amps_sum,2)/max(sum(as_pad_amps_sum,2)));
                dlmwrite(strcat(as_pad_path,'sum.txt'), as_pad_data_sum,'delimiter',' ','precision','%.4f')
                fprintf('Optimized reference back-projection completed at %f seconds.\n',toc)
            end
        end

        %Find the indices of cont_locs in ms_data_rel for all time steps
        cont_locs_idx_tmp = zeros(size(cont_locs,1),1);
        for i = 1:size(cont_locs,1)
            ms_tstep = cont_t(i) + ms_ot;
            cont_locs_idx_tmp(i) = ms_tot_pts*(ms_tstep-1) + ...
                                   find(abs(ms_data_rel(ms_tot_pts*(ms_tstep-1)+1:ms_tot_pts*ms_tstep,1) - cont_t(i))      <= 1e-10 &...
                                        abs(ms_data_rel(ms_tot_pts*(ms_tstep-1)+1:ms_tot_pts*ms_tstep,2) - cont_locs(i,1)) <= 1e-10 & ...
                                        abs(ms_data_rel(ms_tot_pts*(ms_tstep-1)+1:ms_tot_pts*ms_tstep,3) - cont_locs(i,2)) <= 1e-10 & ...
                                        abs(ms_data_rel(ms_tot_pts*(ms_tstep-1)+1:ms_tot_pts*ms_tstep,4) - cont_locs(i,3)) <= 1e-10);
        end
        cont_locs_idx = sort(cont_locs_idx_tmp);
        ms_data_rel_cont = ms_data_rel(cont_locs_idx,:);

        %Set up inversion;
        d = ms_data_rel_cont(:,5);
        %G_idx = zeros(length(d)); %Create earlier to make sure G_idx can be stored in memory
        %G = zeros(length(d)); %Create earlier to make sure G can be stored in memroy
        fprintf('Creating G-matrix of size (%i x %i)...\n',length(G),length(G))

        %Find incices for the first row of G_idx
        i=1;
        for j = 1:size(cont_locs,1)
            as_pad_tstep = ms_data_rel_cont(i,1) - ms_data_rel_cont(j,1) + as_pad_hyp(1);
            idx = as_pad_tot_pts*(as_pad_tstep-1) + ...
                  find(abs(as_data_rel_pad(as_pad_tot_pts*(as_pad_tstep-1)+1:as_pad_tot_pts*as_pad_tstep,1) - (ms_data_rel_cont(i,1) - ms_data_rel_cont(j,1))) <= 1e-10 & ...
                       abs(as_data_rel_pad(as_pad_tot_pts*(as_pad_tstep-1)+1:as_pad_tot_pts*as_pad_tstep,2) - (ms_data_rel_cont(i,2) - ms_data_rel_cont(j,2))) <= 1e-10 & ...
                       abs(as_data_rel_pad(as_pad_tot_pts*(as_pad_tstep-1)+1:as_pad_tot_pts*as_pad_tstep,3) - (ms_data_rel_cont(i,3) - ms_data_rel_cont(j,3))) <= 1e-10 & ...
                       abs(as_data_rel_pad(as_pad_tot_pts*(as_pad_tstep-1)+1:as_pad_tot_pts*as_pad_tstep,4) - (ms_data_rel_cont(i,4) - ms_data_rel_cont(j,4))) <= 1e-10);
            G_idx(i,j) = idx;
        end
        %...then extrapolate to all rows of G_idx (knowing the symmetry of the G-matrix)
        %G(i,:) = as_data_rel_pad(G_idx(i,:),5);
        for i = 2:size(cont_locs,1)
            delta = G_idx(1,i-1)-G_idx(1,i);
            G_idx(i,:) = G_idx(i-1,:)+delta;
        %    G(i,:) = as_data_rel_pad(G_idx(i,:),5);
        end

        %Find amplitudes in G from indicies in G_idx (this is faster for some reason? -HK 11/19/2021)
        for i = 1:length(d)
            for j = 1:length(d)
                G(i,j) = as_data_rel_pad(G_idx(i,j),5);
            end
        end
        fprintf('G-matrix created at %f seconds.\n',toc)

        %Invert
        fprintf('Inverting...\n')
        [m,resnorm,residual,exitflag,output,lambda] = lsqnonneg(G,d);
        %options = optimoptions('lsqlin','Display','iter');
        %m = lsqlin(G,d,[],[],[],[],0*ones(length(d),1),1*ones(length(d),1),[],options);
        A = m/max(m);
        fprintf('Inversion complete at %f seconds.\n',toc)

        %Write inversion result to text files
        idbp_data = ms_data;
        idbp_data(:,5)=0;
        idbp_data(cont_locs_idx,5) = A;
        if save_idbp == 1
            fprintf('Writing inversion result to text files at %f seconds.\n',toc)
            dlmwrite(strcat(ms_path,'hyp.txt'),ms_hyp,'delimiter','\t','precision','%7.4f')
            dlmwrite(strcat(as_path,'hyp.txt'),as_hyp,'delimiter','\t','precision','%7.4f')
            dlmwrite(strcat(idbp_path,'hyp.txt'),ms_hyp,'delimiter','\t','precision','%7.4f')

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
                %No longer need to write all these text files -HK 10/04/2022
                %dlmwrite(strcat(idbp_path,sprintf('3D_%04d.txt',t)), idbp_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,2:5),'delimiter',' ','precision','%.4f')
                idbp_amps_sum(:,t) = idbp_data(ms_tot_pts*(t-1)+1:ms_tot_pts*t,5);
            end
            idbp_data_sum = cat(2,idbp_data(1:ms_tot_pts,2:4),sum(idbp_amps_sum,2)/max(sum(idbp_amps_sum,2)));
            dlmwrite(strcat(idbp_path,'sum.txt'), idbp_data_sum,'delimiter',' ','precision','%.4f')
            dlmwrite(strcat(idbp_path,'cont.txt'), [cont_t cont_locs],'delimiter',' ','precision','%.4f')
        end

        %Clear variables that slow performance
        clear G G_idx

        %Plot IDBP result in MATLAB
        if plot_idbp==1
            fidbp_plot(ms_hyp,cont_t,cont_locs,idbp_data,amp_cutoff)
        end

        %Save run information
        for lll = 1:length(rupv_acuts)
            [misfit,misfit2,vmisfit,vmisfit2,vangmisfit,vangmisfit2,v_loc_misfit,v_loctime_misfit,grdvol_misfit,grdvoltime_misfit,cont_misfit,conttime_misfit,~,mean_vel,vs,az,rup_dur,rup_len,num_idbp,rupv_ratio,~,~,~,~,~,~,~,~,~] = ...
            calc_misfit_rupture_prop(ms_hyp,idbp_data,rupv_acuts(lll),cont_pad_in,strike1,strike2,dip1,dip2,rake1,rake2,c_lat,c_lon,c_dep,c_dur);
            run_data(itr,:) = [itr,tdirs(l),conts(ll),rupv_acuts(lll),size(cont_locs,1),mean_vel/vs,az,rup_dur,rup_len,num_idbp,rupv_ratio,misfit,misfit2,vmisfit,vmisfit2,v_loc_misfit,v_loctime_misfit,grdvol_misfit,grdvoltime_misfit,cont_misfit,conttime_misfit,misfit/misfit2,vmisfit/vmisfit2];
            file = fopen(sprintf('%2.2fs/run_%3.2f_%3.2f.txt',tdirs(l),conts(ll),rupv_acuts(lll)),'w');
            fprintf(file, '%i %3.2f %3.2f %3.2f %i %3.2f %3.2f %3.2f %3.2f %i %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n',...
                          [itr,tdirs(l),conts(ll),rupv_acuts(lll),size(cont_locs,1),mean_vel/vs,az,rup_dur,rup_len,num_idbp,rupv_ratio,misfit,misfit2,vmisfit,vmisfit2,v_loc_misfit,v_loctime_misfit,grdvol_misfit,grdvoltime_misfit,cont_misfit,conttime_misfit,misfit/misfit2,vmisfit/vmisfit2]);
            fclose(file);
            itr = itr + 1;
        end
        
        elapsedtime = toc;
        if save_var == 1
            save(strcat(idbp_path,'variables'))
        end
    end
end

%Calculate cost and save run_data textfile
[best_cost,best_idx,cost] = cost_func(idbp_data,run_data,'run_data.txt');

%Save global variables
save('variables_global')

%Read variables for the best run
tdirs_best = run_data(best_idx,2);
conts_best = run_data(best_idx,3);
ms_path=sprintf('%2.2fs/%s_xcor/',tdirs_best,dirname);
as_path=sprintf('%2.2fs/%s_syn_xcor/',tdirs_best,dirname);
as_pad_path=sprintf('%2.2fs/Pad_%s_%3.2f/',tdirs_best,as_path(7:end-1),conts_best);
idbp_path=sprintf('%2.2fs/IDBP_%s_by_%s_%3.2f/',tdirs_best,ms_path(7:end-1),as_path(7:end-1),conts_best);
load(strcat(idbp_path,'variables.mat'))
%Write textfiles for the best run
if as_pad_out_best==1
    if exist(as_pad_path, 'dir')~=7
        mkdir(as_pad_path)
    elseif isempty(dir(strcat(as_pad_path,'*.txt')))==0 || ...
           isempty(dir(strcat(as_pad_path,'*.pdf')))==0
        fprintf(sprintf('Directory ''%s'' contains protected files.\n',as_pad_path))
        fprintf('Exiting...\n')
        return
    else
        rmdir(as_pad_path,'s')
        mkdir(as_pad_path)
    end
    %Write optimized reference back-projection result to text files
    fprintf('Writing optimized reference back-projection to text files at %f seconds.\n',toc)
    dlmwrite(strcat(as_pad_path,'hyp.txt'),as_pad_hyp,'delimiter','\t','precision','%7.4f')
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
    fprintf('Writing inversion result to text files at %f seconds.\n',toc)
    dlmwrite(strcat(ms_path,'hyp.txt'),ms_hyp,'delimiter','\t','precision','%7.4f')
    dlmwrite(strcat(as_path,'hyp.txt'),as_hyp,'delimiter','\t','precision','%7.4f')
    dlmwrite(strcat(idbp_path,'hyp.txt'),ms_hyp,'delimiter','\t','precision','%7.4f')

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

fprintf('Everything complete at %f seconds.\n',toc)
