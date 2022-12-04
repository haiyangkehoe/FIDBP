%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%17 October 2022
%Modified 17 November 2022

%estimate_smearing.m estimates the smearing direction of a point-source 
%back-projection result

clear;
close all;

%User inputs
run_data = dlmread('run_data.txt');
dirname = 'Aftershock_0';
amp_cut = 0; %Minimum back-projection amplitude to estimate smearing
save_param = 1;
plot_result = 0;

idx = dlmread('best_idxs.txt','\t',[0 1 0 1]);
tdirs_best = run_data(idx,2);
as_path = sprintf('%2.2fs/%s_syn_xcor/',tdirs_best,dirname);

%Find back-projection dimensions
as_tot_pts = length(dlmread(strcat(as_path,'3D_0001.txt')));
as_t_max = length(dir([as_path '3D_*.txt']));
hyp = dlmread(strcat(as_path,'hyp.txt'));

%Load data
%load('variables_global.mat','as_data','as_tot_pts','as_t_max');

%Read back-projection data
as_data = zeros(as_tot_pts*as_t_max,5);
as_data_max_locs = nan(as_t_max,5);
for t = 1:as_t_max
    as_data(as_tot_pts*(t-1)+1:as_tot_pts*t,:) = [t*ones(as_tot_pts,1) dlmread(strcat(as_path,sprintf('3D_%04d.txt',t)),'')];
    
    %Calculate mean location above amplitude cutoff at each time step
    [~,idx_tmp] = max(as_data(as_data(:,1)==t,5));
    idx = as_tot_pts*(t-1) + idx_tmp;
    as_data_max_locs(t,:) = as_data(idx,:);
end

%Determine information from inputs
data_bounds = [min(as_data(:,1))  max(as_data(:,1))...
               min(as_data(:,2))  max(as_data(:,2))...
               min(as_data(:,3))  max(as_data(:,3))...
              -max(as_data(:,4)) -min(as_data(:,4))];

%Calculate potential rupture velocity vectors
idx = find(as_data_max_locs(:,5)>amp_cut);
data_idx = sortrows(sortrows(as_data_max_locs(idx,:),5,'descend'),1);
idxs  = zeros(length(idx)^2,2);
locs  = zeros(length(idx)^2,6);
dists = zeros(length(idx)^2,1);
tsteps = zeros(length(idx)^2,2);
times = zeros(length(idx)^2,1);
amps  = zeros(length(idx)^2,2);
cnt=1;
for i = 1:length(idx)
    for j = 1:length(idx)
        idxs(cnt,:)  = [idx(i) idx(j)];
        locs(cnt,:)  = [data_idx(i,2:4) data_idx(j,2:4)];
        R = 6371-((data_idx(j,4)+data_idx(i,4))/2); %Change other instances of deg2km and distance to the correct radius. -HK 09/09/2022
        dists(cnt,:) = sqrt(deg2km(data_idx(j,2)-data_idx(i,2),R).^2 + ...
                            deg2km(data_idx(j,3)-data_idx(i,3),R).^2 + ...
                           (data_idx(j,4)-data_idx(i,4)).^2);
        tsteps(cnt,:) = [data_idx(i,1) data_idx(j,1)];
        times(cnt,:) = hyp(5)*((1-hyp(1))+(data_idx(j,1))-1) - hyp(5)*((1-hyp(1))+(data_idx(i,1))-1);
        amps(cnt,:)  = [data_idx(i,5) data_idx(j,5)];
        cnt=cnt+1;
    end
end
vels = dists./times;

[vp,vs]=iasp91_lookup(hyp(4));
rupv_idx = find(vels>0);
%[~,ia,~] = unique(idxs(rupv_idx_tmp,2),'first');
%rupv_idx = rupv_idx_tmp(ia);
%mean_loc = mean(locs(rupv_idx,:),1);
%Weight vector locations by amplitude
mean_loc = [sum(amps(rupv_idx,1).*locs(rupv_idx,1:3),1)/sum(amps(rupv_idx,1)) sum(amps(rupv_idx,2).*locs(rupv_idx,4:6),1)/sum(amps(rupv_idx,2))];
mean_vel = mean(vels(rupv_idx,:));
[~,az] = distance(mean_loc(2),mean_loc(1),mean_loc(5),mean_loc(4));
    figure(1)
    ax1 = axes;
    colormap(ax1,haxby_hk);
    colorbar(ax1,'Position',[.92 .11 .03 .815]);
    caxis(ax1,[0 1])
    daspect(ax1,[1 1 100])
    hold on;
    scatter3(ax1,hyp(2),hyp(3),-hyp(4),250,'filled','r','p');
    scatter3(ax1,as_data_max_locs(:,2),...
                 as_data_max_locs(:,3),...
                -as_data_max_locs(:,4),...
                 as_data_max_locs(:,5)*100,...
                 as_data_max_locs(:,5),...
                 'filled',...
                 'MarkerEdgeColor','k');
    % axis(ax1,data_bounds(3:8))
    % for i = 1:length(rupv_idx)
    %     h1 = quiver3(ax1,locs(rupv_idx(i),1),...
    %                      locs(rupv_idx(i),2),...
    %                     -locs(rupv_idx(i),3),...
    %                      locs(rupv_idx(i),4)-locs(rupv_idx(i),1),...
    %                      locs(rupv_idx(i),5)-locs(rupv_idx(i),2),...
    %                     -locs(rupv_idx(i),6)+locs(rupv_idx(i),3),...
    %                      'AutoScale','off','LineWidth',2);
    %     h1.Color = [0.7 0.7 0.7];
    %     h1.MaxHeadSize = 1/norm([h1.UData h1.VData h1.WData]);
    %     text(ax1,locs(rupv_idx(i),1)+((locs(rupv_idx(i),4)-locs(rupv_idx(i),1))/2),...
    %              locs(rupv_idx(i),2)+((locs(rupv_idx(i),5)-locs(rupv_idx(i),2))/2),...
    %             -locs(rupv_idx(i),3)+((-locs(rupv_idx(i),6)+locs(rupv_idx(i),3))/2),...
    %              sprintf('%4.2f km/s; %1.2fVp',vels(rupv_idx(i),:),vels(rupv_idx(i),:)/vp),'Color',h1.Color)
    % end
    h2 = quiver3(ax1,mean_loc(1),...
                     mean_loc(2),...
                    -mean_loc(3),...
                     mean_loc(4)-mean_loc(1),...
                     mean_loc(5)-mean_loc(2),...
                    -mean_loc(6)+mean_loc(3),...
                    'AutoScale','off','LineWidth',4);
    h2.MaxHeadSize = 0.2/norm([h2.UData h2.VData h2.WData]);
    h2.Color = 'k';
    text(ax1,mean_loc(1)+((mean_loc(4)-mean_loc(1))/2),...
             mean_loc(2)+((mean_loc(5)-mean_loc(2))/2),...
            -mean_loc(3)+((-mean_loc(6)+mean_loc(3))/2),...
             sprintf('%4.2f km/s; %1.2fVp',mean_vel,mean_vel/vp),'Color',h2.Color)
    xlabel(ax1,'Lon (°)')
    ylabel(ax1,'Lat (°)')
    zlabel(ax1,'Depth (km)')
%Plot static figure
if plot_result == 1

end

if save_param == 1
    R = 6371-mean_loc(3);
    rupv = [deg2km(mean_loc(1)-mean_loc(4),R);
            deg2km(mean_loc(2)-mean_loc(5),R);
                   mean_loc(3)-mean_loc(6)];
    horiz = [0;0;-1];
    rupv_ang_from_horiz = acosd(dot(rupv,horiz)/(norm(rupv)*norm(horiz))) - 90;
    file = fopen('parameters_smearing.txt', 'w');
    fprintf(file, '%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %3.2f %3.2f %3.2f %3.2f %3.2f\n',...
                  [mean_loc...
                   mean_vel...
                   mean_vel/vs...
                   az...
                   rupv_ang_from_horiz...
                   mean_vel/vp]);
    fclose(file);
    [~,closest_idx] = min(sum(abs(as_data(1:as_tot_pts,2:4)-mean(data_idx(:,2:4))),2));
    smearing_cent = as_data(closest_idx,2:4);
    dlmwrite(strcat(as_path,'smearing_cent.txt'),smearing_cent,'delimiter','\t','precision',7);
end