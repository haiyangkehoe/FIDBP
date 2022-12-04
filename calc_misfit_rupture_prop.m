%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%07 September 2022
%Modified 18 November 2022

function [misfit,misfit2,vmisfit,vmisfit2,vangmisfit,vangmisfit2,v_loc_misfit,v_loctime_misfit,grdvol_misfit,grdvoltime_misfit,cont_misfit,conttime_misfit,mean_loc,mean_vel,vs,az,rup_dur,rup_len,num_idbp,rupv_ratio,data_idx_plane,rupv_idx,locs,vels,rup_len_loc,normal1,normal2,I,J] = ...
         calc_misfit_rupture_prop(hyp,data,rupv_acut,cont_pad_in,strike1,strike2,dip1,dip2,rake1,rake2,c_lat,c_lon,c_dep,c_dur)

%Determine information from inputs
[~,vs]=iasp91_lookup(hyp(4));
data_bounds = [min(data(:,1))  max(data(:,1))...
               min(data(:,2))  max(data(:,2))...
               min(data(:,3))  max(data(:,3))...
              -max(data(:,4)) -min(data(:,4))];

%Calculate GCMT normal vectors
[~,normal1]=SDR2FP(strike1,dip1,rake1);
[~,normal2]=SDR2FP(strike2,dip2,rake2);

%Calculate potential rupture velocity vectors
idx = find(data(:,5)>=rupv_acut);
data_idx = sortrows(sortrows(data(idx,:),5,'descend'),1);
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

%Only keep rupture velocities and IDBP centroids under P-wave speeds
[vp,~]=iasp91_lookup(hyp(4));
rupv_idx = find(vels>0 & vels<vp);
rupv_idx_all = find(vels>0);
if ~isempty(rupv_idx)
    %mean_loc = mean(locs(rupv_idx,:),1);
    mean_loc = [sum(mean(amps(rupv_idx,:),2).*locs(rupv_idx,1:3),1)/sum(mean(amps(rupv_idx,:),2)) sum(mean(amps(rupv_idx,:),2).*locs(rupv_idx,4:6),1)/sum(mean(amps(rupv_idx,:),2))];
    %mean_tstep = mean(tsteps(rupv_idx,:),1);
    mean_tstep = [sum(mean(amps(rupv_idx,:),2).*tsteps(rupv_idx,1),1)/sum(mean(amps(rupv_idx,:),2)) sum(mean(amps(rupv_idx,:),2).*tsteps(rupv_idx,2),1)/sum(mean(amps(rupv_idx,:),2))];
    %mean_vel = mean(vels(rupv_idx,:));
    mean_vel = sum(mean(amps(rupv_idx,:),2).*vels(rupv_idx))/sum(mean(amps(rupv_idx,:),2));
    [~,az] = distance(mean_loc(2),mean_loc(1),mean_loc(5),mean_loc(4));
    rupv_ratio = size(rupv_idx,1)/size(rupv_idx_all,1);

    %Fit GCMT Solution
    data_idx_plane = unique([tsteps(rupv_idx,1) locs(rupv_idx,1:3) amps(rupv_idx,1); tsteps(rupv_idx,2) locs(rupv_idx,4:6) amps(rupv_idx,2)],'rows');
    rup_dur = (max(data_idx_plane(:,1))-min(data_idx_plane(:,1)))*hyp(5);
    data_idx_plane_dists = zeros(size(data_idx_plane,1)^2,1);
    rup_len_idxs = zeros(size(data_idx_plane,1)^2,2);
    cnt=1;
    for i = 1:size(data_idx_plane,1)
        for j = 1:size(data_idx_plane,1)
            R = 6371-((data_idx_plane(j,4)+data_idx_plane(i,4))/2);
            rup_len_idxs(cnt,:) = [i j];
            data_idx_plane_dists(cnt,1) = sqrt(deg2km(data_idx_plane(j,2)-data_idx_plane(i,2),R).^2 + ...
                                          deg2km(data_idx_plane(j,3)-data_idx_plane(i,3),R).^2 + ...
                                          (data_idx_plane(j,4)-data_idx_plane(i,4)).^2);
            cnt=cnt+1;
        end
    end
    [rup_len,max_rup_len_idx] = max(data_idx_plane_dists);
    rup_len_loc = [data_idx_plane(rup_len_idxs(max_rup_len_idx,1),2:4) data_idx_plane(rup_len_idxs(max_rup_len_idx,2),2:4)];
    %num_idbp = size(data_idx_plane((data_idx_plane(:,5)>0.75),:),1); %Number of IDBP points above cutoff value
    num_idbp = size(data_idx_plane(:,5)>rupv_acut,1);
    num_idbp_ratio = num_idbp/length(idx);
    dist1 = zeros(size(data_idx_plane,1),1);
    dist2 = zeros(size(data_idx_plane,1),1);
    for i = 1:size(data_idx_plane,1)
        dist1(i) = sum(abs(normal1(1)*deg2km( data_idx_plane(:,2) - data_idx_plane(i,2)) + ...
                           normal1(2)*deg2km( data_idx_plane(:,3) - data_idx_plane(i,3)) + ...
                           normal1(3)*      (-data_idx_plane(:,4) + data_idx_plane(i,4)))/ ...
                           sqrt(normal1(1)^2 + normal1(2)^2 + normal1(3)^2))/num_idbp;
        dist2(i) = sum(abs(normal2(1)*deg2km( data_idx_plane(:,2) - data_idx_plane(i,2)) + ...
                           normal2(2)*deg2km( data_idx_plane(:,3) - data_idx_plane(i,3)) + ...
                           normal2(3)*      (-data_idx_plane(:,4) + data_idx_plane(i,4)))/ ...
                           sqrt(normal2(1)^2 + normal2(2)^2 + normal2(3)^2))/num_idbp;
    end
    dist12 = [dist1 dist2];
    [~,I] = min(dist12);
    [misfit,J] = min([min(dist12(I,1)) min(dist12(I,2))]);
    best_idx = I(J);
    [misfit2,~] = max([min(dist12(I,1)) min(dist12(I,2))]);
    
    %Find misfit of rupture velocity vector from best-fitting and worst-fitting planes
    if J == 1
        vmisfit = sum(abs(normal1(1)*deg2km( [mean_loc(1); mean_loc(4)] - data_idx_plane(best_idx,2)) + ...
                          normal1(2)*deg2km( [mean_loc(2); mean_loc(5)] - data_idx_plane(best_idx,3)) + ...
                          normal1(3)*      (-[mean_loc(3); mean_loc(6)] + data_idx_plane(best_idx,4)))/ ...
                          sqrt(normal1(1)^2 + normal1(2)^2 + normal1(3)^2))/2;
        vmisfit2 = sum(abs(normal2(1)*deg2km( [mean_loc(1); mean_loc(4)] - data_idx_plane(best_idx,2)) + ...
                           normal2(2)*deg2km( [mean_loc(2); mean_loc(5)] - data_idx_plane(best_idx,3)) + ...
                           normal2(3)*      (-[mean_loc(3); mean_loc(6)] + data_idx_plane(best_idx,4)))/ ...
                           sqrt(normal2(1)^2 + normal2(2)^2 + normal2(3)^2))/2;
        %Find misfit angle of rupture velocity vector from best-fitting and worst-fitting planes
        %If vangmisfit/vangmisfit2 > 1, the best-fitting plane has worse-fitting rupture vector angle
        R = 6371-mean_loc(3);
        rupv = [deg2km(mean_loc(4)-mean_loc(1),R);
                deg2km(mean_loc(5)-mean_loc(2),R);
                      -mean_loc(6)+mean_loc(3)];
        vangmisfit  = abs(90 - acosd(dot(normal1,rupv)/(norm(normal2)*norm(rupv))));
        vangmisfit2 = abs(90 - acosd(dot(normal2,rupv)/(norm(normal1)*norm(rupv))));
    elseif J == 2
        vmisfit = sum(abs(normal2(1)*deg2km( [mean_loc(1); mean_loc(4)] - data_idx_plane(best_idx,2)) + ...
                          normal2(2)*deg2km( [mean_loc(2); mean_loc(5)] - data_idx_plane(best_idx,3)) + ...
                          normal2(3)*      (-[mean_loc(3); mean_loc(6)] + data_idx_plane(best_idx,4)))/ ...
                          sqrt(normal2(1)^2 + normal2(2)^2 + normal2(3)^2))/2;
        vmisfit2 = sum(abs(normal1(1)*deg2km( [mean_loc(1); mean_loc(4)] - data_idx_plane(best_idx,2)) + ...
                           normal1(2)*deg2km( [mean_loc(2); mean_loc(5)] - data_idx_plane(best_idx,3)) + ...
                           normal1(3)*      (-[mean_loc(3); mean_loc(6)] + data_idx_plane(best_idx,4)))/ ...
                           sqrt(normal1(1)^2 + normal1(2)^2 + normal1(3)^2))/2;
        %Find misfit angle of rupture velocity vector from best-fitting and worst-fitting planes
        %If vangmisfit/vangmisfit2 > 1, the best-fitting plane has worse-fitting rupture vector angle
        R = 6371-mean_loc(3);
        rupv = [deg2km(mean_loc(4)-mean_loc(1),R);
                deg2km(mean_loc(5)-mean_loc(2),R);
                      -mean_loc(6)+mean_loc(3)];
        vangmisfit  = abs(90 - acosd(dot(normal2,rupv)/(norm(normal2)*norm(rupv))));
        vangmisfit2 = abs(90 - acosd(dot(normal1,rupv)/(norm(normal1)*norm(rupv))));
    end
    
    %Find misfit of the IDBP radiators from an infinite line defined by the rupture velocity vector
    %(https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html)
    v_loc_misfits = zeros(size(data_idx_plane,1),1);
    proj_data_idx_plane = zeros(size(data_idx_plane,1),3);
    proj_dist_ratio = zeros(size(data_idx_plane,1),1);
    proj_mean_tstep = zeros(size(data_idx_plane,1),1);
    v_loctime_misfits = zeros(size(data_idx_plane,1),1);
    for i = 1:size(data_idx_plane,1)
        R = 6371-data_idx_plane(i,4);
        v_loc_misfits(i) = norm(cross([deg2km(data_idx_plane(i,2)-mean_loc(1),R) deg2km(data_idx_plane(i,3)-mean_loc(2),R) (-data_idx_plane(i,4)+mean_loc(3))],...
                                      [deg2km(data_idx_plane(i,2)-mean_loc(4),R) deg2km(data_idx_plane(i,3)-mean_loc(5),R) (-data_idx_plane(i,4)+mean_loc(6))]))/...
                                 norm([deg2km(mean_loc(4)-mean_loc(1),R) deg2km(mean_loc(5)-mean_loc(2),R) (-mean_loc(6) + mean_loc(3))]);
        %Find coordinates of the IDBP radiator projected onto an infinite line defined by the rupture velocity vector
        vec_a = [deg2km(data_idx_plane(i,2)-mean_loc(1),R) deg2km(data_idx_plane(i,3)-mean_loc(2),R) -data_idx_plane(i,4)+mean_loc(3)];
        vec_b = [deg2km(mean_loc(4)-mean_loc(1),R) deg2km(mean_loc(5)-mean_loc(2),R) -mean_loc(6)+mean_loc(3)];
        proj = dot(vec_a,vec_b)/dot(vec_b,vec_b)*vec_b;
        proj_data_idx_plane(i,:) = [mean_loc(1) mean_loc(2) mean_loc(3)] + [km2deg(proj(1),R) km2deg(proj(2),R) -proj(3)];
        %Find time each IDBP should occur given they all occur on an infinite line defined by the rupture velocity vector
        proj_dist_ratio(i) = dot(vec_a,vec_b)/dot(vec_b,vec_b);
        proj_mean_tstep(i) = mean_tstep(1) + proj_dist_ratio(i)*(mean_tstep(2)-mean_tstep(1));
        %Find misfit between projected time and actual time
        v_loctime_misfits(i) = abs(data_idx_plane(i,1)-proj_mean_tstep(i))*hyp(5);
    end
    %Debugging figure
%     figure(11)
%     hold on;
%     colormap(haxby_hk);
%     colorbar('Position',[.92 .11 .03 .815]);
%     caxis([0 1])
%     scatter3(proj_data_idx_plane(:,1),proj_data_idx_plane(:,2),-proj_data_idx_plane(:,3),100,'filled','b')
%     scatter3(data_idx_plane(:,2),...
%              data_idx_plane(:,3),...
%             -data_idx_plane(:,4),...
%              data_idx_plane(:,5)*100,...
%              data_idx_plane(:,5),...
%              'filled',...
%              'MarkerEdgeColor','k');
%     quiver3(mean_loc(1),...
%             mean_loc(2),...
%            -mean_loc(3),...
%             mean_loc(4)-mean_loc(1),...
%             mean_loc(5)-mean_loc(2),...
%            -mean_loc(6)+mean_loc(3),...
%             'AutoScale','off','LineWidth',4,'Color','k');
%     view(2)
    v_loc_misfit = sum(v_loc_misfits)/size(data_idx_plane,1);
	v_loctime_misfit = sum(v_loctime_misfits)/size(data_idx_plane,1);
    
    %Find distance from each IDBP point to the center of BP grid
    mean_data_loc = mean(data(:,2:4));
    R = 6371-mean_data_loc(3);
    max_grdvol_misfit = max(sqrt(deg2km(mean_data_loc(1)-data(:,2),R).^2 + ...
                                 deg2km(mean_data_loc(2)-data(:,3),R).^2 + ...
                                       (mean_data_loc(3)-data(:,4)).^2));
    mean_data_time = mean(data(:,1));
    max_grdvoltime_misfit = max(abs(mean_data_time-data(:,1))*hyp(5));
    grdvol_misfits = zeros(size(data_idx_plane,1),1);
    grdvoltime_misfits = zeros(size(data_idx_plane,1),1);
    for i = 1:size(data_idx_plane,1)
        R = 6371-((mean_data_loc(3)+data_idx_plane(i,4))/2);
        grdvol_misfits(i) = sum(sqrt(deg2km(mean_data_loc(1)-data_idx_plane(i,2),R).^2 + ...
                                     deg2km(mean_data_loc(2)-data_idx_plane(i,3),R).^2 + ...
                                           (mean_data_loc(3)-data_idx_plane(i,4)).^2));
        grdvoltime_misfits(i) = abs(mean_data_time-data_idx_plane(i,1))*hyp(5);
        %Debugging figure
%         figure(i+20)
%         plot3(mean_data_loc(1),mean_data_loc(2),-mean_data_loc(3),'g','MarkerSize',10,'Marker','+')
%         hold on;
%         scatter3(data_idx_plane(i,2),...
%                  data_idx_plane(i,3),...
%                 -data_idx_plane(i,4),...
%                  data_idx_plane(i,5)*100,...
%                  data_idx_plane(i,5),...
%                  'filled','MarkerEdgeColor','k');
%         axis(data_bounds(3:8))
%         daspect([1 1 100])
%         view(2)
%         colormap(haxby_hk);
%         caxis([0 1])
%         colorbar('Position',[.92 .11 .03 .815]);
    end
    grdvol_misfit = sum(grdvol_misfits)/(max_grdvol_misfit*size(data_idx_plane,1));
    grdvoltime_misfit = sum(grdvoltime_misfits)/(max_grdvoltime_misfit*size(data_idx_plane,1));
    
    %Find distance from each IDBP point to the center of the IDBP contour at each time step
    mean_conttime_pad_in = mean(unique(cont_pad_in(:,1)+hyp(1)));
    max_conttime_misfit = max(abs(mean_conttime_pad_in-unique(cont_pad_in(:,1)+hyp(1)))*hyp(5));
    cont_misfits = zeros(size(data_idx_plane,1),1);
    max_cont_misfits = zeros(size(data_idx_plane,1),1);
    conttime_misfits  = zeros(size(data_idx_plane,1),1);
    num_cont_pad_in_tstep = zeros(size(data_idx_plane,1),1);
    for i = 1:size(data_idx_plane,1)
        cont_pad_in_tstep = cont_pad_in(cont_pad_in(:,1)==data_idx_plane(i,1)-hyp(1),:);
        shp = alphaShape(cont_pad_in_tstep(:,2) + hyp(2),cont_pad_in_tstep(:,3) + hyp(3),-cont_pad_in_tstep(:,4) - hyp(4));
        [~,ID1] = inShape(shp,data_idx_plane(i,2),data_idx_plane(i,3),-data_idx_plane(i,4));
        if isnan(ID1)
            num_cont_pad_in_tstep(i) = 0;
            cont_misfits(i) = inf;
        else
            itr=1;
            cont_pad_in_tstep_region = [];
            for j = 1:size(cont_pad_in_tstep,1)
                [~,ID2] = inShape(shp,cont_pad_in_tstep(j,2) + hyp(2),cont_pad_in_tstep(j,3) + hyp(3),-cont_pad_in_tstep(j,4) - hyp(4));
                if ID1==ID2
                    cont_pad_in_tstep_region(itr,:) = cont_pad_in_tstep(j,:);
                    itr=itr+1;
                end
            end
            mean_cont_pad_in = mean(cont_pad_in_tstep_region(:,2:4)) + hyp(2:4);
            R = 6371-mean_cont_pad_in(3);
            max_cont_misfits(i) = max(sqrt(deg2km(mean_cont_pad_in(1)-cont_pad_in_tstep_region(:,2)-hyp(2),R).^2 + ...
                                           deg2km(mean_cont_pad_in(2)-cont_pad_in_tstep_region(:,3)-hyp(3),R).^2 + ...
                                                 (mean_cont_pad_in(3)-cont_pad_in_tstep_region(:,4)-hyp(4)).^2));
            num_cont_pad_in_tstep(i) = size(cont_pad_in_tstep_region,1);
            R = 6371-((mean_cont_pad_in(3)+data_idx_plane(i,4))/2);
            cont_misfits(i) = sum(sqrt(deg2km(mean_cont_pad_in(1)-data_idx_plane(i,2),R).^2 + ...
                                       deg2km(mean_cont_pad_in(2)-data_idx_plane(i,3),R).^2 + ...
                                             (mean_cont_pad_in(3)-data_idx_plane(i,4)).^2));
            conttime_misfits(i) = abs(mean_conttime_pad_in-data_idx_plane(i,1))*hyp(5);
            %Debugging figure
%             figure(i+10)
%             plot(shp,'FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2)
%             hold on;
%             scatter3(data_idx_plane(i,2),...
%                      data_idx_plane(i,3),...
%                     -data_idx_plane(i,4),...
%                      data_idx_plane(i,5)*100,...
%                      data_idx_plane(i,5),...
%                      'filled','MarkerEdgeColor','k');
%             axis(data_bounds(3:8))
%             daspect([1 1 100])
%             view(2)
%             colormap(haxby_hk);
%             caxis([0 1])
%             colorbar('Position',[.92 .11 .03 .815]);
%             scatter3(cont_pad_in_tstep_region(:,2)+hyp(2),cont_pad_in_tstep_region(:,3)+hyp(3),-cont_pad_in_tstep_region(:,4)-hyp(4),'r','Marker','+')
%             plot3(mean_cont_pad_in(1),mean_cont_pad_in(2),-mean_cont_pad_in(3),'g','MarkerSize',10,'Marker','+')        
        end
    end
    cont_misfit = sum(cont_misfits)/size(data_idx_plane,1);
    conttime_misfit = sum(conttime_misfits)/size(data_idx_plane,1);
else
    misfit = NaN;
    misfit2 = NaN;
    vmisfit = NaN;
    vmisfit2 = NaN;
    vangmisfit = NaN;
    vangmisfit2 = NaN;
    mean_loc = [NaN NaN NaN NaN NaN NaN];
    mean_vel = NaN;
    az = NaN;
    rup_dur = NaN;
    rup_len = NaN;
    num_idbp = NaN;
    %num_idbp_ratio = NaN;
    rupv_ratio = NaN;
    data_idx_plane = [NaN NaN NaN NaN NaN NaN];
    rupv_idx = NaN;
    rup_len_loc = NaN;
    v_loc_misfit = NaN;
    v_loctime_misfit = NaN;
    grdvol_misfit = NaN;
    grdvoltime_misfit = NaN;
    cont_misfit = NaN;
    conttime_misfit = NaN;
    I = [];
    J = [];
end