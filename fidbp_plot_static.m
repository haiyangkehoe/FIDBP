%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%24 June 2022
%Modified 10 November 2022

function [parameters] = fidbp_plot_static(hyp,data,rupv_acut,cont_pad_in,strike1,strike2,dip1,dip2,rake1,rake2,c_lat,c_lon,c_dep,c_dur,misfit_max,misfit2_max,misfit_ratio_max,vmisfit_max,vmisfit2_max,vmisfit_ratio_max)

%Oragnize nodal planes
np = [strike1 strike2; dip1 dip2; rake1 rake2];

%Determine information from inputs
data_bounds = [min(data(:,1))  max(data(:,1))...
               min(data(:,2))  max(data(:,2))...
               min(data(:,3))  max(data(:,3))...
              -max(data(:,4)) -min(data(:,4))];
[xx,yy] = meshgrid(data_bounds(3):0.05:data_bounds(4),data_bounds(5):0.05:data_bounds(6));
[vp,vs]=iasp91_lookup(hyp(4));
[misfit,misfit2,vmisfit,vmisfit2,vangmisfit,vangmisfit2,~,~,~,~,~,~,mean_loc,mean_vel,vs,az,rup_dur,rup_len,~,~,data_idx_plane,rupv_idx,locs,vels,rup_len_loc,normal1,normal2,I,J] = ...
calc_misfit_rupture_prop(hyp,data,rupv_acut,cont_pad_in,strike1,strike2,dip1,dip2,rake1,rake2,c_lat,c_lon,c_dep,c_dur);

%Classify rupture with respect to GCMT nodal planes (fits single plane, both planes, or neither planes)
if vangmisfit2 - vangmisfit > 10 && vangmisfit/vangmisfit2 <= 0.5
    fprintf('Single nodal plane selected.\n')
%     fprintf('misfit = %3.2f\n',misfit)
%     fprintf('misfit_max = %3.2f\n',misfit_max)
%     fprintf('misfit2 = %3.2f\n',misfit2)
%     fprintf('misfit2_max = %3.2f\n',misfit2_max)
%     fprintf('misfit/misfit2 = %3.2f\n',misfit/misfit2)
%     fprintf('misfit_ratio_max = %3.2f\n',misfit_ratio_max)
%     fprintf('vmisfit = %3.2f\n',vmisfit)
%     fprintf('vmisfit_max = %3.2f\n',vmisfit_max)
%     fprintf('vmisfit2 = %3.2f\n',vmisfit2)
%     fprintf('vmisfit2_max = %3.2f\n',vmisfit2_max)
%     fprintf('vmisfit/vmisfit2 = %3.2f\n',vmisfit/vmisfit2)
%     fprintf('vmisfit_ratio_max = %3.2f\n',vmisfit_ratio_max)
    fprintf('vangmisfit = %3.2f\n',vangmisfit)
    fprintf('vangmisfit2 = %3.2f\n',vangmisfit2)
    fprintf('vangmisfit_ratio = %3.2f\n',vangmisfit/vangmisfit2)
    pref_strike  = np(1,J);
    pref_dip     = np(2,J);
    pref_rake    = np(3,J);
    other_strike = NaN;
    other_dip    = NaN;
    other_rake   = NaN;
elseif vangmisfit2 - vangmisfit <= 10 && vangmisfit/vangmisfit2 <= 0.5
    fprintf('Both nodal planes selected.\n')
%     fprintf('misfit = %3.2f\n',misfit)
%     fprintf('misfit_max = %3.2f\n',misfit_max)
%     fprintf('misfit2 = %3.2f\n',misfit2)
%     fprintf('misfit2_max = %3.2f\n',misfit2_max)
%     fprintf('misfit/misfit2 = %3.2f\n',misfit/misfit2)
%     fprintf('misfit_ratio_max = %3.2f\n',misfit_ratio_max)
%     fprintf('vmisfit = %3.2f\n',vmisfit)
%     fprintf('vmisfit_max = %3.2f\n',vmisfit_max)
%     fprintf('vmisfit2 = %3.2f\n',vmisfit2)
%     fprintf('vmisfit2_max = %3.2f\n',vmisfit2_max)
%     fprintf('vmisfit/vmisfit2 = %3.2f\n',vmisfit/vmisfit2)
%     fprintf('vmisfit_ratio_max = %3.2f\n',vmisfit_ratio_max)
    fprintf('vangmisfit = %3.2f\n',vangmisfit)
    fprintf('vangmisfit2 = %3.2f\n',vangmisfit2)
    fprintf('vangmisfit_ratio = %3.2f\n',vangmisfit/vangmisfit2)
    pref_strike  = np(1,J);
    pref_dip     = np(2,J);
    pref_rake    = np(3,J);
    if J == 1
        other_strike = np(1,2);
        other_dip    = np(2,2);
        other_rake   = np(3,2);
    elseif J == 2
        other_strike = np(1,1);
        other_dip    = np(2,1);
        other_rake   = np(3,1);
    end
else
    fprintf('Neither nodal plane selected.\n')
%     fprintf('misfit = %3.2f\n',misfit)
%     fprintf('misfit_max = %3.2f\n',misfit_max)
%     fprintf('misfit2 = %3.2f\n',misfit2)
%     fprintf('misfit2_max = %3.2f\n',misfit2_max)
%     fprintf('misfit/misfit2 = %3.2f\n',misfit/misfit2)
%     fprintf('misfit_ratio_max = %3.2f\n',misfit_ratio_max)
%     fprintf('vmisfit = %3.2f\n',vmisfit)
%     fprintf('vmisfit_max = %3.2f\n',vmisfit_max)
%     fprintf('vmisfit2 = %3.2f\n',vmisfit2)
%     fprintf('vmisfit2_max = %3.2f\n',vmisfit2_max)
%     fprintf('vmisfit/vmisfit2 = %3.2f\n',vmisfit/vmisfit2)
%     fprintf('vmisfit_ratio_max = %3.2f\n',vmisfit_ratio_max)
    fprintf('vangmisfit = %3.2f\n',vangmisfit)
    fprintf('vangmisfit2 = %3.2f\n',vangmisfit2)
    fprintf('vangmisfit_ratio = %3.2f\n',vangmisfit/vangmisfit2)
    pref_strike  = NaN;
    pref_dip     = NaN;
    pref_rake    = NaN;
    other_strike = NaN;
    other_dip    = NaN;
    other_rake   = NaN;
end
plane_cent = data_idx_plane(I(J),:);
zz =  normal1(1)/normal1(3)*(deg2km(plane_cent(2))-deg2km(xx)) + normal1(2)/normal1(3)*(deg2km(plane_cent(3))-deg2km(yy)) - plane_cent(4);
zzz = normal2(1)/normal2(3)*(deg2km(plane_cent(2))-deg2km(xx)) + normal2(2)/normal2(3)*(deg2km(plane_cent(3))-deg2km(yy)) - plane_cent(4);

%Get Slab2 parameters
[s_loc,s_dip,s_str,s_thk] = slab2_lookup(mean_loc(2),mean_loc(1),mean_loc(3));

%Plot static figure
figure(1)
ax1 = axes;
colormap(ax1,haxby_hk);
colorbar(ax1,'Position',[.92 .11 .03 .815]);
caxis(ax1,[0 1])
daspect(ax1,[1 1 100])
hold on;
surf(ax1,xx,yy,zz,'facecolor','r')
surf(ax1,xx,yy,zzz,'facecolor','b')
alpha 0.1
scatter3(ax1,hyp(2),hyp(3),-hyp(4),250,'filled','r','p');
%scatter3(ax1,s_loc(1),s_loc(2),-s_loc(3),250,'filled','b','h');
scatter3(ax1,data_idx_plane(:,2),...
             data_idx_plane(:,3),...
            -data_idx_plane(:,4),...
             data_idx_plane(:,5)*100,...
             data_idx_plane(:,5),...
             'filled',...
             'MarkerEdgeColor','k');
for i = 1:length(rupv_idx)
    h1 = quiver3(ax1,locs(rupv_idx(i),1),...
                     locs(rupv_idx(i),2),...
                    -locs(rupv_idx(i),3),...
                     locs(rupv_idx(i),4)-locs(rupv_idx(i),1),...
                     locs(rupv_idx(i),5)-locs(rupv_idx(i),2),...
                    -locs(rupv_idx(i),6)+locs(rupv_idx(i),3),...
                     'AutoScale','off','LineWidth',2);
    h1.Color = [0.7 0.7 0.7];
    h1.MaxHeadSize = 1/norm([h1.UData h1.VData h1.WData]);
    text(ax1,locs(rupv_idx(i),1)+((locs(rupv_idx(i),4)-locs(rupv_idx(i),1))/2),...
             locs(rupv_idx(i),2)+((locs(rupv_idx(i),5)-locs(rupv_idx(i),2))/2),...
            -locs(rupv_idx(i),3)+((-locs(rupv_idx(i),6)+locs(rupv_idx(i),3))/2),...
             sprintf('%4.2f km/s; %1.2fVs',vels(rupv_idx(i),:),vels(rupv_idx(i),:)/vs),'Color',h1.Color)
end
plot3(ax1,[rup_len_loc(1); rup_len_loc(4);],[rup_len_loc(2); rup_len_loc(5);],-[rup_len_loc(3); rup_len_loc(6);],...
          '--','Color',[0 0 0],'linewidth',3)
h2 = quiver3(ax1,mean_loc(1),...
                 mean_loc(2),...
                -mean_loc(3),...
                 mean_loc(4)-mean_loc(1),...
                 mean_loc(5)-mean_loc(2),...
                -mean_loc(6)+mean_loc(3),...
                'AutoScale','off','LineWidth',4);
h2.MaxHeadSize = 0.2/norm([h2.UData h2.VData h2.WData]);
h2.Color = 'k';
% h3 = quiver3(ax1,plane_cent(2),...
%                  plane_cent(3),...
%                 -plane_cent(4),...
%                  10*km2deg(normal1(1)),...
%                  10*km2deg(normal1(2)),...
%                  10*normal1(3),...
%                 'AutoScale','off','LineWidth',4);
% h3.MaxHeadSize = 0.2/norm([h3.UData h3.VData h3.WData]);
% h3.Color = 'r';
% h4 = quiver3(ax1,plane_cent(2),...
%                  plane_cent(3),...
%                 -plane_cent(4),...
%                  10*km2deg(normal2(1)),...
%                  10*km2deg(normal2(2)),...
%                  10*normal2(3),...
%                 'AutoScale','off','LineWidth',4);
% h4.MaxHeadSize = 0.2/norm([h4.UData h4.VData h4.WData]);
% h4.Color = 'b';
% h5 = quiver3(ax1,plane_cent(2),...
%                  plane_cent(3),...
%                 -plane_cent(4),...
%                  mean_loc(4)-mean_loc(1),...
%                  mean_loc(5)-mean_loc(2),...
%                 -mean_loc(6)+mean_loc(3),...
%                 'AutoScale','off','LineWidth',4);
% h5.MaxHeadSize = 0.2/norm([h2.UData h2.VData h2.WData]);
% h5.Color = 'k';
text(ax1,mean_loc(1)+((mean_loc(4)-mean_loc(1))/2),...
         mean_loc(2)+((mean_loc(5)-mean_loc(2))/2),...
        -mean_loc(3)+((-mean_loc(6)+mean_loc(3))/2),...
         sprintf('%4.2f km/s; %1.2fVs',mean_vel,mean_vel/vs),'Color',h2.Color)
axis(ax1,data_bounds(3:8))
xlabel(ax1,'Lon (°)')
ylabel(ax1,'Lat (°)')
zlabel(ax1,'Depth (km)')

%Slab2 Figure
f2 = figure(2);
f2.Position(1) = f2.Position(1)+600; %Move Figure 2 600 pixels to the right
s_dep_model = dlmread('izu_slab2_dep_02.24.18.xyz',',');
xv = min(s_dep_model(:,1)):0.1:max(s_dep_model(:,1));
yv = min(s_dep_model(:,2)):0.1:max(s_dep_model(:,2));
[X,Y] = meshgrid(xv,yv);
Z = griddata(s_dep_model(:,1),s_dep_model(:,2),s_dep_model(:,3),X,Y);
surf(X,Y,Z,'facecolor','k')
alpha 0.1
hold on;
scatter3(s_loc(1),s_loc(2),-s_loc(3),100,'filled','r')
quiver3(mean_loc(1),...
        mean_loc(2),...
       -mean_loc(3),...
        mean_loc(4)-mean_loc(1),...
        mean_loc(5)-mean_loc(2),...
       -mean_loc(6)+mean_loc(3),...
        'AutoScale','off','LineWidth',4,'Color','k');
data_bounds2 = [min(data(:,1))                    max(data(:,1))...
                min(min(data(:,2)),s_loc(1))-0.5  max(max(data(:,2)),s_loc(1))+0.5...
                min(min(data(:,3)),s_loc(2))-0.5  max(max(data(:,3)),s_loc(2))+0.5...
                min(-max(data(:,4)),-s_loc(3))-50 max(-min(data(:,4)),-s_loc(3))+50];
axis(data_bounds2(3:8))
view(2)

%Pass variables out of function
field1  = 'vp';             value1  = vp;
field2  = 'vs';             value2  = vs;
field3  = 'mean_loc';       value3  = mean_loc;
field4  = 'mean_vel';       value4  = mean_vel;
field5  = 'az';             value5  = az;
field6  = 'rup_dur';        value6  = rup_dur;
field7  = 'rup_len';        value7  = rup_len;
field8  = 'pref_strike';    value8  = pref_strike;
field9  = 'pref_dip';       value9  = pref_dip;
field10 = 'pref_rake';	    value10 = pref_rake;
field11 = 'other_strike';   value11 = other_strike;
field12 = 'other_dip';      value12 = other_dip;
field13 = 'other_rake';     value13 = other_rake;
field14 = 's_loc';          value14 = s_loc;
field15 = 's_dip';          value15 = s_dip;
field16 = 's_str';          value16 = s_str;
field17 = 's_thk';          value17 = s_thk;
field18 = 'vmisfit';        value18 = vmisfit;
field19 = 'vmisfit2';       value19 = vmisfit2;
field20 = 'vangmisfit';     value20 = vangmisfit;
field21 = 'vangmisfit2';    value21 = vangmisfit2;
field22 = 'data_idx_plane'; value22 = data_idx_plane;
parameters = struct(field1, value1,...
                    field2, value2,...
                    field3, value3,...
                    field4, value4,...
                    field5, value5,...
                    field6, value6,...
                    field7, value7,...
                    field8, value8,...
                    field9, value9,...
                    field10,value10,...
                    field11,value11,...
                    field12,value12,...
                    field13,value13,...
                    field14,value14,...
                    field15,value15,...
                    field16,value16,...
                    field17,value17,...
                    field18,value18,...
                    field19,value19,...
                    field20,value20,...
                    field21,value21,...
                    field22,value22);