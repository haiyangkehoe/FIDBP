%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%9 June 2022
%Modified 28 June 2022

function fidbp_plot(hyp,cont_t,cont_locs,data,rupv_acut)

%Determine information from inputs
tsteps = unique(cont_t) + hyp(1);
data_bounds = [min(data(:,1))  max(data(:,1))...
               min(data(:,2))  max(data(:,2))...
               min(data(:,3))  max(data(:,3))...
              -max(data(:,4)) -min(data(:,4))];
tot_pts = size(data,1)/data_bounds(2);

%Plot distance vs. time
f1 = figure(1);
ax3 = axes;
hold on;
if any(data(:,5))
    cont_dists = sqrt(deg2km(cont_locs(:,1)).^2 + deg2km(cont_locs(:,2)).^2 + cont_locs(:,3).^2);
    cont_times = cont_t + hyp(1);
    cont_amps = ones(size(cont_dists,1),1);
    scatter(ax3,hyp(5)*((1-hyp(1))+(cont_times-1)),cont_dists,cont_amps*100,'k','Marker','+')
    idx = find(data(:,5)>=rupv_acut);
    data_idx_rel = [data(idx,1) data(idx,2)-hyp(2) data(idx,3)-hyp(3) -data(idx,4)+hyp(4) data(idx,5)];
    dists = sqrt(deg2km(data_idx_rel(:,2)).^2 + deg2km(data_idx_rel(:,3)).^2 + data_idx_rel(:,4).^2);
    times = data_idx_rel(:,1);
    amps  = data_idx_rel(:,5);
    scatter(ax3,hyp(5)*((1-hyp(1))+(times-1)),dists,amps*100,amps,'filled','MarkerEdgeColor','k');
    colormap(ax3,haxby_hk)
    caxis(ax3,[0 1])
    colorbar(ax3,'Position',[.92 .11 .03 .815])
    axis(ax3,[hyp(5)*((1-hyp(1))+(cont_times(1)-1)) hyp(5)*((1-hyp(1))+(cont_times(end)-1)) min(cont_dists) max(cont_dists)])
    xlabel(ax3,'Time (s)')
    ylabel(ax3,'Distance (km)')
    %Plot vp, vs lines
    plotx = 0:0.1:max(hyp(5)*((1-hyp(1))+(times-1)))+1;
    [vp,vs,depth] = iasp91_lookup(hyp(4));
    ploty_vp = vp*plotx;
    ploty_vs = vs*plotx;
    plot(ax3,plotx,ploty_vp,'b')
    plot(ax3,plotx,ploty_vs,'r')
    legend('Inversion grid','BP Sources',...
           sprintf('vp=%f km/s at %i km depth',vp,depth),...
           sprintf('vs=%f km/s at %i km depth',vs,depth),'Location','northwest')
else
    dists = sqrt(deg2km(cont_locs(:,1)).^2 + deg2km(cont_locs(:,2)).^2 + cont_locs(:,3).^2);
    times = cont_t + hyp(1);
    amps = ones(size(dists,1),1);
    colormap(jet(length(tsteps)));
    scatter(ax3,hyp(5)*((1-hyp(1))+(times-1)),dists,amps*100,hyp(5)*((1-hyp(1))+(times-1)),'Marker','+');
    axis(ax3,[hyp(5)*((1-hyp(1))+(times(1)-1)) hyp(5)*((1-hyp(1))+(times(end)-1)) min(dists) max(dists)])
    xlabel(ax3,'Time (s)')
    ylabel(ax3,'Distance (km)')
    %Plot vp, vs lines
    plotx = 0:0.1:max(hyp(5)*((1-hyp(1))+(times-1)))+1;
    [vp,vs,depth] = iasp91_lookup(hyp(4));
    ploty_vp = vp*plotx;
    ploty_vs = vs*plotx;
    plot(ax3,plotx,ploty_vp,'b')
    plot(ax3,plotx,ploty_vs,'r')
    legend('BP Sources',...
           sprintf('vp=%f km/s at %i km depth',vp,depth),...
           sprintf('vs=%f km/s at %i km depth',vs,depth),'Location','northwest')
end
ButtonHandle = uicontrol('Style', 'togglebutton','String', 'Continue');
while ButtonHandle.Value==0
    fprintf('Click two points to measure slope OR Press "Continue" to play BP movie\n')
    [x1,y1] = ginput(1);
    pause(0.1)
    if ButtonHandle.Value==1
        break
    end
    [x2,y2] = ginput(1);
    pause(0.1)
    if ButtonHandle.Value==1
        break
    end
    slope = (y1-y2)/(x1-x2);
    if exist('h1','var')
        delete(h1)
    end
    h1 = plot(ax3,[x1 x2], [y1 y2],'k');
    legend(ax3,'Inversion grid','BP Sources',...
           sprintf('V_p=%f km/s at %i km depth',vp,depth),...
           sprintf('V_s=%f km/s at %i km depth',vs,depth),...
           sprintf('V_{sel}=%f km/s',slope),...
           'Location','northwest')
    fprintf(sprintf('Slope = %f km/s\n',slope))
end

%Plot back-projection movie
fprintf('Playing BP movie... Close figure to exit.\n')
f2 = figure(2);
f2.Position(1) = f2.Position(1)+600; %Move Figure 2 600 pixels to the right
ax1 = axes;
ax2 = axes;
linkaxes([ax1,ax1])
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax1.ZTick = [];
myColormap = colormap(ax1,jet(length(tsteps)));
colormap(ax2,haxby_hk);
if tsteps(1) < tsteps(end)
    colorbar(ax1,'southoutside','Position',[.13 .04 .78 .03]);
    caxis(ax1,[hyp(5)*((1-hyp(1))+(tsteps(1)-1)) hyp(5)*((1-hyp(1))+(tsteps(end)-1))])
end
colorbar(ax2,'Position',[.92 .11 .03 .815]);
caxis(ax2,[0 1])
daspect([1 1 100])
hold on;
while true
    scatter3(ax2,hyp(2),hyp(3),-hyp(4),250,'filled','r','p');
    for t = 1:length(tsteps)
        if ishghandle(f1)==0 || ishghandle(f2)==0
            close all;
            break
        else
            t_idxs = find(cont_t == tsteps(t)-hyp(1));
            h1 = scatter3(ax2,cont_locs(t_idxs,1)+hyp(2),cont_locs(t_idxs,2)+hyp(3),-cont_locs(t_idxs,3)-hyp(4),...
                          100,myColormap(t,:),'Marker','+');
            idx = find(data(tot_pts*(tsteps(t)-1)+1:tot_pts*tsteps(t),5)>=rupv_acut);
            scatter3(ax2,data(tot_pts*(tsteps(t)-1)+idx,2),...
                         data(tot_pts*(tsteps(t)-1)+idx,3),...
                        -data(tot_pts*(tsteps(t)-1)+idx,4),...
                         data(tot_pts*(tsteps(t)-1)+idx,5)*100,...
                         data(tot_pts*(tsteps(t)-1)+idx,5),...
                         'filled',...
                         'MarkerEdgeColor','k');
            axis(ax2,data_bounds(3:8))
            title(ax2,sprintf('t = %5.2f s',hyp(5)*((1-hyp(1))+(tsteps(t)-1))))
            xlabel(ax2,'Lon (°)')
            ylabel(ax2,'Lat (°)')
            zlabel(ax2,'Depth (km)')
            drawnow
            pause(hyp(5)/10)
            delete(h1)
        end
    end
    pause(5)
    if ishghandle(ax2)
        cla(ax2)
    elseif ishghandle(f1)==0 || ishghandle(f2)==0
        close all;
        break
    end
end

