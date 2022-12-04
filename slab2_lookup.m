function [s_loc,s_dip,s_str,s_thk] = slab2_lookup(in_lat,in_lon,in_dep)

%Find slab2 values at closest point to (in_lat,in_lon,in_dep)
%Previous version found values at a given (in_lat,in_lon) regardless of
%in_dep (i.e., projected to the slab interface)

s_dep_model = dlmread('izu_slab2_dep_02.24.18.xyz',',');
s_dip_model = dlmread('izu_slab2_dip_02.24.18.xyz',',');
s_str_model = dlmread('izu_slab2_str_02.24.18.xyz',',');
s_thk_model = dlmread('izu_slab2_thk_02.24.18.xyz',',');

R=6371-in_dep;
[~,idx] = min(sum([deg2km(abs(s_dep_model(:,1)-in_lon),R) deg2km(abs(s_dep_model(:,2)-in_lat),R) abs(-s_dep_model(:,3)-in_dep)],2));
s_loc = [s_dep_model(idx,1) s_dep_model(idx,2) -s_dep_model(idx,3)];
s_dip = s_dip_model(idx,3);
s_str = s_str_model(idx,3);
s_thk = s_thk_model(idx,3);

% [~,s_dep_idx] = min(abs(s_dep_model(:,1) - in_lon) + abs(s_dep_model(:,2) - in_lat));
% s_dep = -s_dep_model(s_dep_idx,3);
% 
% s_dip_model = dlmread('izu_slab2_dip_02.24.18.xyz',',');
% [~,s_dip_idx] = min(abs(s_dip_model(:,1) - in_lon) + abs(s_dip_model(:,2) - in_lat));
% s_dip = s_dip_model(s_dip_idx,3);
% 
% s_str_model = dlmread('izu_slab2_str_02.24.18.xyz',',');
% [~,s_str_idx] = min(abs(s_str_model(:,1) - in_lon) + abs(s_str_model(:,2) - in_lat));
% s_str = s_str_model(s_str_idx,3);
% 
% s_thk_model = dlmread('izu_slab2_thk_02.24.18.xyz',',');
% [~,s_thk_idx] = min(abs(s_thk_model(:,1) - in_lon) + abs(s_thk_model(:,2) - in_lat));
% s_thk = s_thk_model(s_thk_idx,3);