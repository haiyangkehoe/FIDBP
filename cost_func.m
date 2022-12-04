%Haiyang Kehoe
%University of Arizona
%Department of Geosciences
%13 September 2022
%Modified 02 November 2022

function [best_cost,best_idx,cost] = cost_func(data,run_data,textfilename)

lond = round(mean(diff(unique(data(:,2)))),4);
latd = round(mean(diff(unique(data(:,3)))),4);
depd = round(mean(diff(unique(data(:,4)))),4);
min_rup_len = 2*ones(size(run_data,1),1).*max([deg2km(lond) deg2km(latd) depd]) + 0.0001;
min_rup_dur = 2*run_data(:,2) + 0.0001;
cost = (run_data(:,18)) ./ ...
       ((heaviside(run_data(:,5) - 4999.999)) .* ...
        (heaviside(run_data(:,8) - min_rup_dur)) .* ...
        (heaviside(run_data(:,9) - min_rup_len)));

%Find best run for each BP time-step   
tdirs = unique(run_data(:,2));
best_costs = zeros(size(tdirs,1),1);
best_idxs  = zeros(size(tdirs,1),1);
for i = 1:size(tdirs,1)
    [best_costs(i),best_idxs(i)] = min(cost(:,1)./(run_data(:,2)==tdirs(i)));
    fprintf('Best run for %3.2fs BP timestep at index %i with cost = %3.2e\n',...
            tdirs(i),best_idxs(i),best_costs(i))
end

%Write best idx for each time step
dlmwrite('best_idxs.txt',[tdirs best_idxs],'delimiter','\t')

%Find best run for all BP time-steps
[best_cost,best_idx] = min(cost);

%Save textfile
file = fopen(textfilename, 'w');
for i = 1:size(run_data, 1)
    fprintf(file, '%i %3.2f %3.2f %3.2f %i %3.2f %3.2f %3.2f %3.2f %i %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2e\n',...
                  [run_data(i,:) cost(i)]);
end
fclose(file);