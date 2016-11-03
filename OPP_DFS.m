function [ placement, msg ] = OPP_DFS( A, varargin )
% OPP_DFS Optimal PMU Placement obtained via DFS
%
%   This function uses topological observability without the consideration
%   of Zero Injection buses. 
%
%   OFF_DFS(A) returns the optimal PMU placement with respect to the given
%   adjacent matrix A via DFS method.
%
%   Example:
%   If
%       A = [[1 0 1];[0 1 1];[1 1 1]];
%   Then
%       [placement, ~] = OOP_DFS(A) returns placement := [3]
%   Which means the given power system whose adjacent matrix is A will be
%   fully observable if a PMU is installed at Bus 3 (the bus represented
%   via the 3rd row and column in the A matrix).
%
%
%   See also OPP_GTHM, OPP_GTHN, OPP_SA, OPP_SAB, OPP_RSN
%
%   Author: Pei Xu
%   Copyright: 2015~ Pei Xu
%   License: MIT
%

%% Input Check
narginchk(1, inf);

%% The name of the method currently is employed.
msg.method = 'Depth First Search';

%% Initialize the timer, and the search will begin right now.
DFS_timer = tic;

%% The amount of buses in the given system
num_buses = length(A);

%% Construct X vector
% A vector composed of 1 and 0, where 1 represents a PMU is installed at
% that bus and 0 otherwise.
X = zeros(num_buses, 1);

%% Initialize the vector which records each bus's influence power.
% The bigger the value in the array is, the more buses will be observable 
% if a PMU is installed at the bus which the index of the value represents
in_degree = sum(A, 2);

%% Search starts

% DFS will must success when we do not consider redundancy in observability
msg.success = true;

for k=1:1:num_buses

    % We need to find a bus with the most branches in the unobservable
    % area. Namely, the bus whose in-degree is the maxmium.
    [~, index] = max(in_degree);
%     [val, index] = max(in_degree);

    % This could not happen if we do not consider redundancy in
    % observability
%     if val < 1
%         msg.success = false;
%         break;
%     end

    % Put a PMU at the found bus
    X(index) = 1;

    % Recalcuate the whole system's observability
    F = A*X;
    
    % If the system is fully observable, search ends.
    if ~any(F < 1)
        %msg.success = true;
        break;
    end
    
    % The system is still not fully observable, we should revise in-degree 
    % for each bus influenced via the newly installed PMU.
    in_degree(index) = 0;
    for i=1:1:num_buses
        if A(index, i) == 0 || in_degree(i) == 0
            continue;
        end
        in_degree(i) = 0;
        for j=1:1:num_buses
            if A(i, j) == 0 || F(j) > 0
                continue;
            end
            in_degree(i) = in_degree(i) + 1;
        end
    end
end

% The below part is to consider the situation where more than one candidate
% bus may appear.
% Nevertheless, it would not effective imporve the optimality of DFS,
% although significantly increase its time complexity.
%
% cand_X = eyes(1, num_buses);
% in_degree = transpose(sum(A, 2));
% cX_len = 1;
% X_best = eye(1, num_buses);
% bX_len = 1;
% bX_pmu = num_buses;
% while ~isempty(cand_X)
%     c_X = cand_X(1,:);
%     cand_X(1,:) = [];
%     cX_len = cX_len - 1;
%     c_in_deg = in_degree(1,:);
%     in_degree(1,:) = [];
%     
%     if nnz(c_X) + 1 > bX_pmu
%         continue;
%     end
% 
%     val = max(c_in_deg);
%     max_ind = find(c_in_deg==val);
%     for i=1:1:length(max_ind)
%         tmp_X = c_X;
%         tmp_X(max_ind(i)) = 1;
%         F = A*transpose(tmp_X);
%         if ~any(F < 1)
%            if any(ismember(X_best,tmp_X, 'rows'))
%                continue;
%            end
%            if nnz(tmp_X) < bX_pmu
%                X_best = tmp_X;
%                bX_len = 1;
%                bX_pmu = nnz(tmp_X);
%            else
%                bX_len = bX_len + 1;
%                X_best(bX_len, :) = tmp_X;
%            end
%            continue; 
%         end
%         tmp_in_deg = c_in_deg;
%         tmp_in_deg(max_ind(i)) = 0;
%         for j=1:1:num_buses
%             if A(max_ind(i), j) == 0 || tmp_in_deg(j) == 0
%                 continue;
%             end
%             tmp_in_deg(j) = 0;
%             for k=1:1:num_buses
%                 if A(j, k) == 0 || F(k) > 0
%                     continue;
%                 end
%                 tmp_in_deg(j) = tmp_in_deg(j) + 1;
%             end
%         end
%         
%         cX_len = cX_len+1;
%         cand_X(cX_len, :) = tmp_X;
%         in_degree(cX_len, :) = tmp_in_deg;
%     end
% end

%% Search Ends
% Unable to satisfy the requirement of redundancy
% if msg.success == false
%     placement = [];
% else
    % Arrange the placement vector
placement = find(transpose(X)==1);
% placement = zeros(bX_len, nnz(X_best(1,:)));
% for i=1:1:bX_len
%     placement(i,:) = find(X_best(i,:)==1);
% end
% end

%% Timer should be over
msg.time = toc(DFS_timer);

