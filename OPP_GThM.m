function [ placement, msg] = OPP_GThM( A, ZI_buses, varargin  )
% OPP_GThM Optimal PMU Placement obtained via Graph Theoretic Procedure
%
%   This function uses topological observability with the consideration
%   of Zero Injection buses.
%   This function will try to merge each zero injection bus into one of its
%   adjecent buses, and then use DFS method to deal with the new matrix.
%   A more optimal result may be obtained via OPP_GThN who uses nonlinear
%   inequality to deal with Zero Injection buses.
%
%   OFF_GThM(A, ZI_buses) returns the optimal PMU placement with respect to 
%   the given adjacent matrix A, whose zero injection buses' ids are stroed
%   in the vecotr ZI_buses,through Graph Theoretic Procedure method.
%
%   Example:
%   If
%       A = [[1 1 0];[1 1 1];[0 1 1]]; ZI = [2];
%   Then
%       [placement, ~] = OOP_GThM(A, ZI) returns placement := [2]
%   Which means the given power system whose adjacent matrix is A will be
%   fully observable if a PMU is installed at Bus 2.
%
%
%   See also OPP_DFS, OPP_GTHN, OPP_SA, OPP_SAB, OPP_RSN
%
%   Author: Pei Xu
%   Copyright: 2015~ Pei Xu
%   License: MIT
%

%% Input Check
narginchk(1, inf);

if nargin < 2
    ZI_buses = [];
end
%% Revise the data returned by DFS
msg.method = 'Graphic Theoretic Procedure using Merger Method';

%% Initialize timer
GThM_timer = tic;

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

%% Reconstruct A matrix
% via meraging each ZI bus to one of its adjacent buses
mod_A = A;
%merger_target = zeros(length(ZI_buses), 1);
handled = false;
for i=1:1:length(ZI_buses)
    for j=1:1:num_buses
        if any(ZI_buses == j)
            continue;
        else
            mod_A(j,:) = or(A(i,:), A(j,:));
            mod_A(i,:) = mod_A(j,:);
            %merger_target(i) = j;
            handled = true;
            break;
        end
    end
    if handled == false
        %% ZI bus encirciled by ZI buses
        mod_A(i, :) = zeros(1, num_buses);
    end
    handled = false;
end

%% Search starts

% GTh will must success when we do not consider redundancy in observability
msg.success = true;

for k=1:1:num_buses
    % We need to find a bus with the most branches in the unobservable
    % area. Namely, the bus whose in-degree is the maxmium.
    
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

    X(index) = 1;
    
    % Recalcuate the whole system's observability
    F = mod_A*X;
    
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

%% Search Ends
% Unable to satisfy the requirement of redundancy
% if msg.success == false
%     placement = [];
% else
    % Arrange the placement vector
placement = find(transpose(X)==1);
% end

%% Timer should be over
msg.time = toc(GThM_timer);

