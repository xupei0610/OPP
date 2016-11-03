function [ placement, msg ] = OPP_GThN( A, ZI_buses, varargin )
% OPP_GThN Optimal PMU Placement obtained via Graph Theoretic Procedure and
% nonlinear constraint function method
%
%   This function uses topological observability with the consideration
%   of Zero Injection buses. Its strategy is quite similar to DFS, excpet
%   that it uses nonlinear constraint function to determine a bus's
%   observability, and may try to install PMUs at zero injection buses.
%
%   OFF_GThN(A, ZI_buses) returns the optimal PMU placement with respect to 
%   the given adjacent matrix A, whose zero injection buses' ids are stroed
%   in the vecotr ZI_buses, through Graph Theoretic Procedure method.
%
%   Example:
%   If
%       A = [[1 1 0];[1 1 1];[0 1 1]]; ZI = [2]
%   Then
%       [placement, ~] = OOP_GThN(A, ZI) returns placement := [2]
%   Which means the given power system whose adjacent matrix is A will be
%   fully observable if a PMU is installed at Bus 2. This algorithm may
%   obtained a more optimal solution compared with OPP_GTHM who merges each
%   zero injection bus to one of its adjacent buses when calculating the
%   system's observability. However, this algorithm is slower than than the
%   latter.
%
%   See also OPP_DFS, OPP_GTHM, OPP_SA, OPP_SAB, OPP_RSN
%
%   Author: Pei Xu
%   Copyright: 2015~ Pei Xu
%   License: MIT
%


%% Input Check
narginchk(1, inf);
if nargin() < 2
    ZI_buses = [];
end

%% The name of the method currently is employed.
msg.method = 'Graphic Theoretic Procedure using Nonlinear Constraint Function Method';

%% Initialize Timer
GThB_timer = tic;

%% The amount of buses in the given system
num_buses = length(A);

%% Construct X vector
% A vector composed of 1 and 0, where 1 represents a PMU is installed at
% that bus and 0 otherwise.
X = zeros(num_buses, 1);

%% Initialize the vector which records each bus's influence power.
% The bigger the value is, the more buses will be observable if a PMU is
% installed at the bus the value represents.
in_degree = full(sum(A, 2));

%% Change the adjacent matrix into the vector format for searching speed
adj = cell(num_buses, 1);
for i=1:1:num_buses
    A(i,i)=0;
    adj{i} = find(A(i,:)~=0);
    A(i,i)=1;
end

%% Construct Topological Observability Formula
% Construct the linear function without consideration of ZI buses
equ = cell(1,num_buses);
for i=1:1:num_buses
    equ{i} = sprintf('f(%d)=X(%d)%s;', i, i, sprintf('+X(%d)', adj{i}));
end
% Revise the linear constraint function into a nonlinear function.
for i=1:1:length(ZI_buses)
    equ{ZI_buses(i)} = sprintf('%s+~any([%s]==0);', equ{ZI_buses(i)}(1:length(equ{ZI_buses(i)})-1), sprintf('f(%d) ', adj{ZI_buses(i)}));
    
    num_of_incidents = length(adj{ZI_buses(i)});
    for j=1:1:num_of_incidents
        incidents = adj{ZI_buses(i)};
        incidents(j) = [];
        incidents(num_of_incidents) = ZI_buses(i);
        
        equ{adj{ZI_buses(i)}(j)} = sprintf('%s+~any([%s]==0);', equ{adj{ZI_buses(i)}(j)}(1: length(equ{adj{ZI_buses(i)}(j)})-1), sprintf('f(%d) ', incidents)); 
        
    end
end
% Construct the function
constraint_function = strjoin(equ, '\n');
function f = F(X)
    f=zeros(num_buses,1);
    c_f = f;
    % due to the coupling, constraint function should be calculated several
    % times until its value could not change anymore.
    while true
        eval(constraint_function);
        if c_f == f
            break;
        else
            c_f = f;
        end
    end
end

%% Search Begins
% Its procedure is same to DFS except that it employs a different formula
% to calculate observability

% This algorithm will must success when not consider redundancy in observability
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
    ob = F(X);
    
    % If the system is fully observable, search ends.
    if ~any(ob < 1 )
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
            if A(i, j) == 0 || ob(j) > 0
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
msg.time = toc(GThB_timer);

end