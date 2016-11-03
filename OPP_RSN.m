function [ placements, msg ] = OPP_RSN( A, ZI_buses, varargin )
% OPP_RSN Optimal PMU Placement obtained via Recursive Security N Algorithm
%
%   This function uses topological observability with the consideration
%   of Zero Injection buses. It uses the idea of the minimum spanning tree.
%   This algorithm is proposed in A Security Oriented Approach to PMU
%   Positioning for Advanced Monitoring of a Transmission Grid, 2002.
%
%   OFF_RSN(A) returns the optimal PMU placement with respect to the given
%   adjacent matrix A via RSN method.
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
%   See also OPP_DFS, OPP_GTHM, OPP_GTHN, OPP_SA, OPP_SAB
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

%% Initialize the timer, and the search will begin right now.
RSN_timer = tic;

%% The name of the method currently is employed.
msg.method = 'Recursive Security N Algorithm';

%% The amount of buses in the given system
num_buses = length(A);

%% Reverse the adjacent matrix
% via the symmetric reverse Cuthill-McKee permutation
CM_p = symrcm(A);
A = A(CM_p,CM_p);
for i=1:1:length(ZI_buses)
    ZI_buses(i) = find(CM_p==ZI_buses(i));
end

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
constraint_function = strjoin(equ, '\n');

function f = F(X)
    f=zeros(num_buses,1);
    af=zeros(num_buses,1);
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

%% Obtain initial solutions
X_best = ones(1, num_buses);
xb_len = 1;
xb_used = num_buses;
cand_X = nan(num_buses, num_buses);
cX_len = 0;

tmp = eye(num_buses);
while ~isempty(tmp)
    c_X = transpose(tmp(1,:));
    tmp(1,:)=[];
    obd = find(A*c_X>0);
    if length(obd) == num_buses
        cX_len = cX_len + 1;
        cand_X(cX_len, :) = c_X;
        continue;
    end
    in_degree = sum(A, 2);
    for i=1:1:length(obd)
        in_degree(obd(i)) = 0;
    end
    
    [~, max_index] = max(in_degree);
    %fprintf('handled: %d; remains: %d\n', nnz(tmp(1,:)), size(tmp, 1));
    tmp_X = transpose(c_X);
    tmp_X(max_index) = 1;
    if ~any(ismember(tmp, tmp_X, 'rows'))
    	tmp(size(tmp, 1)+1,:) = tmp_X;
    end
end

cand_X(cX_len+1:end,:) = [];
cand_X = unique(cand_X, 'rows');


%% Generate Alternative Patterns
Xs = zeros(0, num_buses);
Xs_len = 0;
tic
while ~isempty(cand_X)
    c_X = cand_X(1,:);
    cand_X(1,:) = [];
    os = find(c_X==1);
    add = false;
    all_continue = true;
    for i=1:1:length(os)
        for j=1:1:length(adj{os(i)})
            if c_X(adj{os(i)}(j)) == 1
                continue;
            end
            all_continue = false;
            tmp = c_X;
            tmp(os(i)) = 0;
            tmp(adj{os(i)}(j)) = 1;
            if ~any(F(tmp) < 1)
                %Xs_len = Xs_len + 1;
                %Xs(Xs_len,:)= tmp;
                add = true;
                break;
            end
        end
        if add==true
            break;
        end
    end
    if all_continue == false && add == false;
        continue;
    else
        Xs_len = Xs_len + 1;
        Xs(Xs_len,:) = c_X;
    end
end
%Xs = unique(Xs, 'rows');

%% Remove Pure Transit Node
%cand_X = Xs;
cand_X = Xs(find(sum(Xs,2) == min(sum(Xs, 2))),:);
i=0;

while ~isempty(cand_X)

    % Each Iteration, the placements use maximum PMUs will be dealt firstly
    % So, we can avoid keeping a large explored or closed list for graph search
%     maxs = find(sum(cand_X, 2)==max(sum(cand_X, 2)));
%     Xs = cand_X(maxs, :);
%     cand_X(maxs, :) = [];
    Xs = cand_X;
    cand_X = zeros(0,num_buses);

%     i = i+1;
%     fprintf('iter: %d\t', i);
%     fprintf('length_cx: %d\t', nnz(Xs(1,:)));
%     fprintf('number_cs: %d\t', size(Xs, 1));
%     fprintf('remains: %d\n', size(cand_X, 1))

    for k=1:1:size(Xs, 1)
        os = find(Xs(k,:)==1);
        for l=1:1:length(os)
            ttmp_X = Xs(k, :);
            ttmp_X(os(l)) = 0;
            
            % check if the placement is already in the candidate list
            repeat = find(ismember(cand_X, ttmp_X, 'rows'), 1);
            if ~isempty(repeat)
                cand_X(repeat,:) = [];
            end            
            
            observability = F(ttmp_X);
            
            % check if it is real a pure tansit node
            if any(observability < 1)
                % discard if it in fact is not
                continue;
            end
            % check if better than the current best placement
            c_used = nnz(ttmp_X);
            if xb_used > c_used
                X_best = ttmp_X;
                xb_len = 1;
                xb_used = c_used;
            elseif xb_used == nnz(ttmp_X)
                xb_len = xb_len + 1;
                X_best(xb_len,:) = ttmp_X;
            end
            % store the placement for continuing checking if it still has
            % any pure transit node
            cand_X(size(cand_X,1)+1,:) = ttmp_X;
        end
    end
    
end

X_best = unique(X_best, 'rows');

%% Arrange solutions
placements = nan(size(X_best, 1), nnz(X_best(1,:)));
for i=1:1:size(X_best, 1)
    os = find(X_best(i,:)==1);
    for j=1:1:length(os)
        placements(i, j) = CM_p(os(j));
    end
end

%% Timer should be over
msg.time = toc(RSN_timer);

end

