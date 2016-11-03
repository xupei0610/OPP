function [ placement, msg ] = OPP_SAB( A, ZI_buses, varargin )
% OPP_SAB Optimal PMU Placement obtained via the modified version of Simulated
% Annealing method who increasing the randomness of the original version
% and who reducing the time needed.
%
%   This function uses topological observability with the consideration
%   of Zero Injection buses. It uses the basic idea proposed in Power System 
%   Observability with Minimal Phasor Measurement Placement, 1993, but
%   different stop cerition and way to generate random placement.
%
%   OFF_SAB(A, ZI_buses) returns the optimal PMU placement with respect to 
%   the given adjacent matrix A, whose zero injection buses' ids are stroed
%   in the vecotr ZI_buses, through a modified Simulated Annealing method.
%
%   Attention: Due to the nature of Simulated Anealing method as a
%   stochastic algorithm, the solution or placement found by this method is
%   not fixed, which means it may be the optimal but also may not.
%
%   See also OPP_DFS, OPP_GTHM, OPP_GTHN, OPP_SAB, OPP_RSN
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
msg.method = 'Simulated Annealing Method (Modified Version proposed in the report)';

%% Initialize the timer, and the search will begin right now.
SAB_timer = tic;

%% The amount of buses in the given system
num_buses = length(A);

%% Initialize X vector
% A vector composed of 1 and 0, where 1 represents a PMU is installed at
% that bus and 0 otherwise.
% In the orignal version, a solution obtained via Graph Theoretical
% Procedure is used as the initial solution.
[p, ~] = OPP_GThN(A, ZI_buses);
X = zeros(1, num_buses);
X(1, p) = 1;

%% Define Energy Function
% E = the number of buses - the number of observable buses
% Namely, E = the amount of unobservable buses
E = @(observability) length(find(observability<1));

%% Initialize Temperature Parameter
T0 = 15;

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

%% Search Starts

% Consider the initial solution as the best solution found so far
X_best = X;

T = T0;
T_diminished = 0;
while T_diminished < num_buses
    X = X_best;
    E_X = 0;
    for i = 1:1:num_buses
        % Achieve a new placement
        X_new = perturb(X);
        
        % the observabiltiy of the newly generated placement
        ob = F(X_new);
        
        % the the number of unobservable buses when the newly generated
        % placement is applied
        E_nX = E(ob);
        
        % If fully observable when the newly generated placment is applied
        if E_nX == 0 && nnz(X_new) < nnz(X_best)
            % better solution
            X_best = X_new;
            T_diminished = -1;
            break;
        end
        
        delta_E = E_nX - E_X;

        % reject the worse, new solution in probability
        if delta_E > 0 && rand > exp(-1 * delta_E / T)
            continue;
        end
        
        X = X_new;
        E_X = E_nX;

    end
    % The times of continuously failing 
    T_diminished = T_diminished + 1;
    % Annealing
    T = 0.879 * T;
end

%% Arrange the result
placement = find(X_best == 1);

%% Timer should be over
msg.time = toc(SAB_timer);

end

%% Randomly generate a new configuration from the given configuration
function new_configuration = perturb(X)
    os = find(X==1);
    zs = find(X==0);
    if or(or(isempty(zs), isempty(os)), rand > 0.5)
        % Flip
        new_configuration = X;
        tar_index = unidrnd(length(X));
        if X(tar_index) == 0
            new_configuration(tar_index) = 1;
        else
            new_configuration(tar_index) = 0;
        end
    else
        % Exchange
        candidate1 = zs(unidrnd(length(zs)));
        candidate2 = os(unidrnd(length(os)));
        new_configuration = X;
        new_configuration(candidate1) = 1;
        new_configuration(candidate2) = 0;
    end
end
 