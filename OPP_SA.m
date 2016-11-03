function [ placement, msg ] = OPP_SA( A, ZI_buses, varargin )
% OPP_SA Optimal PMU Placement obtained via the original version of Simulated
% Annealing method.
%
%   This function uses topological observability with the consideration
%   of Zero Injection buses. It uses the strategy proposed in Power System 
%   Observability with Minimal Phasor Measurement Placement, 1993. 
%
%   OFF_SA(A, ZI_buses) returns the optimal PMU placement with respect to 
%   the given adjacent matrix A, whose zero injection buses' ids are stroed
%   in the vecotr ZI_buses, through Simulated Annealing method.
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
msg.method = 'Simulated Annealing Method (Original Version proposed in 1993)';

%% Initialize the timer, and the search will begin right now.
SA_timer = tic;

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

% Determine the inital upper and lower bound
upper_bound = nnz(X);
lower_bound = 0;

% Warning message may appear due to a too huge value of M
warning('off','all');

% Search really starts
while upper_bound-lower_bound > 1
    
    % Calculate v_test
    if lower_bound <= 0
        mid_point_coefficient = 0.85;
    else
        mid_point_coefficient = 0.5;
    end
    v_test = floor(mid_point_coefficient*(upper_bound-lower_bound))+lower_bound;

    % Load initial temperature
    T = T0;
    
    % Calculate M
    M = ceil(0.002 * nchoosek(num_buses, v_test));
    
    % Flag for if a new acceptable solution is found
   	observable = false;
    
    % Randomly generate a set of placement in which the amount of PMU
    % installed is equal to v_test
    c_X = zeros(1, num_buses);
    c_X(1, randperm(num_buses, v_test)) = 1;
    E_cX = E(c_X);

    if E_cX == 0 % This is almost 'unhappenable'
        X_best = c_X;
        observable = true;
    else
        % Execute SA
        for i=1:1:40
            X = c_X;
            E_X = E_cX;
            for j=1:1:M
                %fprintf('i: %d; j: %d; v_test: %d\n', i, j, v_test)
                % Achieve a new placement via randomly
                % moving a PMU to a bus without PMU     
                os = find(X==1);
                zs = find(X==0);
                candidate1 = zs(unidrnd(length(zs)));
                candidate2 = os(unidrnd(length(os)));
                X_new = X;
                X_new(candidate1) = 1;
                X_new(candidate2) = 0;

                % the observabiltiy of the newly generated placement
                ob = F(X_new);
                
                % the the number of unobservable buses when the newly generated
                % placement is applied
                E_nX = E(ob);

                % If fully observable, return
                if E_nX == 0
                    X_best = X_new;
                    observable = true;
                    break;
                end

                delta_E = E_nX - E_X;

                % reject the worse, new solution in probability
                if delta_E > 0 && rand > exp(-1 * delta_E / T)
                    break;
                end
                
                % accept the new solution
                X = X_new;
                E_X = E_nX;
            end
            if observable == true;
                break;
            end
            % Annealing
            T = 0.879 * T;
        end
    end
    % Adjust bound
    if observable == true
        upper_bound = v_test;
    else
        lower_bound = v_test;
    end

end


warning('on', 'all');

%% Arrange the result
placement = find(X_best == 1);


%% Timer should be over
msg.time = toc(SA_timer);

end