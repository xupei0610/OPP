%
% Run Me for single test
%

%%
disp('Which Test Case do you want to use?');
disp('1. IEEE 14-bus');
disp('2. IEEE 30-bus');
disp('3. IEEE 57-bus');
disp('4. IEEE 118-bus');
disp('5. New England 39-bus');
while true
    c = input('Pleas input 1, 2, 3, 4, or 5: ');
    if c == 1
        data = case14();
        test_case = 'IEEE 14-bus System';
        break;
    elseif c == 2
        data = case30();
        test_case = 'IEEE 30-bus System';
        break;
    elseif c == 3
        data = case57();
        test_case = 'IEEE 57-bus System';
        break;
    elseif c == 4
        data = case118();
        test_case = 'IEEE 118-bus System';
        break;
    elseif c == 5
        data = case39();
        test_case = 'New England 39-bus System';
        break;
    end
end
fprintf('\n')
disp('Which Algorithm do you want to use?');
disp('1. Depth First Search');
disp('2. Graphic Theoretic Procedure using Merger Method');
disp('3. Graphic Theoretic Procedure using Nonlinear Constraint Function Method');
disp('4. Original Simulated Annealing Method');
disp('5. Modified Simulated Annealing Method');
disp('6. Recursive Security N Algorithm');
while true
    c = input('Pleas input 1, 2, 3, 4, 5 or 6: ');
    if c == 1
        alg = 'DFS';
        break;
    elseif c == 2
        alg = 'GThM';
        break;
    elseif c == 3
        alg = 'GThN';
        break;
    elseif c == 4
        alg = 'SA';
        break;
    elseif c == 5
        alg = 'SAB';
        break;
    elseif c == 6
        alg = 'RSN';
        break;
    end
end
clear c;

%% Load case

% The number of buses in the system
[num_buses, ~] = size(data.bus);

% Construct A matrix
% This is not the admittance matrix but the adjacent matrix where the
% diagonal elements are all 1.
[branch_len, ~] = size(data.branch);
r1 = zeros(branch_len, 1);
r2 = r1; c1 = r1; c2 = r2; v1 = r1; v2 = r2;
for i=1:1:branch_len
    r1(i) = data.branch(i, 1);
    c1(i) = data.branch(i, 2);
    v1(i) = 1;
    r2(i) = data.branch(i, 2);
    c2(i) = data.branch(i, 1);
    v2(i) = 1;
end
diag_vector = transpose(1:1:num_buses);
r = [r1; r2; diag_vector];
c = [c1; c2; diag_vector];
v = [v1; v2; diag(eye(num_buses, num_buses))];
A = sparse(r, c, v, num_buses, num_buses);
clear r r1 r2 c c1 c2 v v1 v2 diag_vector branch_len i;

% Identify zero injection buses
ZI_buses = zeros(num_buses, 1);
% Find all buses with generators
gen_buses = data.gen(:, 1);
for i=1:1:num_buses
    % Find all buses with no load and no generator
    if data.bus(i, 3) == 0 && data.bus(i, 4) == 0 && ~any(gen_buses == i)
        ZI_buses(i) = i;
    end
end
ZI_buses = nonzeros(ZI_buses);
clear gen_buses;

%% Begin Solving the OPP
fprintf('\nSolver starts work...\n');

eval(sprintf('[p, m] = OPP_%s(A, ZI_buses);', alg));
fprintf('\n=================================================================\n')
fprintf('Test Case: %s', test_case)
fprintf('\n=================================================================\n')
fprintf('\nZero Injection Buses:');
fprintf(' %d', ZI_buses);
fprintf('\n');

fprintf('Method: %s\n\n', m.method);
fprintf('No. of PMUs: %d\nPlacement (bus no.):\n', length(p));
for i=1:1:size(p,1)
    fprintf('  ')
    fprintf('%d ', p(i,:));
    fprintf('\n')
end
fprintf('\nTime taken:%fs\n\n', m.time)
clear p m alg test_case;

%% Clear global system variables
clear data num_buses A ZI_buses;
