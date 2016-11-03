%% Description
%
% This file will run test for all the five systems using all the six
% algorithms for 20 times.
% It may take quite a long time to finish all tests.
% 
% Run main.m for single test.
%

%% Set up output mode
disp('Please Choose output mode:');
disp('1. Command-line');
disp('2. File');
while true
    c = input('Please input 1 or 2: ');
    if c == 1
        mode = 1;
        break;
    elseif c == 2
        mode = 2;
        break;
    end
end

if mode == 2
    c = input('Please input a file name for the output file: ', 's');
    [file, ~] = fopen(c, 'w');
    fprintf(file, 'DFS: Depth First Search');
    fprintf(file, '\nGThM: Graphic Theoretic Procedure using Merger Method');
    fprintf(file, '\nGThN: Graphic Theoretic Procedure using Nonlinear Constraint Function Method');
    fprintf(file, '\nSA: Simulated Annealing Method (Original Version proposed in 1993)');
    fprintf(file, '\nSAB: Simulated Annealing Method (Modified Version proposed in the report)');
    fprintf(file, '\nRSN: Recursive Security N Algorithm');
    fprintf(file, '\n\nTest System\tMethod Name\tNo. of PMUs\tTime Taken (second)\tPlacement (Bus No.)\n');
end

clear c;

%% Simulation Begins
disp('');
disp('Simulation Begins...')

cases = [14, 30, 39, 57, 118];
cases_name = {'IEEE 14-bus System', 'IEEE 30-bus System', 'IEEE 39-bus System', 'New England 57-bus System', 'IEEE 118-bus System'};
cases_name_simple = {'14-bus', '30-bus', '39-bus', '57-bus', '118-bus'};
algs = {'DFS', 'GThM', 'GThN', 'RSN', 'SA', 'SAB', 'RSN'};
for cas=1:1:5
    data = eval(sprintf('case%d()', cases(cas)));
    test_case = cases_name{cas};
    
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
    
    % Begin Solving the OPP
    fprintf('\n=================================================================\n')
    fprintf('Test Case: %s', test_case)
    fprintf('\n=================================================================\n')
    fprintf('\nZero Injection Buses:');
    fprintf(' %d', ZI_buses);
    fprintf('\n');
        
    for alg=1:1:3
        for test=1:1:20
            eval(sprintf('[p, m] = OPP_%s(A, ZI_buses);', algs{alg}));
            if mode == 1
                if test == 1
                    fprintf('Method: %s\n\n', m.method);
                end
                fprintf('  Test Case: %d\n', test);
                fprintf('    No. of PMUs: %d\n    Placement (bus no.):\n', length(p(1,:)));
                for i=1:1:size(p,1)
                    fprintf('      ');
                    fprintf('%d ', p(i,:));
                    fprintf('\n')
                end
                fprintf('    Time taken:%fs\n\n', m.time);
            else
                if test == 1
                    fprintf('\nMethod: %s\n\n', m.method);
                end
                fprintf('  Test Case %d Finished\n', test);
                placement = cell(1, size(p,1));
                for i=1:1:size(p,1)
                    placement{i} = strcat(sprintf('%d,', p(i,1:length(p(i,:))-1)), sprintf('%d', p(i,length(p(i,:)))));
                end
                fprintf(file, '%s\t%s\t%d\t%f\t%s\n', cases_name_simple{cas}, algs{alg}, length(p(i,:)),  m.time, strjoin(placement, ' / '));
                
            end
                
            clear p m
        end
    end
    clear alg test_case;

    % Clear system variables
    clear data num_buses A ZI_buses;
 end
fclose(file);
clear cas cases cases_name_simple algs mode file