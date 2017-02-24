% MATH 446: Project 05
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 05\nWritten by Zachary Ferguson\n\n');

    fprintf('Gaussian Elimination:\n\n');

    % Q1a
    A = [2 -2 -1 ; 4 1 -2 ; -2 1 -1];
    b = [-2 ; 1 ; -3];
    x = gaussian_elimination(A, b);
    fprintf('Q1a:\n');
    print_Axb(A, x, b);
    
    % Q1b
    A = [1 2 -1 ; 0 3 1 ; 2 -1 1];
    b = [2 ; 4 ; 2];
    x = gaussian_elimination(A, b);
    fprintf('Q1b:\n');
    print_Axb(A, x, b);
    
    % Q1c
    A = [2 1 -4 ; 1 -1 1 ; -1 3 -2];
    b = [-7 ; -2 ; 6];
    x = gaussian_elimination(A, b);
    fprintf('Q1c:\n');
    print_Axb(A, x, b);

    % Q2a
    n = 2;
    H = hilbert_matrix(n);
    b = ones(n, 1);
    x = gaussian_elimination(H, b);
    fprintf('Q2a:\n\tn = %d\n', n);
    print_Axb(H, x, b, 'H');

    % Q2b
    n = 5;
    H = hilbert_matrix(n);
    b = ones(n, 1);
    x = gaussian_elimination(H, b);
    fprintf('Q2b:\n\tn = %d\n', n);
    print_Axb(H, x, b, 'H');

    % Q2a
    n = 10;
    H = hilbert_matrix(n);
    b = ones(n, 1);
    x = gaussian_elimination(H, b, 1e-10);
    fprintf('Q2c:\n\tn = %d\n', n);
    print_Axb(H, x, b, 'H');

    fprintf('\nLU Decomposition:\n\n');

    % Q1a
    A = [3 1 2 ; 6 3 4 ; 3 1 5];
    [L, U] = lu_decomposition(A);
    fprintf('Q1a:\n');
    print_ALU(A, L, U);

    % Q1b
    A = [4 2 0 ; 4 4 2 ; 2 2 3];
    [L, U] = lu_decomposition(A);
    fprintf('Q1b:\n');
    print_ALU(A, L, U);

    % Q1c
    A = [1 -1 1 2 ; 0 2 1 0 ; 1 3 4 4 ; 0 2 1 -1];
    [L, U] = lu_decomposition(A);
    fprintf('Q1c:\n');
    print_ALU(A, L, U);

    % Q2a
    A = [3 1 2 ; 6 3 4 ; 3 1 5];
    b = [0 ; 1 ; 3];
    [L, U] = lu_decomposition(A);
    x = lu_solve(L, U, b);
    fprintf('Q2a:\n');
    print_AbLUx(A, b, L, U, x);

    % Q2b
    A = [4 2 0 ; 4 4 2 ; 2 2 3];
    b = [2 ; 4 ; 6];
    [L, U] = lu_decomposition(A);
    x = lu_solve(L, U, b);
    fprintf('Q2b:\n');
    print_AbLUx(A, b, L, U, x);
end


function print_Axb(A, x, b, nameA)
    if nargin < 4
        nameA = 'A';
    end
    
    fprintf('\t%s = \n', nameA);
	disp(A);
    fprintf('\tb = \n');
    disp(b);
    fprintf('\tx = \n');
    disp(x);
    fprintf('\tBackwards Error = ||%sx - b|| = %.10f\n\n', nameA, ...
        max(abs(A*x - b)));
end

function print_ALU(A, L, U, nameA)
    if nargin < 4
        nameA = 'A';
    end
    
    fprintf('\t%s =\n', nameA);
    disp(A);
    fprintf('\tL =\n');
    disp(L);
    fprintf('\tU =\n');
    disp(U);
    diff = abs(A - L*U);
    fprintf('\tBackwards Error = ||%s - LU|| = %.10f\n\n', nameA, ...
        max(diff(:)));
end

function print_AbLUx(A, b, L, U, x, nameA)
    if nargin < 6
        nameA = 'A';
    end
  
    fprintf('\t%s =\n', nameA);
    disp(A);
    fprintf('\tb =\n');
    disp(b);
    fprintf('\tL =\n');
    disp(L);
    fprintf('\tU =\n');
    disp(U);
    fprintf('\tx =\n');
    disp(x);
    fprintf('\tBackwards Error = ||%sx - b|| = %.10f\n\n', nameA, ...
        max(abs(A*x - b)));
end
