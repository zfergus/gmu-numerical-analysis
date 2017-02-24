% MATH 446: Project 05
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 05\nWritten by Zachary Ferguson\n\n');

    fprintf('Gaussian Elimination:\n\n');

    % Q1a
    A = [2 -2 -1 ; 4 1 -2 ; -2 1 -1];
    b = [-2 ; 1 ; -3];
    x = gaussian_elimination(A, b);
    print('Q1a:\n\tA = %s\n\tb = %s\n\tx = %s\n\n', A, b, x);
    print('\tAx = %s\n', A*x);

    % Q1b
    A = [1 2 -1 ; 0 3 1 ; 2 -1 1];
    b = [2 ; 4 ; 2];
    x = gaussian_elimination(A, b);
    print('Q1b:\n\tA = %s\n\tb = %s\n\tx = %s\n\n', A, b, x);
    print('\tAx = %s\n', A*x);

    % Q1c
    A = [2 1 -4 ; 1 -1 1 ; -1 3 -2];
    b = [-7 ; -2 ; 6];
    x = gaussian_elimination(A, b);
    print('Q1c:\n\tA = %s\n\tb = %s\n\tx = %s\n\n', A, b, x);
    print('\tAx = %s\n', A*x);

    % Q2a
    n = 2;
    H = hilbert_matrix(n);
    b = ones(n);
    x = gaussian_elimination(H, b);
    print('Q2a:\n\tn = %d\n\tH = %s\n\tb = %s\n\tx = %s\n\n', n, H, b, x);
    print('\tHx = %s\n', H*x);

    % Q2b
    n = 5;
    H = hilbert_matrix(n);
    b = ones(n);
    x = gaussian_elimination(H, b);
    print('Q2b:\n\tn = %d\n\tH = %s\n\tb = %s\n\tx = %s\n\n', n, H, b, x);
    print('\tHx = %s\n', H*x);

    % Q2a
    n = 10;
    H = hilbert_matrix(n);
    b = ones(n);
    x = gaussian_elimination(H, b);
    print('Q2c:\n\tn = %d\n\tH = %s\n\tb = %s\n\tx = %s\n\n', n, H, b, x);
    print('\tHx = %s\n', H*x);

    fprintf('LU Decomposition:\n\n');

    % Q1a
    A = [3 1 2 ; 6 3 4 ; 3 1 5];
    [L, U] = lu_decomposition(A);
    print('Q1a:\n\tA = %s\n\tL = %s\n\tU = %s\n\n', A, L, U);
    print('\tLU = %s\n', L*U);

    % Q1b
    A = [4 2 0 ; 4 4 2 ; 2 2 3];
    [L, U] = lu_decomposition(A);
    print('Q1b:\n\tA = %s\n\tL = %s\n\tU = %s\n\n', A, L, U);
    print('\tLU = %s\n', L*U);

    % Q1c
    A = [1 -1 1 2 ; 0 2 1 0 ; 1 3 4 4 ; 0 2 1 -1];
    [L, U] = lu_decomposition(A);
    print('Q1c:\n\tA = %s\n\tL = %s\n\tU = %s\n\n', A, L, U);
    print('\tLU = %s\n', L*U);

    % Q2a
    A = [3 1 2 ; 6 3 4 ; 3 1 5];
    b = [0 ; 1 ; 3];
    [L, U] = lu_decomposition(A);
    x = lu_solve(L, U, b);
    print('Q2a:\n\tA = %s\n\tb = %s\n\tL = %s\n\tU = %s\n\tx = %s\n\t\n\n', ...
        A, b, L, U, x);
    print('\tLU = %s\n', L*U);
    print('\tAx = %s\n', A*x);
    print('\tLUx = %s\n', L*U*x);

    % Q2b
    A = [4 2 0 ; 4 4 2 ; 2 2 3];
    b = [2 ; 4 ; 6];
    [L, U] = lu_decomposition(A);
    x = lu_solve(L, U, b);
    print('Q2b:\n\tA = %s\n\tb = %s\n\tL = %s\n\tU = %s\n\tx = %s\n\t\n\n', ...
        A, b, L, U, x);
    print('\tLU = %s\n', L*U);
    print('\tAx = %s\n', A*x);
    print('\tLUx = %s\n', L*U*x);
end
