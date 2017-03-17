% MATH 446: Project 07
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 07\nWritten by Zachary Ferguson\n\n');

    fprintf('Jacobi Method:\n\n')
    n = 100;
    fprintf('Q1a:\n\tn=%d\n', n);
    [A, b] = build_system(n);
    x = ones(n, 1);
    xc = jacobi_method(A, b, zeros(n, 1));
    print_errors(A, b, x, xc);

    m = 100000;
    fprintf('Q1b:\n\tn=%d\n', m);
    [C, d] = build_system(m);
    y = ones(m, 1);
    yc = jacobi_method(C, d, zeros(m, 1));
    print_errors(C, d, y, yc);

    fprintf('\nGauss-Seidel Method:\n\n');
    fprintf('Q5a:\n\tn=%d\n', n);
    xc = gauss_seidel_method(A, b, zeros(n, 1));
    print_errors(A, b, x, xc);

    fprintf('\nSuccessive Over Relaxation:\n\n');
    fprintf('Q5b:\n\tn=%d\n', n);
    xc = successive_over_relaxation(A, b, zeros(n, 1), 1.2);
    print_errors(A, b, x, xc);
end

function [A, b] = build_system(n)
    diag_elements = [[-1 * ones(n-1, 1); 0], 3*ones(n, 1), ...
        [[0; -1 * ones(n-1, 1)]]];
    diag_indices = [-1; 0; 1];
    A = spdiags(diag_elements, diag_indices, n, n);
    b = [2; ones(n-2, 1); 2];
end

function print_errors(A, b, x, xc)
    BE = norm(b - A*xc, inf); % infiniry norm
    FE = norm(x - xc, inf);

    fprintf('\tBackwards Error = %g\n', BE);
    fprintf('\tForwards Error = %g\n', FE);
end
