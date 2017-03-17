% Solve a system of linear equations, Ax = b, using the Guass-Seidel Method
% Written by Zachary Ferguson

function xc = gauss_seidel_method(A, b, x0, eps)
    % Solve the equation Ax = b using the Guass-Seidel Method
    % Input:
    %   A - matrix of coefficients to the linear equations
    %   b - Right hand side of the linear equations
    %   x0 - intial guess for solution vector
    %   eps - tolerance of forward error
    % Output:
    %   xc - computed solution to a eps tolerance
    if nargin < 4
        eps = 1e-6;
    end

    U = triu(A, 1);
    L_plus_D_inv = (A - U)^-1; % (L + D = A - U)^-1

    x_prev = x0;
    xc = L_plus_D_inv*(b - U * x0);
    n_steps = 1;
    while (norm(xc - x_prev, inf) >= 0.5 * eps)
        x_prev = xc;
        xc = L_plus_D_inv*(b - U * xc);
        n_steps = n_steps + 1;
    end
    fprintf('\tNumber of steps to find solution: %d\n', n_steps);
end
