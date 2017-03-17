% Solve a system of linear equations, Ax = b, using the Jacobi Method
% Written by Zachary Ferguson

function xc = jacobi_method(A, b, x0, eps)
    % Solve the equation Ax = b using the Jacobi Method
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

    n = size(A, 1);
    D = spdiags(spdiags(A, 0), 0, n, n);
    L_plus_U = A - D; % L + U = A - D
    D_inv = spdiags(spdiags(A,0).^-1, 0, n, n); % = D^-1

    x_prev = x0;
    xc = D_inv*(b - L_plus_U * x0);
    n_steps = 1;
    while (norm(xc - x_prev, inf) >= 0.5 * eps)
        x_prev = xc;
        xc = D_inv*(b - (L_plus_U) * xc);
        n_steps = n_steps + 1;
    end
    fprintf('\tNumber of steps to find solution: %d\n', n_steps);
end
