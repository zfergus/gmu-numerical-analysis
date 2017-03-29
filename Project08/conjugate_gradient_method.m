% Solve a system of linear equations, Ax = b, using the Conjugate Gradient
% Method
% Written by Zachary Ferguson

function xc = conjugate_gradient_method(A, b, x0, M, disp_figure, disp_prop)
    % Solve the equation Ax = b using the Conjugate Gradient Method
    % Input:
    %   A - matrix of coefficients to the linear equations
    %   b - Right hand side of the linear equations
    %   x0 - intial guess for solution vector
    %   M - optional preconditioner matrix (Default: Identity)
    % Output:
    %   xc - computed solution to a eps tolerance
    n = size(A, 1);
    if nargin < 4
        M = speye(n);
    end
    if nargin < 5
        disp_figure = false;
    end

    xc = x0; % Computed solution
    r = b - A*x0; % Residual
    d = M \ r; % Direction
    z = d;
    errors = [];
    for k = 0 : (n-1)
        val_r = norm(r, inf);
        errors = [errors val_r];
        if val_r <= 1e-16
            break
        end
        Ad = A * d;
        alpha = (r' * z) / (d' * Ad);
        xc = xc + alpha * d;
        r_prev = r;
        r = r - alpha * Ad;
        z_prev = z;
        z = M \ r; % M^-1 * r
        beta = (r' * z) / (r_prev' * z_prev);
        d = z + beta * d;
    end
    if disp_figure
        semilogy(1:size(errors, 2), errors, disp_prop)
    end
    fprintf('\tNumber of steps to find solution: %d\n', k+1);
end
