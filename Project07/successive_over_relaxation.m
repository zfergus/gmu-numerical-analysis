% Solve a system of linear equations, Ax = b, using the Successive Over
% Relaxation.
% Written by Zachary Ferguson

function xc = successive_over_relaxation(A, b, x0, omega, eps)
    % Solve the equation Ax = b using the Successive Over Relaxation
    % Input:
    %   A - matrix of coefficients to the linear equations
    %   b - Right hand side of the linear equations
    %   x0 - intial guess for solution vector
    %   ω - parameter for how much to relax
    %   ε - tolerance of forward error
    % Output:
    %   xc - computed solution to a eps tolerance
    if nargin < 5
        eps = 1e-6;
    end

    % A = L + D + U
    U = triu(A, 1);
    L = tril(A, -1);
    D = A - L - U;
    omegaL_plus_D_inv = (omega*L + D)^-1; % (ωL + D)⁻¹
    omega_rhs = omega*omegaL_plus_D_inv*b; % ω(ωL + D)⁻¹b

    x_prev = x0;
    % x_{k+1} = (ωL + D)⁻¹[(1-ω)Dx_k - ωUx_k]+ω(ωL + D)⁻¹b
    xc = omegaL_plus_D_inv*((1-omega)*D*x0 - omega*U*x0) + omega_rhs;
    n_steps = 1;
    while (norm(xc - x_prev, inf) >= 0.5 * eps)
        x_prev = xc;
        xc = omegaL_plus_D_inv*((1-omega)*D*xc - omega*U*xc) + omega_rhs;
        n_steps = n_steps + 1;
    end
    fprintf('\tNumber of steps to find solution: %d\n', n_steps);
end
