% Solve a system of linear equations, Ax = b, using naive Guassian Elimination
% Written by Zachary Ferguson

function x = gaussian_elimination(A, b, eps)
    % Solve the equation Ax = b
    % Input:
    %   A - matrix of coefficients to the linear equations
    %   b - Right hand side of the linear equations
    %   eps - tolerance of a zero pivot
    % Output:
    %   x - solved value
    if nargin < 3
            eps = 1e-9;
    end

    n = size(A, 1);
    
    % Elimination step (O(2/3 * n^3))
    for j = 1 : n-1
        if abs(A(j, j)) < eps
            error('Zero Pivot encountered.');
        end
        for i = j+1 : n
            mult = A(i, j)/A(j, j);
            for k = j+1 : n
                A(i, k) = A(i, k) - mult * A(j, k); % Row operation
            end
            b(i) = b(i) - mult * b(j);
        end
    end

    x = zeros(size(b));

    % Perform Back Substitution
    for i = n : -1 : 1
        for j = i+1 : n
            b(i) = b(i) - A(i, j) * x(j);
        end
        x(i) = b(i) / A(i, i);
    end
end
