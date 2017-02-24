% Generates the n x n Hilbert matrix where H(i, j) = 1 / (i+j-1)
% Written by Zachary Ferguson

function H = hilbert_matrix(n)
    % Generates the n x n Hilbert matrix where H(i, j) = 1 / (i+j-1)
    % Input:
    %   n - size of matrix
    % Output:
    %   H - nxn Hilbert matrix
    H = zeros(n, n);

    for i = 1 : n
        for j = 1 : n
            H(i, j) = 1 / (i + j - 1);
        end
    end
end
