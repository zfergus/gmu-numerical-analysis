% MATH 446: Project 06
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 06\nWritten by Zachary Ferguson\n\n');

    % Q1a
    n = 6;
    fprinf('Q1a:\n');
    question1(n);

    % Q1b
    n = 10;
    fprinf('Q1b:\n');
    question1(n);
end


function A = build_matrix(n)
    % Generates the nxn matrix where A(i, j) = 5 / (i+2j-1)
    % Input:
    %   n - size of matrix
    % Output:
    %   A - nxn Hilbert matrix
    A = zeros(n, n);

    for i = 1 : n
        for j = 1 : n
            A(i, j) = 5 / (i + 2*j - 1);
        end
    end
end


function question1(n)
    % Prints out appropriate information to answer question 1.
    fprinf('\tn = %d\n', n);
    A = build_matrix(n);
    x = ones(n);
    b = A * x;
    fprintf('\tA =\n'); disp(A);
    fprintf('\tx =\n'); disp(x);
    fprintf('\tb = Ax =\n'); disp(b);
    xc = A \ b; % Solve for xc
    fprintf('\txc = \n'); disp(xc);
    BE = norm(b - A*xc, inf); % infiniry norm
    FE = norm(x - xc, inf);
    RBE = BE / norm(b, inf);
    RFE = FE / norm(x, inf);
    EMF = RFE / RBE;
    condA = cond(A, inf); % Condition number of A
    fprintf('Backwards Error = %g\n\tForwards Error = %g\n
        \tRelative BE = %g\n\tRelative FE = %g\n
        \tError Magnification Factor = %g\n\tcond(A) = %g', ...
        BE, FE, RBE, RFE, EMF, condA);
end
