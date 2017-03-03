# MATH 446: Project 06
### Zachary Ferguson

## Contents

1. Questions
2. Code
    1. Main
3. Output

## Questions

### Question 5

Values of n greater than 12 result in a solution with no significant digits.
The Condition Number of the 12x12 matrix is ~2.2e17. This implies even a
machine precision error of 1e-16 will be magnified to affect all digits.
Notably, the A matrix of size 12x12 is singular to the machine precision.

## Code

### Main

```matlab
% MATH 446: Project 06
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 06\nWritten by Zachary Ferguson\n\n');

    % Q1a
    n = 6;
    fprintf('Q1a:\n');
    question1(n);

    % Q1b
    n = 10;
    fprintf('\nQ1b:\n');
    question1(n);

    % Q5
    n = 12;
    fprintf('\nQ5:\n');
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
    fprintf('\tn = %d\n', n);
    A = build_matrix(n);
    x = ones(n, 1);
    b = A * x;

    fprintf('\tA =\n'); disp(A);
    fprintf('\n\tx =\n'); disp(x);
    fprintf('\n\tb = Ax =\n'); disp(b);

    xc = A \ b; % Solve for xc
    fprintf('\n\txc = \n'); disp(xc);

    BE = norm(b - A*xc, inf); % infiniry norm
    FE = norm(x - xc, inf);
    RBE = BE / norm(b, inf);
    RFE = FE / norm(x, inf);
    EMF = RFE / RBE;
    condA = cond(A, inf); % Condition number of A

    fprintf('\n\tBackwards Error = %g\n', BE);
    fprintf('\tForwards Error = %g\n', FE);
    fprintf('\tRelative BE = %g\n', RBE);
    fprintf('\tRelative FE = %g\n', RFE);
    fprintf('\tError Magnification Factor = %g\n', EMF);
    fprintf('\tcond(A) = %g\n', condA);
end
```

## Output

```
MATH 446: Project 06
Written by Zachary Ferguson

Q1a:
        n = 6
        A =
   2.50000   1.25000   0.83333   0.62500   0.50000   0.41667
   1.66667   1.00000   0.71429   0.55556   0.45455   0.38462
   1.25000   0.83333   0.62500   0.50000   0.41667   0.35714
   1.00000   0.71429   0.55556   0.45455   0.38462   0.33333
   0.83333   0.62500   0.50000   0.41667   0.35714   0.31250
   0.71429   0.55556   0.45455   0.38462   0.33333   0.29412

        x =
   1
   1
   1
   1
   1
   1

        b = Ax =
   6.1250
   4.7757
   3.9821
   3.4423
   3.0446
   2.7365

        xc =
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000

        Backwards Error = 8.88178e-16
        Forwards Error = 3.22728e-10
        Relative BE = 1.45009e-16
        Relative FE = 3.22728e-10
        Error Magnification Factor = 2.22557e+06
        cond(A) = 7.0342e+07

Q1b:
        n = 10
        A =
   2.50000   1.25000   0.83333   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000
   1.66667   1.00000   0.71429   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810
   1.25000   0.83333   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727
   1.00000   0.71429   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739
   0.83333   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727   0.20833
   0.71429   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739   0.20000
   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727   0.20833   0.19231
   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739   0.20000   0.18519
   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727   0.20833   0.19231   0.17857
   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739   0.20000   0.18519   0.17241

        x =
   1
   1
   1
   1
   1
   1
   1
   1
   1
   1

        b = Ax =
   7.3224
   5.9044
   5.0497
   4.4551
   4.0080
   3.6551
   3.3670
   3.1260
   2.9206
   2.7429

        xc =
   1.00000
   1.00000
   0.99999
   1.00006
   0.99961
   1.00133
   0.99743
   1.00285
   0.99832
   1.00041

        Backwards Error = 8.88178e-16
        Forwards Error = 0.00284718
        Relative BE = 1.21296e-16
        Relative FE = 0.00284718
        Error Magnification Factor = 2.34731e+13
        cond(A) = 1.31346e+14

Q5:
        n = 12
        A =
 Columns 1 through 10:

   2.50000   1.25000   0.83333   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000
   1.66667   1.00000   0.71429   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810
   1.25000   0.83333   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727
   1.00000   0.71429   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739
   0.83333   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727   0.20833
   0.71429   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739   0.20000
   0.62500   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727   0.20833   0.19231
   0.55556   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739   0.20000   0.18519
   0.50000   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727   0.20833   0.19231   0.17857
   0.45455   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739   0.20000   0.18519   0.17241
   0.41667   0.35714   0.31250   0.27778   0.25000   0.22727   0.20833   0.19231   0.17857   0.16667
   0.38462   0.33333   0.29412   0.26316   0.23810   0.21739   0.20000   0.18519   0.17241   0.16129

 Columns 11 and 12:

   0.22727   0.20833
   0.21739   0.20000
   0.20833   0.19231
   0.20000   0.18519
   0.19231   0.17857
   0.18519   0.17241
   0.17857   0.16667
   0.17241   0.16129
   0.16667   0.15625
   0.16129   0.15152
   0.15625   0.14706
   0.15152   0.14286

        x =
   1
   1
   1
   1
   1
   1
   1
   1
   1
   1
   1
   1

        b = Ax =
   7.7580
   6.3218
   5.4503
   4.8403
   4.3789
   4.0127
   3.7122
   3.4597
   3.2435
   3.0557
   2.8905
   2.7440
warning: matrix singular to machine precision, rcond = 3.70333e-18

        xc =
   1.00000
   1.00000
   0.99998
   1.00031
   0.99774
   1.00872
   0.98170
   1.01833
   0.99955
   0.98278
   1.01513
   0.99575
warning: inverse: matrix singular to machine precision, rcond = 3.70333e-18

        Backwards Error = 8.88178e-15
        Forwards Error = 0.0183309
        Relative BE = 1.14485e-15
        Relative FE = 0.0183309
        Error Magnification Factor = 1.60116e+13
        cond(A) = 2.16381e+17
```
