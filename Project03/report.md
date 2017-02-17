# MATH 446: Project 03
### Zachary Ferguson

## Contents

1. Questions
2. Code
    1. Newton's Method
    2. Main
3. Output

## Questions

### Q1.
#### a.
\begin{equation}
\begin{matrix}
f(x) = x^3 - 2x - 2 = 0 \\
r = 1.76929235
\end{matrix}
\end{equation}

#### b.
\begin{equation}
\begin{matrix}
f(x) = e^x + x - 7 = 0 \\
r = 1.67282169
\end{matrix}
\end{equation}

#### c.
\begin{equation}
\begin{matrix}
f(x) = e^x + sin(x) - 4 = 0 \\
r = 1.12998049
\end{matrix}
\end{equation}

### Q3.
#### a.
\begin{equation}
f(x) = 27x^3 + 54x^2 + 36x + 8 = 0
\end{equation}
\begin{equation}
f'(x) = 81x^2 + 108x + 36
\end{equation}
\begin{equation}
f''(x) = 162x + 108
\end{equation}
\begin{equation}
f'''(x) = 162 \neq 0
\end{equation}
\begin{equation}
\begin{matrix}
r = -{2 \over 3} \\
f(r) = f'(r) = f''(r) = 0 \rightarrow \text{multiplicity of } r \text{ is } 3
\end{matrix}
\end{equation}

### Q9.
\begin{equation}
f(x)  = 14xe^{x-2} - 12e^{x-2} - 7x^3 + 20x^2 - 26x + 12
\end{equation}
\begin{equation}
f'(x) = 14xe^{x-2} + 2e^{x-2} - 21x^2 + 40x - 26
\end{equation}
\begin{equation}
f''(x) = 14xe^{x-2} + 16e^{x-2} - 42x + 40
\end{equation}
\begin{equation}
f'''(x) = 14xe^{x-2} + 30e^{x-2} - 42
\end{equation}

\begin{equation}
\begin{matrix}
r_1 = 0.85714285 \\
f'(r1) \approx -2.67817 \neq 0 \rightarrow M = lim_{i \rightarrow \infty}
    {e_{i+1} \over {e_i}^2} = \left| {f''(r_1) \over {2f'(r_1)}} \right|
    \approx 1.69939188
\end{matrix}
\end{equation}

\begin{equation}
\begin{matrix}
r_2 = 2.0 \\
f(r_2) = f'(r_2) = f''(r_2) = 0 \rightarrow \text{multiplicity of } r_2
    \text{ is } 3 \\
\therefore S = lim_{i \rightarrow \infty} {e_{i+1} \over e_i} = {m-1 \over m}
    = {2 \over 3}
\end{matrix}
\end{equation}

## Code

### Newton's Method

```matlab
% Computes the roots of a function using the Newton's Method.
% Written by Zachary Ferguson

function xc = newtons_method(f, fp, x0, tol, m, print_ei)
    % Compute the root to f(x) using Newton's Method
    % Input:
    %   f - function to find the roots of
    %   fp - first dirivative of f(x)
    %   x0 - intial guess
    %   tol - tolerance for the root
    %   m - multiplicity of the root
    %   print_ei - which e_i limit should be printed
    % Output:
    %   xc - computed root to the function f(x).
    if nargin < 4
        tol = 1e-9;
    end
    if nargin < 5
        m = 1;
    end
    if nargin < 6
        print_ei = 0;
    end

    r = fzero(f, x0);

    n = 0;
    x = x0;
    ei = 1;
    ei_1 = 1;
    while (abs(f(x)) >= 0.5 * tol)
        x = x - m * f(x) / fp(x);
        n = n + 1;
        ei_1 = ei;
        ei = abs(r - x);
        if print_ei == 1
            fprintf('\te_i = %.8f; e_(i+1)/e_i = %.8f\n', ei, ei/ei_1);
        elseif print_ei == 2
            fprintf('\te_i = %.8f; e_(i+1)/(e_i)^2 = %.8f\n', ei, ei/(ei_1^2));
        end
    end
    fprintf('\tn = %d\n', n)
    xc = x;
end
```

### Main

```matlab
% MATH 446: Project 03
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 03\nWritten by Zachary Ferguson\n\n');

    % Q1a
    fprintf('Q1a:\n\tf(x) = x^3 - 2x - 2 = 0\n');
    f  = @(x) x^3 - 2*x - 2;
    fp = @(x) 3*x^2 - 2;
    x0 = 2;
    fprintf('\tx0 = %g\n', x0);
    fprintf('\tr = %.10f\n', newtons_method(f, fp, x0));

    % Q1b
    fprintf('Q1b:\n\tf(x) = e^x + x - 7 = 0\n');
    f  = @(x) exp(x) + x - 7;
    fp = @(x) exp(x) + 1;
    x0 = 0;
    fprintf('\tx0 = %g\n', x0);
    fprintf('\tr = %.10f\n', newtons_method(f, fp, x0));

    % Q1c
    fprintf('Q1c:\n\tf(x) = e^x + sin(x) - 4 = 0\n');
    f  = @(x) exp(x) + sin(x) - 4;
    fp = @(x) exp(x) + cos(x);
    x0 = 2;
    fprintf('\tx0 = %g\n', x0);
    fprintf('\tr = %.10f\n', newtons_method(f, fp, x0));

    % Q3a
    fprintf('Q3a:\n\tf(x) = 27x^3 + 54x^2 + 36x + 8 = 0\n');
    f  = @(x) 27*x^3 + 54*x^2 + 36*x + 8;
    fp = @(x) 81*x^2 + 108*x + 36;
    x0 = 0.0;
    r = -2/3;
    fprintf('\tx0 = %g\n', x0);
    xc = newtons_method(f, fp, x0, 1e-16);
    fprintf('\txc = %.16f\n', xc);
    fprintf('\tForward Error = |r - xc| = %.16f\n', abs(r-xc));
    fprintf('\tBackward Error = f(xc) = %.16f\n', f(xc));
    fprintf('\tmultiplicity of r is 3\n');
    xc = newtons_method(f, fp, x0, 1e-16, 3);
    fprintf('\txc = %.16f\n', xc);
    fprintf('\tForward Error = |r - xc| = %.16f\n', abs(r-xc));
    fprintf('\tBackward Error = f(xc) = %.16f\n', f(xc));

    % Q9
    fprintf('Q9:\n\tf(x) = 14xe^(x-2) - 12e^(x-2) - 7x^3 + 20x^2 - 26x + 12\n');
    f  = @(x) 14*x*exp(x-2) - 12*exp(x-2) - 7*x^3 + 20*x^2 - 26*x + 12;
    fp = @(x) 14*x*exp(x-2) + 2*exp(x-2) - 21*x^2 + 40*x - 26;
    fpp = @(x) 14*x*exp(x-2) + 4*exp(x-2) - 42*x + 40;
    x0 = 0;
    r = newtons_method(f, fp, x0, 1e-9, 1, 2);
    fprintf('\tr1 = %.10f\n', r);
    fprintf('\tM = lim i->inf (e_(i+1)/(e_i)^2) = %.10f\n\n', ...
        abs(fpp(r)/(2*fp(r))));

    x0 = 3.0;
    r = newtons_method(f, fp, x0, 1e-9, 1, 1);
    fprintf('\tr2 = %.10f\n', r);
    fprintf('\tmultiplicity of r2 is 3 -> ');
    fprintf('S = lim i->inf (e_(i+1)/e_i) = %.10f\n', 2/3);
end
```

## Output

```
MATH 446: Project 03
Written by Zachary Ferguson

Q1a:
	f(x) = x^3 - 2x - 2 = 0
	x0 = 2
	n = 4
	r = 1.7692923542
Q1b:
	f(x) = e^x + x - 7 = 0
	x0 = 0
	n = 7
	r = 1.6728216986
Q1c:
	f(x) = e^x + sin(x) - 4 = 0
	x0 = 2
	n = 5
	r = 1.1299804987
Q3a:
	f(x) = 27x^3 + 54x^2 + 36x + 8 = 0
	x0 = 0
	n = 31
	xc = -0.6666638081419596
	Forward Error = |r - xc| = 0.0000028585247071
	Backward Error = f(xc) = 0.0000000000000000
	multiplicity of r is 3
	n = 1
	xc = -0.6666666666666666
	Forward Error = |r - xc| = 0.0000000000000000
	Backward Error = f(xc) = 0.0000000000000000
Q9:
	f(x) = 14xe^(x-2) - 12e^(x-2) - 7x^3 + 20x^2 - 26x + 12
	e_i = 0.45386858; e_(i+1)/(e_i)^2 = 0.45386858
	e_i = 0.19642053; e_(i+1)/(e_i)^2 = 0.95351302
	e_i = 0.05608698; e_(i+1)/(e_i)^2 = 1.45374528
	e_i = 0.00639065; e_(i+1)/(e_i)^2 = 2.03151817
	e_i = 0.00009651; e_(i+1)/(e_i)^2 = 2.36321178
	e_i = 0.00000002; e_(i+1)/(e_i)^2 = 2.41307008
	e_i = 0.00000000; e_(i+1)/(e_i)^2 = 1.75788724
	n = 7
	r1 = 0.8571428571
	M = lim i->inf (e_(i+1)/(e_i)^2) = 1.6993918897

	e_i = 0.73383688; e_(i+1)/e_i = 0.73383688
	e_i = 0.52975850; e_(i+1)/e_i = 0.72190225
	e_i = 0.37665244; e_(i+1)/e_i = 0.71098896
	e_i = 0.26413220; e_(i+1)/e_i = 0.70126241
	e_i = 0.18301754; e_(i+1)/e_i = 0.69290129
	e_i = 0.12555240; e_(i+1)/e_i = 0.68601292
	e_i = 0.08544813; e_(i+1)/e_i = 0.68057740
	e_i = 0.05780157; e_(i+1)/e_i = 0.67645219
	e_i = 0.03892465; e_(i+1)/e_i = 0.67341853
	e_i = 0.02612762; e_(i+1)/e_i = 0.67123586
	e_i = 0.01749716; e_(i+1)/e_i = 0.66968057
	e_i = 0.01169797; e_(i+1)/e_i = 0.66856388
	e_i = 0.00781113; e_(i+1)/e_i = 0.66773373
	e_i = 0.00521055; e_(i+1)/e_i = 0.66706772
	e_i = 0.00347263; e_(i+1)/e_i = 0.66646198
	e_i = 0.00231214; e_(i+1)/e_i = 0.66581839
	e_i = 0.00153764; e_(i+1)/e_i = 0.66502955
	e_i = 0.00102094; e_(i+1)/e_i = 0.66396188
	e_i = 0.00067630; e_(i+1)/e_i = 0.66242890
	e_i = 0.00044646; e_(i+1)/e_i = 0.66015724
	n = 20
	r2 = 2.0004598418
	multiplicity of r2 is 3 -> S = lim i->inf (e_(i+1)/e_i) = 0.6666666667
```
