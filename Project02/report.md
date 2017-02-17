# MATH 446: Project 02
#### Zachary Ferguson

## Contents

1. Questions
2. Code
	1. Fixed Point Iterative Method
	2. Main
3. Output

## Questions

1. $f(x) = 3x^3 - 7x^2 + 3x - e^x + 2 = 0$

\begin{equation}
\begin{matrix}
g1(x) = {{e^x - 2}\over{3x^2 - 7x + 3}} = x \\
r1 = -0.24789639 \\
x_0 = -1 \\
number of steps = 17
\end{matrix}
\end{equation}

\begin{equation}
\begin{matrix}
g_2(x) = {{3x^4 - 7x^3 + 3x^2 + 2x} \over {e^x}} = x \\
r2 = 0.62616943 \\
x_0 = 1 \\
number of steps = 22
\end{matrix}
\end{equation}

\begin{equation}
\begin{matrix}
g_3(x) = \left({{-7x^2 + 3x - e^x + 2} \over -3.0}\right)^{1 \over 3} = x \\
r3 = 2.46222248 \\
x_0 = 2 \\
number of steps = 84
\end{matrix}
\end{equation}

\begin{equation}
\begin{matrix}
g_4(x) = \ln(3x^3 - 7x^2 + 3x + 2) = x \\
r4 = 6.07305409 \\
x_0 = 6 \\
number of steps = 34
\end{matrix}
\end{equation}

2. $S = |g'(r)|$

\begin{equation}
\begin{matrix}
{d \over dx}g_1(x) = {e^x(3x^2-13x+10)+2(6x-7) \over (3x^2-7x+3)^2} \\
|{d \over dx}g_1(r1)| = 0.269034
\end{matrix}
\end{equation}

\begin{equation}
\begin{matrix}
{d \over dx}g_2(x) = {-3x^4+19x^3-24x^2+4x+2 \over e^x} \\
|{d \over dx}g_2(r2)| = 0.375249
\end{matrix}
\end{equation}

\begin{equation}
\begin{matrix}
{d \over dx}g_3(x) = {1 \over 3}\left(-7x^2+3x-e^x+2 \over 3 \right )^{-2\over3} \left(-14x+3-e^x \over -3 \right ) \\
|{d \over dx}g_3(r3)| = 0.791783
\end{matrix}
\end{equation}

\begin{equation}
\begin{matrix}
{d \over dx}g_4(x) = {9x^2-14x+3 \over 3x^3-7x^2+3x+2} \\
|{d \over dx}g_4(r4)| = 0.575836
\end{matrix}
\end{equation}

3. $\lim_{k \to \infty} {e_{k+1} \over e_k} = S$
	1. For r1: $\lim_{k \to \infty} {e_{k+1} \over e_k} \approx 0.2690342181$
	2. For r2: $\lim_{k \to \infty} {e_{k+1} \over e_k} \approx 0.3752496237$
	3. For r3: $\lim_{k \to \infty} {e_{k+1} \over e_k} \approx 0.7917836336$
	4. For r4: $\lim_{k \to \infty} {e_{k+1} \over e_k} \approx 0.5758353631$

## Code

### Fixed Point Iterative Method

```Matlab
% Computes the fixed point of a function using the FPI.
% Written by Zachary Ferguson

function xc = fixed_point_iteration(g, x0, f, tol)
    % Compute the fixed point of g(x).
    % Input:
    %   g   - function to solve for the fixed point.
    %   x0  - initial guess
    %   f   - f(x) = g(x) - x
    %   tol - solution tolerance
    % Output:
    %   xc - computed root to the function g(x) = x.
    if nargin < 4
        tol = 1e-9;
    end

    r = fzero(f, x0);
    fprintf('r = %f\n', r);
    ei = 0;

    prev_x = x0;
    x = g(x0);
    n = 1;
    while (abs(prev_x - x) > 0.5 * tol)
        prev_x = x;
        x = g(x);
        n = n + 1;
        ei1 = abs(x - r);
        if (abs(prev_x - x) <= 0.5 * tol)
             fprintf('e_(i+1)/e_i = %.10f\n', ei1 / ei);
        end
        ei = ei1;
    end
    fprintf('n = %d\n', n);
    xc = x;
end
```

### Main

```Matlab
% MATH 446: Project 02
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 02\nWritten by Zachary Ferguson\n\n');

    fprintf('f(x) = 3*x^3 - 7*x^2 + 3*x - e^x + 2 = 0\n\n');
    f = @(x) 3*x^3 - 7*x^2 + 3*x - exp(x) + 2;

    fprintf('g1(x) = (e^x - 2) / (3x^2 - 7x + 3) = x\n')
    g1 = @(x) (exp(x) - 2) / (3*x^2 - 7*x + 3);
    fprintf('r1 = %.10f\n\n', fixed_point_iteration(g1, -1, f));

    fprintf('g2(x) = (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / e^x = x\n');
    g2 = @(x) (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / exp(x);
    fprintf('r2 = %.10f\n\n', fixed_point_iteration(g2, 1, f));

    fprintf('g3(x) = ((-7*x^2 + 3*x - e^x + 2) / -3.0)^(1/3) = x\n');
    g3 = @(x) ((-7*x^2 + 3*x - exp(x) + 2) / -3.0)^(1/3);
    fprintf('r3 = %0.10f\n\n', fixed_point_iteration(g3, 2, f));

    fprintf('g4(x) = ln(3*x^3 - 7*x^2 + 3*x + 2) = x\n');
    g4 = @(x) log(3*x^3 - 7*x^2 + 3*x + 2);
    fprintf('r4 = %.10f\n', fixed_point_iteration(g4, 6, f));
end
```

## Output

```
MATH 446: Project 02
Written by Zachary Ferguson

f(x) = 3*x^3 - 7*x^2 + 3*x - e^x + 2 = 0

g1(x) = (e^x - 2) / (3x^2 - 7x + 3) = x
r = -0.247896
e_(i+1)/e_i = 0.2690342181
n = 17
r1 = -0.2478963963

g2(x) = (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / e^x = x
r = 0.626169
e_(i+1)/e_i = 0.3752496237
n = 22
r2 = 0.6261694387

g3(x) = ((-7*x^2 + 3*x - e^x + 2) / -3.0)^(1/3) = x
r = 2.462222
e_(i+1)/e_i = 0.7917836336
n = 84
r3 = 2.4622224868

g4(x) = ln(3*x^3 - 7*x^2 + 3*x + 2) = x
r = 6.073054
e_(i+1)/e_i = 0.5758353631
n = 34
r4 = 6.0730540924
```
