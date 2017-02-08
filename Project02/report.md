# MATH 446: Project 02
#### Zachary Ferguson
---

## Contents

1. Questions
2. Code
3. Output

## Questions

1. $f(x) = 3x^3 - 7x^2 + 3x - e^x + 2 = 0$
    1. $g1(x) = {{e^x - 2}\over{3x^2 - 7x + 3}} = x$ <br>
        $r1 = -0.24789639$<br>
        $x_0 = -1$<br>
        $number of steps = 17$
    2. $g_2(x) = {{3x^4 - 7x^3 + 3x^2 + 2x} \over {e^x}} = x$ <br>
        $r2 = 0.62616943$<br>
        $x_0 = 1$<br>
        $number of steps = 22$
    3. $g_3(x) = \left({{-7x^2 + 3x - e^x + 2} \over -3.0}\right)^{1 \over 3} = x$ <br>
        $r3 = 2.46222248$<br>
        $x_0 = 2$<br>
        $number of steps = 84$
    4. $g_4(x) = \ln(3x^3 - 7x^2 + 3x + 2) = x$ <br>
        $r4 = 6.07305409$<br>
        $x_0 = 6$<br>
        $number of steps = 34$
2. $S = |g'(r)|$
    1. ${d \over dx}g_1(x) = {e^x(3x^2-13x+10)+2(6x-7) \over (3x^2-7x+3)^2}$<br>
$|{d \over dx}g_1(r1)| = 0.269034$
    2. ${d \over dx}g_2(x) = {-3x^4+19x^3-24x^2+4x+2 \over e^x}$<br>
$|{d \over dx}g_2(r2)| = 0.375249$
    3. ${d \over dx}g_3(x) = {1 \over 3}\left(-7x^2+3x-e^x+2 \over 3 \right )^{-2\over3} \left(-14x+3-e^x \over -3 \right )$
$|{d \over dx}g_3(r3)| = 0.791783$
    4. ${d \over dx}g_4(x) = {9x^2-14x+3 \over 3x^3-7x^2+3x+2}$
$|{d \over dx}g_4(r4)| = 0.575836$

## Code

```Matlab
% Computes the fixed point of a function using the FPI.
% Written by Zachary Ferguson


function fixed_point_iteration
    fprintf('Fixed Point Iteration\nWrtten by Zachary Ferguson\n\n');

    fprintf('f(x) = 3*x^3 - 7*x^2 + 3*x - e^x + 2 = 0\n\n');

    fprintf('g1(x) = (e^x - 2) / (3x^2 - 7x + 3) = x\n')
    g1 = @(x) (exp(x) - 2) / (3*x^2 - 7*x + 3);
    fprintf('r1 = %.10f\n\n', compute_fixed_point(g1, -1));

    fprintf('g2(x) = (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / e^x = x\n');
    g2 = @(x) (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / exp(x);
    fprintf('r2 = %.10f\n\n', compute_fixed_point(g2, 1));

    fprintf('g3(x) = ((-7*x^2 + 3*x - e^x + 2) / -3.0)^(1/3) = x\n');
    g3 = @(x) ((-7*x^2 + 3*x - exp(x) + 2) / -3.0)^(1/3);
    fprintf('r3 = %0.10f\n\n', compute_fixed_point(g3, 2));

    fprintf('g4(x) = ln(3*x^3 - 7*x^2 + 3*x + 2) = x\n');
    g4 = @(x) log(3*x^3 - 7*x^2 + 3*x + 2);
    fprintf('r4 = %.10f\n', compute_fixed_point(g4, 6));
end


% Compute the fixed point of g(x).
function xc = compute_fixed_point(g, x0, tol)
    if nargin < 3
        tol = 1e-9;
    end

    prev_x = x0;
    x = g(x0);
    n = 1;
    while (abs(prev_x - x) > 0.5 * tol)
        prev_x = x;
        x = g(x);
        n = n + 1;
    end
    fprintf('n = %d\n', n);
    xc = x;
end
```

## Output

```
Fixed Point Iteration
Wrtten by Zachary Ferguson

f(x) = 3*x^3 - 7*x^2 + 3*x - e^x + 2 = 0

g1(x) = (e^x - 2) / (3x^2 - 7x + 3) = x
n = 17
r1 = -0.2478963963

g2(x) = (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / e^x = x
n = 22
r2 = 0.6261694387

g3(x) = ((-7*x^2 + 3*x - e^x + 2) / -3.0)^(1/3) = x
n = 84
r3 = 2.4622224868

g4(x) = ln(3*x^3 - 7*x^2 + 3*x + 2) = x
n = 34
r4 = 6.0730540924
```
