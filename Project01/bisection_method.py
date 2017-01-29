"""
Computes the root of a function using the bisection method.
Written by Zachary Ferguson
"""


def compute_root(f, a, b, tol=1e-7):
    """
    Finds the root of the function, f, in the interval [a, b] within an
    absolute tolerance.
    """
    assert f(a) * f(b) <= 0
    # import pdb; pdb.set_trace()

    n = 0
    while((abs(b - a) / 2.0) > (0.5 * tol)):
        c = (a + b) / 2.0
        if(f(c) == 0):
            break

        if(f(a) * f(c) <= 0):
            b = c
        else:
            a = c
        n += 1
    print(n)
    return (a + b) / 2.0


def print_errors(f, r, y):
    """ Print the forward and backward error of r. """
    print("Forward error: %.10f" % abs(y - r))
    print("Backward error: %.10f" % abs(f(r)))

if __name__ == "__main__":
    ######
    # Q1a
    print("f1(x) = x^3 - 9:")
    f = lambda x: x**3 - 9
    r = compute_root(f, 2, 3)
    print("r = %.10f" % r)
    y = 3**(2. / 3.)
    print_errors(f, r, y)

    ######
    # Q1b
    print("\nf2(x) = 3x^3 + x^2 - x - 5:")
    f = lambda x: 3 * x**3 + x**2 - x - 5
    r = compute_root(f, 1, 2)
    print("r = %.10f" % r)
    y = (1. / 9. * (-1 + (593 - 27 * (481**0.5))**(1.0 / 3.0) + (593 + 27 *
        (481**0.5))**(1.0 / 3.0)))
    print_errors(f, r, y)

    ######
    # Q1c
    print("\nf3(x) = cos^2(x) - x + 6:")
    from math import cos
    f = lambda x: cos(x) * cos(x) - x + 6
    r = compute_root(f, 6, 7)
    print("r = %.10f" % r)
    y = 6.7760923163195023262
    print_errors(f, r, y)

    ######
    # Q3a
    print("\nf4(x) = 2x^3 - 6x - 1:")
    f = lambda x: 2.0 * x**3 - 6 * x - 1

    r = compute_root(f, -2, -1)
    print("r1 = %.10f" % r)
    y = -1.64178352745293
    print_errors(f, r, y)

    r = compute_root(f, -1, 0)
    print("r2 = %.10f" % r)
    y = -0.168254401781027
    print_errors(f, r, y)

    r = compute_root(f, 1, 2)
    print("r3 = %.10f" % r)
    y = 1.81003792923395
    print_errors(f, r, y)

    ######
    # Q3b
    print("\nf5(x) = e^(x-2) + x^3 - x:")
    from math import e
    f = lambda x: e**(x - 2) + x**3 - x

    r = compute_root(f, -2, -1)
    print("r1 = %.10f" % r)
    y = -1.0234821948582364944
    print_errors(f, r, y)

    r = compute_root(f, -0.5, 0.5)
    print("r2 = %.10f" % r)
    y = 0.16382224325010849634
    print_errors(f, r, y)

    r = compute_root(f, 0.5, 1.5)
    print("r3 = %.10f" % r)
    y = 0.78894138905554556637
    print_errors(f, r, y)

    ######
    # Q3c
    print("\nf6(x) = 1 + 5x - 6x^3 - e^(2x):")
    f = lambda x: 1 + 5 * x - 6 * x**3 - e**(2 * x)

    r = compute_root(f, -1.5, -0.5)
    print("r1 = %.10f" % r)
    y = -0.81809373448119542124
    print_errors(f, r, y)

    r = compute_root(f, -0.6, 0.4)
    print("r2 = %.10f" % r)
    y = 0.0
    print_errors(f, r, y)

    r = compute_root(f, 0.5, 1.5)
    print("r3 = %.10f" % r)
    y = 0.50630828634622119599
    print_errors(f, r, y)

    ######
    # Q4a
    print("\nf7(x) = x^2 - A:")
    print("A = 2, (a, b) = (1, 2), n = 22")
    A = 2
    f = lambda x: x**2 - A
    r = compute_root(f, 1, 2)
    print("r = %.10f" % r)
    y = 2**0.5
    print_errors(f, r, y)

    print("\nA = 3, (a, b) = (1, 2), n = 22")
    A = 3
    r = compute_root(f, 1, 2)
    print("r = %.10f" % r)
    y = 3**0.5
    print_errors(f, r, y)

    print("\nA = 5, (a, b) = (2, 3), n = 22")
    A = 5
    r = compute_root(f, 2, 3)
    print("r = %.10f" % r)
    y = 5**0.5
    print_errors(f, r, y)
