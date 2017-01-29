"""
Computes the fixed point of a function using the FPI.
Written by Zachary Ferguson
"""


def compute_fixed_point(g, x0, tol=1e-8):
    """ Compute the fixed point of g(x). """
    prev_x = x0
    x = g(x0)
    while(abs(prev_x - x) > 0.5 * tol):
        prev_x = x
        x = g(x)
    return x

if __name__ == "__main__":
    print("g1(x) = (2x+2)^(1/3) = x:")
    print("xc = %.10f" % compute_fixed_point(lambda x: (2 * x + 2)**(1 / 3.),
        2))

    from math import log
    print("\ng2(x) = ln(7-x) = x:")
    print("xc = %.10f" % compute_fixed_point(lambda x: log(7 - x), 2))

    from math import sin
    print("\ng3(x) = ln(4-sin(x)) = x:")
    print("xc = %.10f" % compute_fixed_point(lambda x: log(4 - sin(x)), 2))

    A = 3.
    g = lambda x: (x + A / x) / 2.
    print("\ng4(x) = (x + 3 / x) / 2:")
    print("xc = %.10f" % compute_fixed_point(g, 2))

    A = 5.
    print("\ng5(x) = (x + 5 / x) / 2:")
    print("xc = %.10f" % compute_fixed_point(g, 2))

    from math import cos
    print("\ng6(x) = cos^2(x):")
    print("xc = %.10f" %
        compute_fixed_point(lambda x: (cos(x))**2, 1, tol=1e-6))
