"""
Computes the fixed point of a function using the FPI.
Written by Zachary Ferguson
"""


def compute_fixed_point(g, x0, tol=1e-8, print_n=False):
    """ Compute the fixed point of g(x). """
    prev_x = x0
    x = g(x0)
    n = 1
    while(abs(prev_x - x) > 0.5 * tol):
        prev_x = x
        x = g(x)
        n += 1
    if(print_n):
        print("\tn = %d" % n)
    return x

if __name__ == "__main__":
    print("Fixed Point Iteration\nWrtten by Zachary Ferguson\n")

    print("Q1a:\n\tg(x) = (2x+2)^(1/3) = x")
    print("\txc = %.10f" % compute_fixed_point(lambda x: (2 * x + 2)**(1 / 3.),
        2))

    from math import log
    print("Q1b:\n\tg(x) = ln(7-x) = x")
    print("\txc = %.10f" % compute_fixed_point(lambda x: log(7 - x), 2))

    from math import sin
    print("Q1c:\n\tg(x) = ln(4-sin(x)) = x")
    print("\txc = %.10f" % compute_fixed_point(lambda x: log(4 - sin(x)), 2))

    A = 3.
    g = lambda x: (x + A / x) / 2.
    print("Q3a:\n\tg(x) = (x + 3 / x) / 2")
    print("\tx0 = 2")
    print("\txc = %.10f" % compute_fixed_point(g, 2, print_n=True))

    A = 5.
    print("Q3b:\n\tg(x) = (x + 5 / x) / 2")
    print("\tx0 = 2")
    print("\txc = %.10f" % compute_fixed_point(g, 2, print_n=True))

    from math import cos
    print("Q5:\n\tg(x) = cos^2(x)")
    xc = compute_fixed_point(lambda x: (cos(x))**2, 1, tol=1e-6, print_n=True)
    print("\txc = %.10f" % xc)
    print("\n\tg'(x) = -2*cos(x)*sin(x)")
    print("\t|g'(xc)| = %.10f" % abs(-2 * cos(xc) * sin(xc)))
    print("\tTherefore g(x) is locally convergent to xc.")
