from typing import Iterable

from diofant import sympify, Rational, Poly, prod, Symbol, symbols, oo, Max, polylog, factorial, product, gamma
from diofant.stats import Laplace, E
from scipy.stats import norm
from math import sqrt
import re

LOG_NOTHING = 0
LOG_ESSENTIAL = 10
LOG_VERBOSE = 20
LOG_LEVEL = LOG_ESSENTIAL

class Update:
    # parse updates
    # takes string "x = P @ p; Q @ q" or x = RV(d, a, b)
    # creates class to deal with substituing powers of variables and moments
    def __init__(self, var, update_string=None, program_variables=None, is_random_var=False, random_var=None):
        self.is_random_var = is_random_var
        self.random_var = random_var
        self.var = var
        self.is_probabilistic = True

        if update_string is None:
            return

        # check if this is a RV or expression update
        rv = re.search(r"RV\((?P<params>.+)\)", update_string)
        if rv is not None:
            self.is_random_var = True
            dist, *params = map(str.strip, rv.group('params').split(','))
            params = list(map(sympify, params))
            params = [make_symbols_positive(p, program_variables) for p in params]
            self.random_var = RandomVar(dist, params, var_name=str(self.var))

        # here: if not is_random_var == else
        if not self.is_random_var:
            self.branches = []
            branches = update_string.split(";")
            for update in branches:
                if '@' in update:
                    exp, prob = update.split("@")
                else:
                    exp, prob = update, 1-sum([b[1] for b in self.branches])
                prob = sympify(prob)
                if not prob.is_zero:
                    exp = make_symbols_positive(sympify(exp), program_variables)
                    exp = make_floats_rational(exp)
                    prob = make_symbols_positive(prob, program_variables)
                    prob = make_floats_rational(prob)
                    self.branches.append((sympify(exp), prob))
            if sum([b[1] for b in self.branches]) != 1:
                raise Exception(f"Branch probabilities for {self.var} update do not sum up to 1. Terminating.")

    def update_term(self, term, pow):
        if self.is_random_var:
            return term.subs({self.var**pow: self.random_var.compute_moment(pow) })
        else:
            return term.subs({self.var**pow: self.power(pow)})

    def power(self, k):
        return sum(prob * (exp ** k) for exp, prob in self.branches)


class RandomVar:
    def __init__(self, distribution, parameters, var_name=None):
        self.distribution = distribution
        self.parameters = parameters
        self.var_name = var_name

    def get_support(self, k=1):
        if self.distribution == 'bernoulli':
            return sympify(0), sympify(1)

        if self.distribution == 'geometric':
            return sympify(1), oo

        if self.distribution == 'exponential':
            return sympify(0), oo

        if self.distribution == 'beta':
            return sympify(0), sympify(1)

        if self.distribution == 'uniform':
            l, u = self.parameters
            return interval_to_power(l, u, k)

        if self.distribution == 'chi-squared':
            l, u = self.parameters
            return sympify(0), oo

        if self.distribution == 'rayleigh':
            l, u = self.parameters
            return sympify(0), oo

        if self.distribution == 'symbolic-support':
            l, u = self.parameters
            return interval_to_power(l, u, k)

        if self.distribution == 'gauss':
            return interval_to_power(-oo, oo, k)

        if self.distribution == 'laplace':
            return interval_to_power(-oo, oo, k)

    def compute_moment(self, k):
        if self.distribution == 'finite':
            return sum([p * (b ** k) for b, p in self.parameters])

        if self.distribution == 'uniform':
            l, u = self.parameters
            return (u**(k+1)-l**(k+1))/((k+1)*(u-l))

        if self.distribution == 'gauss' or self.distribution == 'normal':
            mu, sigma_squared = self.parameters
            # For low moments avoid scipy.stats.moments as it does not support
            # parametric parameters. In the future get all moments directly,
            # using the following properties:
            # https://math.stackexchange.com/questions/1945448/methods-for-finding-raw-moments-of-the-normal-distribution
            if k == 0:
                return 1
            elif k == 1:
                return mu
            elif k == 2:
                return mu**2 + sigma_squared
            elif k == 3:
                return mu*(mu**2 + 3*sigma_squared)
            elif k == 4:
                return mu**4 + 6*mu**2*sigma_squared + 3*sigma_squared**2
            moment = norm(loc=mu, scale=sqrt(sigma_squared)).moment(k)
            return Rational(moment)

        if self.distribution == 'bernoulli':
            return sympify(self.parameters[0])

        if self.distribution == 'geometric':
            p = sympify(self.parameters[0])
            return p*polylog(-k, 1-p)

        if self.distribution == 'exponential':
            lambd = sympify(self.parameters[0])
            return factorial(k) / (lambd ** k)

        if self.distribution == 'beta':
            alpha, beta = self.parameters
            alpha = sympify(alpha)
            beta = sympify(beta)
            r = symbols('r')
            return product((alpha + r) / (alpha + beta + r), (r, 0, k-1))

        if self.distribution == 'chi-squared':
            n = sympify(self.parameters[0])
            i = symbols('i')
            return product(n + 2*i, (i, 0, k - 1))

        if self.distribution == 'rayleigh':
            s = sympify(self.parameters[0])
            return (2**(k / 2)) * (s**k) * gamma(1 + k/2)

        if self.distribution == 'unknown':
            return sympify(f"{self.var_name}(0)^{k}")

        if self.distribution == 'laplace':
            mu, b = self.parameters
            mu = sympify(mu)
            b = sympify(b)
            x = Laplace("x", mu, b)
            return E(x**k)


def EV(expression):
    if issubclass(type(expression), RandomVar):
        return expression.compute_moment(1)
    else:
        return expression


def get_exponent_of(var, mono):
    monoms = mono.as_poly([var]).monoms()
    if len(monoms) > 0 and len(monoms[0]) > 0:
        return monoms[0][0]
    return 0


def get_monoms(poly: Poly):
    """
    Returns the list of monoms for a given polynomial
    """
    monoms = []
    for powers in poly.monoms():
        m = prod(var ** power for var, power in zip(poly.gens, powers))
        if m != 1:
            monoms.append(m.as_poly(poly.gens))
    return monoms


def monomial_is_constant(monomial: Poly):
    """
    Returns true iff the given monomial is constant
    """
    if monomial.is_zero:
        return True
    powers = monomial.monoms()[0]
    return all(p == 0 for p in powers)


def is_independent_from_all(program, x, ys):
    """
    Returns true iff x is statistially independent from all ys, where x is from the current iteration and the ys
    are from the previous iteration.
    """
    if x not in program.ancestors[x]:
        return True

    for y in ys:
        if y in program.dependencies[x]:
            return False

    return True


def all_are_independent_from_all(program, xs, ys):
    """
    Returns true iff all xs are statistially independent from all ys, where the xs are from the current iteration
    and the ys are from the previous iteration.
    """
    for x in xs:
        if not is_independent_from_all(program, x, ys):
            return False
    return True


def get_powers_of_variable_in_polynomial(variable: Symbol, polynomial: Poly):
    """
    Returns the set of all powers p for which variable ** p occurs in the polynomial
    """
    monoms = get_monoms(polynomial)
    all_powers = []
    for monomial in monoms:
        powers = monomial.monoms()[0]
        powers_for_var = {var: power for var, power in zip(monomial.gens, powers) if power > 0}
        if variable in powers_for_var.keys():
            all_powers.append(powers_for_var[variable])
    all_powers.sort(reverse=True)
    return all_powers


def set_log_level(log_level):
    global LOG_LEVEL
    LOG_LEVEL = log_level


def log(message, level):
    """
    Logs a message depending on the log level
    """
    if level <= LOG_LEVEL:
        print(message)


def without_piecewise(expr):
    """
    Removes the Piecewise from an expression by assuming that all restricting assumptions are false.
    """
    if not expr.args:
        return expr

    if expr.is_Piecewise:
        return without_piecewise(expr.args[-1].expr)

    return expr.func(*[without_piecewise(a) for a in expr.args])


def make_symbols_positive(expr, exclude_symbols=None):
    """
    Ensures that all symbols in a given expression (with the exception of excluded symbols) are assumed to be positive
    """
    if expr.is_Symbol:
        if exclude_symbols is None or expr not in exclude_symbols:
            return symbols(expr.name, positive=True)

    if not expr.args:
        return expr

    return expr.func(*[make_symbols_positive(a, exclude_symbols) for a in expr.args])


def make_floats_rational(expr):
    """
    Converts all floats in an expression to rationals
    """
    if expr.is_Float:
        return Rational(expr)

    if not expr.args:
        return expr

    return expr.func(*[make_floats_rational(a) for a in expr.args])


def interval_to_power(low, high, power):
    """
    If x in [low, high] returns an interval [l, h] s.t. x**power in [l, h]
    """
    l = low ** power
    h = high ** power
    if power % 2 == 0:
        if high < 0:
            h, l = l, h
        elif low < 0:
            h = Max(h, l)
            l = sympify(0)
    return l, h
