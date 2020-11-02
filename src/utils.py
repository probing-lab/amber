import math
from enum import Enum, auto
from diofant import *

from mora.core import Program


LOG_NOTHING = 0
LOG_ESSENTIAL = 10
LOG_VERBOSE = 20
LOG_LEVEL = LOG_ESSENTIAL


class Answer(Enum):
    FALSE = auto()
    TRUE = auto()
    UNKNOWN = auto()

    def is_true(self):
        return self is Answer.TRUE

    def is_known(self):
        return self is not Answer.UNKNOWN

    def __str__(self):
        if self is Answer.TRUE:
            return "Yes"
        if self is Answer.FALSE:
            return "No"
        if self is Answer.UNKNOWN:
            return "Maybe"


__COUNTER = 0


def unique_symbol(s: str, **args):
    """
    Returns a symbol which every time has a different name
    """
    global __COUNTER
    s = symbols(s + str(__COUNTER), **args)
    __COUNTER += 1
    return s


def get_max_0(expression: Expr, n: Symbol):
    """
    Returns the maximum positive 0 of a given expression or 0 if it does not exist
    """
    try:
        exp_zeros = solve(expression, n)
        if exp_zeros == [{}]:
            return 0
        exp_zeros = [z[n] for z in exp_zeros if z[n].is_real]
    except NotImplementedError:
        exp_zeros = []
    exp_zeros = [math.ceil(float(z)) for z in exp_zeros] + [0]
    return max(exp_zeros)


def get_monoms(poly: Poly):
    """
    Returns for the list of monoms for a given polynomial
    """
    monoms = []
    for powers in poly.monoms():
        m = prod(var ** power for var, power in zip(poly.gens, powers))
        if m != 1:
            monoms.append(m)
    return monoms


def get_polarity(expression: Expr, n: Symbol):
    """
    Given an expression in n, returns whether or not the expression is positive and negative for some values of n
    """
    expression = simplify(expression)
    if expression.is_number:
        return expression > 0, expression < 0

    max_0 = get_max_0(expression, n)
    if max_0 > 0:
        pos = True
        neg = True
    else:
        expr_1 = expression.subs({n: 1})
        pos = expr_1 > 0
        neg = expr_1 < 0

    return pos, neg


def get_signums_in_expression(expression: Expr) -> [Expr]:
    """
    Given an expression it returns all expressions which occur within the signum function.
    E.g for sign(x+2)*sign(y**3) it returns [x+1, y**3]
    """
    if isinstance(expression, sign):
        return [expression.args[0]]

    signums = []
    for arg in expression.args:
        signums += get_signums_in_expression(arg)

    return signums


def get_all_monom_powers(monom: Expr) -> [Number]:
    """
    Returns the degrees of all variables in a given monomial in a list
    """
    monom = monom.as_poly(monom.free_symbols)
    return list(monom.degree_list())


def monom_is_deterministic(monom: Expr, program: Program):
    """
    Returns true iff a given monomial is deterministic, that means all variables in the monomial are deterministic
    """
    variables_deterministic = [not program.updates[m].is_probabilistic for m in monom.free_symbols]
    return all(variables_deterministic)


def set_log_level(log_level):
    global LOG_LEVEL
    LOG_LEVEL = log_level


def log(message, level):
    """
    Logs a message depending on the log level
    """
    if level <= LOG_LEVEL:
        print(message)
