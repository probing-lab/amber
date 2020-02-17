"""
This module contains functions deciding whether or not a given expression is an invariant of the program,
More precisely, it decides whether expression <= 0 is eventually invariant.
The methods are of course not complete in general.
"""
import math

from diofant import Expr, sympify, symbols, solve, Symbol
from mora.core import Program


def is_invariant(expression: Expr, program: Program):
    """
    Main function deciding whether expression <= 0 is eventually invariant
    """
    expression = strap_expression(expression)
    n = symbols('n')
    is_deterministic = len(expression.free_symbols.difference({n})) == 0
    if is_deterministic:
        return is_deterministic_invariant(expression)
    else:
        return is_probabilistic_invariant(expression, program)


def is_deterministic_invariant(expression: Expr):
    """
    Checks whether an expression only containing n eventually becomes <= 0
    """
    n = symbols('n')
    max_0 = get_max_0(expression, n)
    return expression.subs({n: max_0 + 1}) <= 0


def is_probabilistic_invariant(expression: Expr, program: Program):
    """
    For more complex expressions the decision whether expression <= 0 is invariant gets dispatched.
    """
    # TODO: dispatch to Z3
    raise NotImplementedError()


def get_max_0(expression: Expr, n: Symbol):
    """
    Returns the maximum 0 of a given expression or 0.
    """
    try:
        zeros = solve(expression, n)
    except NotImplementedError:
        zeros = []
    zeros = [math.ceil(float(z[n])) for z in zeros] + [0]
    return max(zeros)


def strap_expression(expression: Expr):
    expression = expression.args[0] if len(expression.args) > 0 else expression
    return sympify(str(expression))