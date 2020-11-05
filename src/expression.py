"""
This module contains functions which compute for for a given expression M_{i+1} all possible M_i
which could be the predecessor of M_{i+1} before executing the loop body together with the associated
probabilities.
"""

from diofant import Expr, Symbol, simplify, Rational, symbols, Number
from mora.core import Program, RandomVar, Update

# Type aliases to improve readability
from src.utils import unique_symbol

Probability = Rational
Case = (Expr, Probability)


def get_cases_for_expression(expression: Expr, program: Program) -> [Case]:
    """
    The main function computing all possible expression_{i+1} together with the associated probabilities
    """
    result = [(expression, 1)]

    for symbol in reversed(program.variables):
        if hasattr(program.updates[symbol], "branches"):
            result = split_expressions_on_symbol(result, symbol, program)
        result = combine_expressions(result)

    return to_polynomials(result, program.variables)


def get_initial_value_for_expression(expression: Expr, program: Program) -> Number:
    """
    For a given expression returns its initial value
    """
    result = expression.xreplace({symbols("n", integer=True, positive=True): 0})
    for var, update in program.initial_values.items():
        if hasattr(update, 'branches') and len(update.branches) > 0:
            result = result.subs({var: update.branches[0][0]})
    return simplify(result)


def split_expressions_on_symbol(expressions: [Case], symbol: Symbol, program: Program):
    """
    Splits all given expressions on the possibilities of updating a given symbol
    """
    result = []
    for expr, prob in expressions:
        if symbol in program.updates.keys() and symbol in expr.free_symbols:
            for u, p in program.updates[symbol].branches:
                new_expression = simplify(expr.subs({symbol: u}))
                new_prob = prob * p
                result.append((new_expression, new_prob))
        else:
            result.append((expr, prob))

    return result


def split_expressions_on_rvs(expressions: [Case], program: Program):
    """
    Splits given expressions on all random variables
    """
    variables = program.variables.copy()
    for var in variables:
        if program.updates[var].is_random_var and not hasattr(program.updates[var], "branches"):
            expressions = split_expressions_on_rv(expressions, var, program)
    return expressions


def split_expressions_on_rv(expressions: [Case], rv: Symbol, program: Program):
    """
    Splits expressions on a given random variable rv. If there is a probability > 0 that the rv can be positive and
    the same for negative, then every expression containing rv splits into 3 expressions:
        1. rv is replaced by a random variable with only negative support
        2. rv is replaced by a random variable with support ranging over 0
        3. rv is replaced by a random variable with support only positive
    """
    low, high = program.updates[rv].random_var.get_support()
    if low > 0 or high < 0:
        return expressions

    split_rvs = []
    epsilon1 = unique_symbol("eps", real=True, positive=True)
    epsilon2 = unique_symbol("eps", real=True, positive=True)
    split_rvs.append(RandomVar("symbolic-support", (low, -epsilon1)))
    split_rvs.append(RandomVar("symbolic-support", (-epsilon1, epsilon2)))
    split_rvs.append(RandomVar("symbolic-support", (epsilon2, high)))

    cases = []
    for split_rv in split_rvs:
        var = unique_symbol("var")
        update = Update(var)
        update.is_random_var = True
        update.random_var = split_rv
        program.updates[var] = update
        program.variables.append(var)
        for expression, _ in expressions:
            cases.append((expression.xreplace({rv: var}), unique_symbol("p")))

    return cases



def combine_expressions(expressions: [Case]) -> [Case]:
    """
    In a given list of expressions with probabilities, combines equal expressions and their probabilities
    """
    tmp_map = {}
    for e, p in expressions:
        tmp_map[e] = tmp_map[e] + p if e in tmp_map else p
    return list(tmp_map.items())


def to_polynomials(expressions: [Case], variables) -> [Case]:
    """
    Converts all expressions in a list of cases to polynomials in the program variables
    """
    result = []
    for e, p in expressions:
        e = e.as_poly(variables)
        result.append((e, p))

    return result
