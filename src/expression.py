"""
This module contains functions which compute for for a given expression M_{i+1} all possible M_i
which could be the predecessor of M_{i+1} before executing the loop body together with the associated
probabilities.
"""
from typing import Iterable, Tuple

from diofant import Expr, Symbol, simplify, Rational, symbols, Number, Min, Max
from mora.core import Program, RandomVar, Update

# Type aliases to improve readability
from src.utils import unique_symbol, get_monoms, flatten_substitution_choices

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


def get_initial_supports_for_variable_powers(var_powers: Iterable[Tuple[Symbol, Number]], program: Program):
    """
    Returns lower and upper bounds for the initial support of variable powers
    """
    supports = {}
    for v, p in var_powers:
        if program.initial_values[v].is_random_var:
            supports[v ** p] = program.initial_values[v].random_var.get_support(p)
        else:
            values = [b[0] ** p for b in program.initial_values[v].branches]
            supports[v ** p] = Min(*values), Max(*values)

    return supports


def get_initial_polarity_for_expression(expression: Expr, program: Program) -> (bool, bool):
    """
    Returns a sound estimate whether a given expression can initial be positive and negative. It does so
    by substituting the variable power in the given expression with all possible combination of lower and upper
    bounds of their respective supports.
    """
    variables = expression.free_symbols.intersection(program.variables)
    # First we can replace n and all variables which are deterministic initially (meaning they have exactly one branch)
    expression = expression.xreplace({symbols("n", integer=True, positive=True): 0})
    for v in variables:
        if not expression.free_symbols & variables:
            continue
        if hasattr(program.initial_values[v], "branches") and len(program.initial_values[v].branches) == 1:
            expression = expression.subs({v: program.initial_values[v].branches[0][0]})

    variables = expression.free_symbols.intersection(program.variables)
    if not variables:
        maybePos = bool(expression > 0) if (expression > 0).is_Boolean else True
        maybeNeg = bool(expression < 0) if (expression < 0).is_Boolean else True
        return maybePos, maybeNeg

    expression = expression.as_poly(variables)
    var_powers = set()
    monoms = get_monoms(expression)
    for m in monoms:
        m = m.as_poly(variables)
        powers = m.monoms()[0]
        var_powers.update([(v, p) for v, p in zip(m.gens, powers) if p > 0])
    initial_supports = get_initial_supports_for_variable_powers(var_powers, program)
    possible_substitutions = flatten_substitution_choices(initial_supports)
    expression = expression.as_expr()
    possible_initial_polarities = []
    for ps in possible_substitutions:
        pos = expression.subs(ps) > 0
        if pos.is_Boolean:
            possible_initial_polarities.append(bool(pos))
        else:
            possible_initial_polarities.append(None)

    allPositive = all([v is True for v in possible_initial_polarities])
    allNegative = all([v is False for v in possible_initial_polarities])

    maybePos = not allNegative
    maybeNeg = not allPositive

    return maybePos, maybeNeg


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
    epsilon = unique_symbol("eps", real=True, positive=True)
    split_rvs.append(RandomVar("symbolic-support", (low, -epsilon)))
    split_rvs.append(RandomVar("symbolic-support", (-epsilon, epsilon)))
    split_rvs.append(RandomVar("symbolic-support", (epsilon, high)))

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
