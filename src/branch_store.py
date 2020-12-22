"""
This modules contains functions providing the branches of a given monomial.

Example:
while true:
    x = x - 1 @ 1/2; x + 1
    y = y - 1 @ 1/2; y + 1

Branches of x*y:
- (x - 1)(y - 1) @ 1/4
- (x + 1)(y - 1) @ 1/4
- (x - 1)(y + 1) @ 1/4
- (x + 1)(y + 1) @ 1/4

The branches of monomials are computed just in time and stored so they can be reused.
"""

from diofant import *
from mora.core import Program
from .expression import get_cases_for_expression, get_initial_polarity_for_expression


class Branch:
    """
    The class actually representing a branches of a monomial. It represents:
    monom_n = recurrence_constant * monom_{n - 1} + inhom_part_{n - 1} @ probability
    """
    monom: Expr
    recurrence_constant: Number
    inhom_part: Poly
    probability: Number
    initial_value: Number


store = {}
initial_value_store = {}
program: Program = None


def set_program(p: Program):
    """
    Set the program and initialize the store. This function needs to be called before the store is used.
    """
    global program, store, initial_value_store
    program = p
    store = {}
    initial_value_store = {}


def get_branches_of_monom(monom: Expr) -> [Branch]:
    """
    Lazily computes the branches of a given monomial and returns them.
    """
    monom = sympify(monom)
    if monom not in store:
        __compute_branches(monom)
    return store[monom]


def get_initial_polarity_of_monom(monom: Expr) -> (bool, bool):
    """
    Lazily computes the initial value of a given monomial and returns them.
    """
    global program, initial_value_store
    monom = sympify(monom)
    if monom not in initial_value_store:
        initial_value_store[monom] = get_initial_polarity_for_expression(monom, program)
    return initial_value_store[monom]


def __compute_branches(monom: Expr):
    global program
    monom = sympify(monom)
    cases = get_cases_for_expression(monom, program)
    branches = __cases_to_branches(cases, monom)
    store[monom] = branches


def __cases_to_branches(cases, monom):
    result = []
    for case in cases:
        branch = Branch()
        branch.monom = monom
        branch.recurrence_constant = case[0].coeff_monomial(monom)
        branch.inhom_part = case[0] - (branch.recurrence_constant * monom)
        branch.probability = case[1]
        result.append(branch)
    return result
