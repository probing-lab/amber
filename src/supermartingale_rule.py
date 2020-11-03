"""
This module implements the general proof rule for AST
"""

from diofant import symbols, sympify

from . import bound_store
from .asymptotics import is_dominating_or_same, Direction, Answer
from .expression import get_cases_for_expression
from .invariance import is_invariant
from .rule import Rule, Result, Witness
from .utils import amber_limit


class SupermartingaleRule(Rule):
    def is_applicable(self):
        n = symbols("n", integer=True, positive=True)
        lim = amber_limit(self.loop_guard_change, n)
        return lim <= 0

    def run(self, result: Result):
        if result.AST.is_known():
            return result

        # Martingale expression has to be <= 0 eventually
        if not is_invariant(self.martingale_expression, self.program):
            return result

        # Eventually one branch of LG_{i+1} - LG_i has to decrease more or equal than constant
        branches = get_cases_for_expression(sympify(self.program.loop_guard), self.program)
        for branch, prob in branches:
            bounds = bound_store.get_bounds_of_expr(branch - sympify(self.program.loop_guard))
            n = symbols("n", integer=True, positive=True)
            if is_dominating_or_same(bounds.upper, sympify(-1), n, direction=Direction.NegInf):
                result.AST = Answer.TRUE
                result.add_witness(ASTWitness(
                    self.program.loop_guard,
                    self.martingale_expression,
                    branch,
                    bounds.upper,
                    prob
                ))
                return result

        return result


class ASTWitness(Witness):

    def __init__(self, martingale, martingale_expression, decreasing_branch, bound, prob):
        super(ASTWitness, self).__init__("AST")
        martingale = sympify(martingale).as_expr()
        martingale_expression = sympify(martingale_expression).as_expr()
        decreasing_branch = sympify(decreasing_branch).as_expr()
        bound = sympify(bound).as_expr()
        self.data = {
            "SM": martingale,
            "SM expression": martingale_expression,
            "Decreasing branch": decreasing_branch,
            "Branch change bound": bound,
            "Probability": prob
        }
        self.explanation = f"Eventually, '{martingale}' is a supermartingale. Also eventually, taking the branch\n" \
                           f"'{decreasing_branch}' (which happens with probability {prob}) " \
                           f"changes the supermartingale by at least {bound}."
