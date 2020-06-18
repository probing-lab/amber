"""
This module implements the geometric proof rule
"""

from diofant import sympify, symbols
from . import bound_store
from .asymptotics import is_dominating_or_same, Direction
from .expression import get_cases_for_expression
from .rule import Rule, Result, Witness
from .utils import Answer


class GeometricRule(Rule):
    def is_applicable(self):
        return True

    def run(self, result: Result):
        if result.PAST.is_known():
            return result

        branches = get_cases_for_expression(sympify(self.program.loop_guard), self.program)
        branches = [branch for branch, _ in branches]
        bounds = [bound_store.get_bounds_of_expr(branch) for branch in branches]
        n = symbols('n')

        for bound in bounds:
            if is_dominating_or_same(bound.upper, sympify(-1), n, direction=Direction.NegInf):
                result.PAST = Answer.TRUE
                result.AST = Answer.TRUE
                result.add_witness(GeometricWitness(bound.expression, bound.upper))
                return result

        return result


class GeometricWitness(Witness):

    def __init__(self, branch, bound):
        super(GeometricWitness, self).__init__("PAST")
        self.data = {
            "Loop Guard Branch": branch,
            "Bound": bound,
        }
        self.explanation = f"The loop guard has the branch '{branch}' which is eventually upper bounded by '{bound}'.\n"\
                           f"Therefore, the termination time follows a geometric distribution after the point from\n"\
                           f"which the bound holds."
