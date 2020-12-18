"""
This module implements the simple rule that checks whether the loop terminates immediately
because of the initial condition
"""
from diofant import sympify

from .expression import get_initial_polarity_for_expression
from .rule import Rule, Result
from .utils import Answer


class InitialStateRule(Rule):
    def is_applicable(self):
        return True

    def run(self, result: Result):
        loop_guard = sympify(self.program.loop_guard)
        maybePos, _ = get_initial_polarity_for_expression(loop_guard, self.program)

        if maybePos:
            return result

        result.PAST = Answer.TRUE
        result.AST = Answer.TRUE
        result.NONTERM = Answer.FALSE

        return result
