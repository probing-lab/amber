"""
This module implements the simple rule that checks whether the loop terminates immediately
because of the initial condition
"""
from diofant import sympify

from .rule import Rule, Result, Witness
from .utils import Answer


class InitialStateRule(Rule):
    def is_applicable(self):
        return True

    def run(self, result: Result):
        loop_guard = sympify(self.program.loop_guard)
        for var, update in self.program.initial_values.items():
            if hasattr(update, "branches"):
                loop_guard = loop_guard.subs({var: update.branches[0][0]})

        if not loop_guard.is_number or bool(loop_guard > 0):
            return result

        result.PAST = Answer.TRUE
        result.AST = Answer.TRUE
        result.add_witness(InitialStateWitness(loop_guard))

        return result


class InitialStateWitness(Witness):

    def __init__(self, loop_guard):
        super(InitialStateWitness, self).__init__("PAST")
        self.data = {
            "Loop guard in initial state": loop_guard,
        }
        self.explanation = f"The value of the loop guard in the initial state is '{loop_guard}'. " \
                           f"Therefore the loop is not entered."
