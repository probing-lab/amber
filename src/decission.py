"""
This is the module which handles the decision on what proof-rule to apply to a program in order
to get something about its termination behavior. Then the proof-rule gets applied
"""

from mora.core import Program, get_solution as get_expected, get_recurrence, reset_mora
from mora.input import LOOP_GUARD_VAR
from diofant import Expr, sympify, symbols, expand

from . import branch_store, bound_store
from .initial_state_rule import InitialStateRule
from .supermartingale_rule import SupermartingaleRule
from .ranking_sm_rule import RankingSMRule
from .repulsing_sm_rule import RepulsingSMRule
from .rule import Result
from .utils import LOG_ESSENTIAL, log


def decide_termination(program: Program):
    """
    The main function, gathering all the information, deciding on and calling a proof-rule
    """
    reset_mora()
    branch_store.set_program(program)
    bound_store.set_program(program)
    lgc = get_loop_guard_change(program)
    me_pos = create_martingale_expression(program)
    me_neg = expand(me_pos * (-1))
    log(f"Martingale expression: {me_pos.as_expr()}", LOG_ESSENTIAL)
    rules = [
        InitialStateRule(lgc, me_pos, program),
        RankingSMRule(lgc, me_pos, program),
        SupermartingaleRule(lgc, me_pos, program),
        RepulsingSMRule(lgc, me_neg, program)
    ]
    result = Result()

    for rule in rules:
        if rule.is_applicable():
            result = rule.run(result)
            if result.all_known():
                break

    return result


def create_martingale_expression(program: Program):
    """
    Creates the martingale expression E(M_{i+1} - M_i | F_i). Also deterministic variables get substituted
    with their representation in n.
    """
    lg = symbols(LOOP_GUARD_VAR).as_poly(program.variables)
    expected_guard = get_recurrence(program, lg)
    lg = program.updates[symbols(LOOP_GUARD_VAR)].branches[0][0]
    expression = expand(expected_guard - lg).as_expr()
    expression = substitute_deterministic_variables(expression, program)
    return expand(expression)


def get_loop_guard_change(program: Program):
    """
    Returns E[LG_{n+1} - LG_{n}]
    """
    n = symbols("n", integer=True, positive=True)
    lg = sympify(LOOP_GUARD_VAR).as_poly(program.variables)
    expected_lg = get_expected(program, lg)
    expected_lg_plus = expected_lg.xreplace({n: n+1})
    return expand(expected_lg_plus - expected_lg)


def substitute_deterministic_variables(expr: Expr, program: Program):
    """
    Substitutes deterministic variables in a given expression with their representation in n.
    """
    for symbol, update in program.updates.items():
        if str(symbol) != LOOP_GUARD_VAR and not update.is_probabilistic:
            closed_form = get_expected(program, symbol.as_poly(program.variables))
            expr = expand(expr.xreplace({symbol: closed_form}))
    return expr
