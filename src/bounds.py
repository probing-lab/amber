from mora.mora import mora
from . import branch_store, bound_store
from .utils import log, LOG_ESSENTIAL


def bounds(benchmark, expression):
    program = mora(benchmark, goal=1)
    branch_store.set_program(program)
    bound_store.set_program(program)
    bounds = bound_store.get_bounds_of_expr(expression)
    log(f"Expression: {bounds.expression}", LOG_ESSENTIAL)
    log(f"Lower bound: {bounds.lower}", LOG_ESSENTIAL)
    log(f"Upper bound: {bounds.upper}", LOG_ESSENTIAL)
    log(f"Absolute upper bound: {bounds.absolute_upper}", LOG_ESSENTIAL)
    log(f"Maybe positive: {bounds.maybe_positive}", LOG_ESSENTIAL)
    log(f"Maybe negative: {bounds.maybe_negative}", LOG_ESSENTIAL)
