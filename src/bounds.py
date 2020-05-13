from mora.mora import mora
from . import branch_store, bound_store


def bounds(benchmark, expression):
    program = mora(benchmark, goal=1)
    branch_store.set_program(program)
    bound_store.set_program(program)
    bounds = bound_store.get_bounds_of_expr(expression)
    print("Expression: ", bounds.expression)
    print("Lower bound: ", bounds.lower)
    print("Upper bound: ", bounds.upper)
    print("Absolute upper bound: ", bounds.absolute_upper)
    print("Maybe positive: ", bounds.maybe_positive)
    print("Maybe negative: ", bounds.maybe_negative)
