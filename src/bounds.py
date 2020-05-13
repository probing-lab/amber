from mora.input import InputParser
from . import branch_store, bound_store


def bounds(benchmark, expression):
    ip = InputParser()
    ip.set_source(benchmark)
    program = ip.parse_source()
    branch_store.set_program(program)
    bound_store.set_program(program)
    bounds = bound_store.get_bounds_of_expr(expression)
    print("Expression: ", bounds.expression)
    print("Lower bound: ", bounds.lower)
    print("Upper bound: ", bounds.upper)
    print("Absolute upper bound: ", bounds.absolute_upper)
    print("Maybe positive: ", bounds.maybe_positive)
    print("Maybe negative: ", bounds.maybe_negative)
