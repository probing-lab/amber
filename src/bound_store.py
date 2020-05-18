"""
This modules contains functions providing the bounds of given monomials and polynomial expressions.
"""

from diofant import *
from mora.core import Program
from .utils import *
from .asymptotics import *
from . import branch_store

store = {}
program = None


class Bounds:
    expression: Poly
    lower: Expr
    upper: Expr
    maybe_positive: bool
    maybe_negative: bool
    __absolute_upper__: Expr = None

    @property
    def absolute_upper(self):
        if self.__absolute_upper__ is None:
            n = symbols('n')
            self.__absolute_upper__ = get_eventual_bound([self.upper, self.lower * -1], n)
        return self.__absolute_upper__


def set_program(p: Program):
    """
    Set the program and initialize the store. This function needs to be called before the store is used.
    """
    global program, store
    program = p
    store = {}


def get_bounds_of_expr(expression: Expr) -> Bounds:
    """
    Computes the bounds of a polynomial over the program variables. It does so by substituting the bounds of the monomials.
    """
    expression = sympify(expression)
    variables = set(program.variables).difference({symbols('n')})
    expression = expression.as_poly(variables)
    expr_bounds = __initialize_bounds_for_expression(expression)
    monoms = get_monoms(expression)
    for monom in monoms:
        monom_bounds = __get_bounds_of_monom(monom)
        __replace_monom_in_expr_bounds(monom, monom_bounds, expression, expr_bounds)

    upper_candidates = __split_on_signums(expr_bounds.upper.as_expr())
    lower_candidates = __split_on_signums(expr_bounds.lower.as_expr())

    expr_bounds.upper = get_eventual_bound(upper_candidates, symbols('n'), direction=Direction.PosInf)
    expr_bounds.lower = get_eventual_bound(lower_candidates, symbols('n'), direction=Direction.NegInf)
    return expr_bounds


def __replace_monom_in_expr_bounds(monom, monom_bounds: Bounds, expression: Poly, expr_bounds: Bounds):
    """
    Helper function which replaces a single monomial with its bounds. Which bound to take depends on the coefficient
    of the monomial.
    """
    coeff = expression.coeff_monomial(monom)
    if coeff > 0:
        upper = monom_bounds.upper
        lower = monom_bounds.lower
        pos = monom_bounds.maybe_positive
        neg = monom_bounds.maybe_negative
    else:
        upper = monom_bounds.lower
        lower = monom_bounds.upper
        pos = monom_bounds.maybe_negative
        neg = monom_bounds.maybe_positive

    expr_bounds.upper = expr_bounds.upper.subs({monom: upper})
    expr_bounds.lower = expr_bounds.lower.subs({monom: lower})
    # Rough estimate of whether the expression is positive/negative
    expr_bounds.maybe_positive = expr_bounds.maybe_positive or pos
    expr_bounds.maybe_negative = expr_bounds.maybe_negative or neg


def __initialize_bounds_for_expression(expression: Poly) -> Bounds:
    """
    Initializes the bounds object for an expression by setting the lower and upper bounds equal to the expression.
    """
    bounds = Bounds()
    bounds.expression = expression.as_expr()
    bounds.lower = expression.copy()
    bounds.upper = expression.copy()

    # Initialize the polarity of the expression by the polarity of the deterministic part
    n_expr = expression.coeff_monomial(1)
    pos, neg = get_polarity(n_expr, symbols('n'))

    bounds.maybe_positive = pos
    bounds.maybe_negative = neg
    return bounds


def __get_bounds_of_monom(monom: Expr) -> Bounds:
    """
    Computes the bounds of a monomial in a lazy way
    """
    monom = sympify(monom).as_expr()
    if monom not in store:
        __compute_bounds_of_monom(monom)
    return store[monom]


def __compute_bounds_of_monom(monom: Expr):
    global program
    print(f"Computing bounds for {monom.as_expr()}")
    if monom_is_deterministic(monom, program):
        __compute_bounds_of_deterministic_monom(monom)
    elif len(monom.free_symbols) == 1 and get_all_monom_powers(monom)[0] > 2:
        variable = monom.free_symbols.pop()
        power = get_all_monom_powers(monom)[0]
        if power % 2 == 0:
            monom = variable ** 2
            power = int(power / 2)
        else:
            monom = variable
            power = power
        __compute_bounds_of_monom_power(monom, power)
    else:
        __compute_bounds_of_monom_recurrence(monom)
    print(f"Found bounds for {monom.as_expr()}")


def __compute_bounds_of_deterministic_monom(monom):
    """
    Computes the bounds of a deterministic monomial by replacing its variables by their first moments, which are
    their exact closed-form representations
    """
    global program
    bound = monom
    for variable in monom.free_symbols:
        moment = program.moments[str(variable) + "^1"]
        bound = bound.subs({variable: moment})

    n_int = symbols('n', integer=True)
    n = symbols('n')
    bound = bound.subs({n_int: n})
    pos, neg = get_polarity(bound, n)
    bound = simplify_asymptotically(bound, n)

    bounds = Bounds()
    bounds.expression = monom
    bounds.upper = bound
    bounds.lower = bound
    bounds.maybe_positive = pos
    bounds.maybe_negative = neg

    store[bounds.expression] = bounds


def __compute_bounds_of_monom_power(monom: Expr, power: Number):
    """
    Computes the bounds of monom**odd_power by just taking the bounds of monom and raising it to the given power.
    This is only sound if the given power is odd or the monom is always positive
    """
    n = symbols('n')
    monom_bounds = __get_bounds_of_monom(monom)
    upper_bound = simplify_asymptotically(monom_bounds.upper ** power, n)
    lower_bound = simplify_asymptotically(monom_bounds.lower ** power, n)

    bounds = Bounds()
    bounds.expression = (monom ** power).as_expr()
    bounds.upper = upper_bound
    bounds.lower = lower_bound
    bounds.maybe_positive = monom_bounds.maybe_positive
    bounds.maybe_negative = monom_bounds.maybe_negative

    store[bounds.expression] = bounds


def __compute_bounds_of_monom_recurrence(monom: Expr):
    """
    Computes the bounds of a monomial by representing it as a recurrence relation
    """
    n = symbols('n')
    branches = branch_store.get_branches_of_monom(monom)
    inhom_parts_bounds = [get_bounds_of_expr(b.inhom_part) for b in branches]
    initial_value = branch_store.get_initial_value_of_monom(monom)
    maybe_pos, maybe_neg = __get_monom_polarity(monom, inhom_parts_bounds, initial_value)

    inhom_parts_bounds_lower = [b.lower.subs({n: n - 1}) for b in inhom_parts_bounds]
    inhom_parts_bounds_upper = [b.upper.subs({n: n - 1}) for b in inhom_parts_bounds]

    max_upper = get_eventual_bound(inhom_parts_bounds_upper, n, direction=Direction.PosInf)
    min_lower = get_eventual_bound(inhom_parts_bounds_lower, n, direction=Direction.NegInf)
    min_rec = min([b.recurrence_constant for b in branches])
    max_rec = max([b.recurrence_constant for b in branches])
    starting_values = __get_starting_values(maybe_pos, maybe_neg)

    upper_candidates = __compute_bound_candidates({min_rec, max_rec}, {max_upper}, starting_values)
    lower_candidates = __compute_bound_candidates({min_rec, max_rec}, {min_lower}, starting_values)

    max_upper_candidate = get_eventual_bound(upper_candidates, n, direction=Direction.PosInf)
    min_lower_candidate = get_eventual_bound(lower_candidates, n, direction=Direction.NegInf)

    # If monom is negative upper bound cannot be larger than 0
    if not maybe_pos:
        max_upper_candidate = get_eventual_bound([max_upper_candidate, sympify(0)], n, direction=Direction.NegInf)
    # If monom is positive lower bound cannot be smaller than 0
    if not maybe_neg:
        min_lower_candidate = get_eventual_bound([min_lower_candidate, sympify(0)], n, direction=Direction.PosInf)

    bounds = Bounds()
    bounds.expression = monom.as_expr()
    bounds.upper = max_upper_candidate.as_expr()
    bounds.lower = min_lower_candidate.as_expr()
    bounds.maybe_positive = maybe_pos
    bounds.maybe_negative = maybe_neg

    store[bounds.expression] = bounds


def __get_monom_polarity(monom: Expr, inhom_parts_bounds: [Bounds], initial_value: Number) -> (bool, bool):
    """
    Returns a rough but sound estimate of whether or not a given monomial can be positive and negative
    """
    # If all powers of the variables are even, the monomial is only positive
    powers = get_all_monom_powers(monom)
    all_powers_even = all([p % 2 == 0 for p in powers])
    if all_powers_even:
        return True, False

    # Otherwise estimate monom polarity by polarity of initial condition and polarity of inhomogenous parts
    initial_pos = sympify(initial_value) > 0
    if initial_pos.is_Relational:
        initial_pos = True
    else:
        initial_pos = bool(initial_pos)

    initial_neg = sympify(initial_value) < 0
    if initial_neg.is_Relational:
        initial_neg = True
    else:
        initial_neg = bool(initial_neg)

    maybe_pos = initial_pos or any([b.maybe_positive for b in inhom_parts_bounds])
    maybe_neg = initial_neg or any([b.maybe_negative for b in inhom_parts_bounds])
    return maybe_pos, maybe_neg


def __get_starting_values(maybe_pos: bool, maybe_neg: bool) -> [Expr]:
    """
    Returns the possible values of a monomial after which it is within a certain bound. This gets used
    to solve the recurrence relations for the bounds candidates.
    """
    values = []
    if maybe_pos:
        values.append(unique_symbol('d', positive=True, real=True))
    if maybe_neg:
        values.append(unique_symbol('d', positive=True, real=True) * -1)
    if not maybe_pos and not maybe_neg:
        values.append(sympify(0))

    return values


def __compute_bound_candidates(coefficients: [Number], inhom_parts: [Expr], starting_values: [Expr]) -> [Expr]:
    """
    Computes functions which could potentially be bounds
    """
    candidates = []

    for c in coefficients:
        for part in inhom_parts:
            c0 = symbols('c0')
            solution = __compute_bound_candidate(c, part, c0)
            for v in starting_values:
                solution = solution.subs({c0: v})
                # If a candidate contains signum functions, we have to split the candidate into more candidates
                new_candidates = __split_on_signums(solution)
                candidates += new_candidates

    return candidates


def __compute_bound_candidate(c: Number, inhom_part: Expr, starting_value: Expr) -> Expr:
    """
    Computes a single function which is potentially a bound be solving a recurrence relation
    """
    f = Function('f')
    n = symbols('n')
    solution = rsolve(f(n) - c*f(n-1) - inhom_part, f(n), init={f(0): starting_value})
    solution = solution[0][f](n)
    return solution


def __split_on_signums(expression: Expr) -> [Expr]:
    """
    For a given expression returns all expression resulting from splitting it on signum functions occurring in
    its limit, e.g. for c*sign(d - 1) returns [c*e, c*(-e)]
    """
    exps = [expression]
    n = symbols('n')
    expression_limit = limit(expression, n, oo)
    signums = get_signums_in_expression(expression_limit)
    for s in signums:
        # Choose an arbitrary symbol from the signum expression
        assert len(s.free_symbols) >= 1
        symbol = list(s.free_symbols)[0]
        new_exps = []
        for exp in exps:
            # Get rid of the signum expression by replacing it by a positve and an negative constant
            # This is done by substituting the arbitrary symbol by just the right expression s.t. things cancel out
            constant = unique_symbol('e', positive=True, real=True)
            solutions = solve(s - constant, [symbol])
            assert len(solutions) >= 1
            solution_pos = solutions[0][symbol]
            solution_neg = solution_pos.subs({constant: constant * -1})
            new_exps.append(exp.subs({symbol: solution_pos}))
            new_exps.append(exp.subs({symbol: solution_neg}))
        exps = new_exps
    return exps
