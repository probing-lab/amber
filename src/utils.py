import math
from enum import Enum, auto
from diofant import *

from mora.core import Program, get_solution as get_expected
from mora.input import LOOP_GUARD_VAR

LOG_NOTHING = 0
LOG_ESSENTIAL = 10
LOG_VERBOSE = 20
LOG_LEVEL = LOG_ESSENTIAL


class Answer(Enum):
    FALSE = auto()
    TRUE = auto()
    UNKNOWN = auto()

    def is_true(self):
        return self is Answer.TRUE

    def is_known(self):
        return self is not Answer.UNKNOWN

    def __str__(self):
        if self is Answer.TRUE:
            return "Yes"
        if self is Answer.FALSE:
            return "No"
        if self is Answer.UNKNOWN:
            return "Maybe"


__COUNTER = 0


def unique_symbol(s: str, **args):
    """
    Returns a symbol which every time has a different name
    """
    global __COUNTER
    s = symbols(s + str(__COUNTER), **args)
    __COUNTER += 1
    return s


def get_max_0(expression: Expr, n: Symbol):
    """
    Returns the maximum positive 0 of a given expression or 0 if it does not exist
    """
    n_real = symbols("n", real=True)
    try:
        exp_zeros = solve(expression.xreplace({n: n_real}), n_real)
        if exp_zeros == [{}]:
            return 0
        exp_zeros = [z[n_real] for z in exp_zeros if z[n_real].is_real]
    except NotImplementedError:
        exp_zeros = []
    exp_zeros = [math.ceil(float(z)) for z in exp_zeros] + [0]
    return max(exp_zeros)


def get_monoms(poly: Poly):
    """
    Returns for the list of monoms for a given polynomial
    """
    monoms = []
    for powers in poly.monoms():
        m = prod(var ** power for var, power in zip(poly.gens, powers))
        if m != 1:
            monoms.append(m)
    return monoms


def get_polarity(expression: Expr, n: Symbol):
    """
    Given an expression in n, returns whether or not the expression is positive and negative for some values of n
    """
    expression = simplify(expression)
    if expression.is_number:
        return expression > 0, expression < 0

    max_0 = get_max_0(expression, n)
    if max_0 > 0:
        pos = True
        neg = True
    else:
        expr_1 = expression.subs({n: 1})
        pos = expr_1 > 0
        pos = bool(pos) if pos.is_Boolean else True
        neg = expr_1 < 0
        neg = bool(neg) if neg.is_Boolean else True

    return pos, neg


def get_signums_in_expression(expression: Expr) -> [Expr]:
    """
    Given an expression it returns all expressions which occur within the signum function.
    E.g for sign(x+2)*sign(y**3) it returns [x+1, y**3]
    """
    if isinstance(expression, sign):
        return [expression.args[0]]

    signums = []
    for arg in expression.args:
        signums += get_signums_in_expression(arg)

    return signums


def get_all_monom_powers(monom: Expr) -> [Number]:
    """
    Returns the degrees of all variables in a given monomial in a list
    """
    monom = monom.as_poly(monom.free_symbols)
    return list(monom.degree_list())


def monom_is_deterministic(monom: Expr, program: Program):
    """
    Returns true iff a given monomial is deterministic, that means all variables in the monomial are deterministic
    """
    variables_deterministic = [not program.updates[m].is_probabilistic for m in monom.free_symbols]
    return all(variables_deterministic)


def divide_monom_powers_by(monom: Expr, divisor):
    """
    Returns the given monom where all powers a divided by divisor
    """
    monom = monom.as_poly(monom.free_symbols)
    powers = monom.monoms()[0]
    vars_with_powers = [v ** (p // divisor) for v, p in zip(monom.gens, powers)]
    return prod(vars_with_powers)


def separate_rvs_from_monom(monom: Expr, program: Program):
    """
    Given a monomial returns a list of all random variables (together with their powers)
    it contains and the remaining monomial
    """
    if not program.contains_rvs:
        return [], monom

    monom = monom.as_poly(monom.free_symbols)
    powers = monom.monoms()[0]
    vars_with_powers = [(v, p) for v, p in zip(monom.gens, powers)]
    m = sympify(1)
    rvs = []
    for v, p in vars_with_powers:
        if program.updates[v].is_random_var and not hasattr(program.updates[v], "branches"):
            rvs.append((v, p))
        else:
            m *= v ** p
    return rvs, m


def set_log_level(log_level):
    global LOG_LEVEL
    LOG_LEVEL = log_level


def log(message, level):
    """
    Logs a message depending on the log level
    """
    if level <= LOG_LEVEL:
        print(message)


def amber_limit(expr, n):
    if n not in expr.free_symbols:
        return expr

    return limit(expr, n, oo)


def flatten_substitution_choices(subs_choices):
    """
    For a given dict {expr: (expr1, expr2)} returns a list of all possible substitution arising from choosing to subs
    expr by expr1 or expr2.
    """
    subs_choices = subs_choices.copy()
    if not subs_choices:
        return [{}]

    result = []
    expr = next(iter(subs_choices.keys()))
    choice1, choice2 = subs_choices.pop(expr)
    remaining_choices_flat = flatten_substitution_choices(subs_choices)
    for c in remaining_choices_flat:
        c1 = c.copy()
        c1[expr] = choice1
        result.append(c1)
        if choice1 != choice2:
            c2 = c.copy()
            c2[expr] = choice2
            result.append(c2)
    return result


def substitute_deterministic_variables(expr, program: Program):
    """
    Substitutes deterministic variables in a given expression with their representation in n.
    """
    for symbol, update in program.updates.items():
        if str(symbol) is not LOOP_GUARD_VAR and not update.is_probabilistic:
            closed_form = get_expected(program, symbol.as_poly(program.variables))
            expr = expr.xreplace({symbol: closed_form})
    return expr
