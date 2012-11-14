from sympy.core.function import Function


class InterpolatingFunction(Function):
    """
    A symbolic wrapper to approximate a function.

    This probably shouldn't be constructed by a user, but instead returned from
    a function and used to symbolically 'mask' numerical data.

    When called with symbolic arguements, a SymPy Function is returned. But
    when called with all numerical arguements, a numerical value is returned.

    Examples
    ========

    >>> # somehow it has been created
    >>> f
    InterpolatingFunction(t)
    >>> f.subs(t, 5)
    2.341
    >>> 2 * f
    2 * InterpolatingFunction(t)
    >>> (2 * f).subs(t, 5)
    4.682

    """

    # I'm not sure where the numerical data should be stored - in args probably?

    def __init__(self, variables, data):
        pass

    def __call__(self, args):
        defined = list(self._args[0])
        passed = args
        out = [passed[i] if passed[i] != defined[i] else defined[i] for i in
                range(len(defined))]
        if min([i.is_Number for i in out]):
            # still symbolic, return new InterpolatingFunction with mix of
            # symbolic and numberic arguements, e.g. f(x, 2, v)
            pass
        else:
            # no longer symbolic, return numerical value
            pass


