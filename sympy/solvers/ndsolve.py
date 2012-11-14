

def ndsolve(equations, func, indep_var, initial_conditions, indep_var_domain,
        params={}, solver=None, solver_options=None, return_data=False,
        interpolation_options=None):
    """
    Solve a (system of) differential equation(s) numerically.

    Returns a dictionary of the state variables, represented by
    InterpolatingFunctions, and possibly the numerical data.

    Independent variable domain - the independent variable (usually time)
    should have its domain (initial through final values) specified as either:
    [t_init, t_final, number_points]
    [t_0, t_1, t_2, ..., t_n]

    Parameters
    ==========

    equations : list
        A list of differential equations, all equal to 0.
    func : list
        A list of the functions which are dependent (state variables).
    indep_var : Symbol
        The independent variable in the system description.
    initial_conditions : dict
        A dictionary containing all of the initial conditions for each state,
        and their derivatives if necessary.
    indep_var_domain:
        See above for description.
    params : dict
        A dictionary of all the symbolic quantities which are not the
        independent variable or dependent functions. These must be substituted
        for numerical values.
    solver : str
        The desired solver to be used - the Python bindings must be installed.
    solver_options : dict
        Options to be passed to the chosen solver.
    return_data : bool
        Declare whether numerical data should be returned or not.
    interpolation_options : dict
        Options for initializing the InterpolatingFunctions.

    Examples
    ========

    >>> from sympy import symbols, Function
    >>> from sympy.solvers.ndsolve import ndsolve
    >>> a, c, t = symbols('a, c, t')
    >>> x, y, z = symbols('x, y, z', cls=Function)
    >>> x_ = Derivative(x(t), t)
    >>> y_ = Derivative(y(t), t)
    >>> z_ = Derivative(z(t), t)
    >>> system = [x_ - z(t) + y(t), y_ - 2 * a * x(t), z_ - c * x(t)]
    >>> states = [x(t), y(t), z(t)]
    >>> params = {a : 1, c : 2}
    >>> init_conds = [0, 2, 0]
    >>> ndsolve(system, states, t, init_conds, [0, 10, 100], params, 'SciPy')
    {x(t) : InterpolatingFunction(t), y(t) : InterpolatingFunction(t), z(t) :
    InterpolatingFunction(t)}

    """
    # Order - first, sanitize input
    # Next, massage equations
    #       get states, system_dict
    # Then, prep_num_ode
    #       get filename
    # call driver(filename, ics, times, opts)
    #       get data
    # create InterpolatingFunctions, and return then

    pass

def prep_num_ode(system, states, solver, output_file, fixed_args={}, args={}):
    """
    Writes a file for numerical intergration.

    Returns 'True' if a file which contains the right-hand side (rhs) function
    necessary to solve the differential equations for a system has been
    successfully been written.

    Defining the problem
        system - A dictionary which has the system defined in either a
        mass-matrix/rhs format or as a list of inseperable equations (which
        will be solved using rootfinding).
        Needs keys 'mass matrix' and 'rhs', or key 'inseperable equations'.

    Supported solvers
        solver = 'SciPy'
                 'pygsl'

    Parameters
    ==========

    system : dict
        This contains the information about the system of differential
        equations. It should have the keys outlined above.
    states : iterable
        This contains the states (and implicitly, the order of the differential
        equations if in a first order form).
    solver : str
        A string which identifies which numerical ODE solver will be used; see
        the supported solvers listed above.
    output_file : str
        The name of the file the rhs function will be written to.
    fixed_args : dict
        A dictionary of all non-independent/non-state variables which will be
        'hard-coded' into the file generated.
    args : iterable
        A container with all of the other quantities which will be defined at
        every timestep (and implicitly, the order they will be supplied in).

    Examples
    ========

    >>> from sympy import symbols, eye, Function
    >>> from sympy.solvers.ndsolve import prep_num_ode
    >>> a, c, t = symbols('a, c, t')
    >>> x, y, z = symbols('x, y, z', cls=Function)
    >>> states = [x(t), y(t)]
    >>> system = {'mass matrix' : eye(2), 'rhs' : [-y(t), 2 * x(t)]}
    >>> prep_num_ode(system, states, 'SciPy', 'my_rhs_file')
    True
    >>> states = [x(t), y(t), z(t)]
    >>> system = {'mass matrix' : eye(3), 'rhs' : [z(t) - y(t), 2 * a * x(t), c * x(t)]}
    >>> prep_num_ode(system, states, 'SciPy', 'my_rhs_file', {a : 1}, [c])
    True
    >>> # Note that this example is not 'inseperable'
    >>> x_ = Derivative(x(t), t)
    >>> y_ = Derivative(y(t), t)
    >>> z_ = Derivative(z(t), t)
    >>> system = {'inseperable equations' : [x_ - z(t) + y(t), y_ - 2 * a * x(t), z_ - c * x(t)]}
    >>> prep_num_ode(system, states, 'SciPy', 'my_rhs_file', {a : 1}, [c])

    """
    pass

def _massage_eqs(equations, states):
    """
    Private function for 'rearranging' a system of equations.

    This function takes in a list of equations and states. The equations are
    not required to be in any particular form, other than that they are all
    equal to 0. The system of equations is rewritten as:

    M * dx/dt = rhs

    which will be passed to 'prep_num_ode' as a system.

    Parameters
    ==========

    equations : iterable
        The equations which describe this system - all must be equal to 0.
    states : iterable
        The dependent functions (states) which are going to be solved for.

    """

    # system_dict needs either:
    #                          a key for 'inseperable equations'
    #                          - or -
    #                          keys 'mass matrix' and 'rhs'

    # Note that if equations are not in first order, some work might have to be
    # done with second derivatives / recasting problem in first order form,
    # e.g. if we have one equation - x''(t) + x(t) = 0, we need to rewrite that
    # as 2 1st order equations and re-identify states

    # Also, if things don't go well, 'inseperable' equations will be found -
    # see prep_num_ode

    return (states, system_dict)

def _solver_driver_scipy(filename, initial_conditions, times, solver_options):
    """
    Private function which numerically integrates equations with: SciPy

    The parameter times can be provided as:
        a beginning time, end time, and the number of evenly spaced timesteps.
        -or-
        a list of all timesteps to be evaluated at.

    Returns the times that each state has been evaluated at, as well as all the
    states at those times.
    Returns: time   - 1 list of length n
             states - 1 list of (lists of length n)

    Parameters
    ==========

    filename : str
        Name of file containing rhs function.
    intial_conditions : iterable
        The values of the states at the start of integration.
    times : iterable
        See above for description.
    solver_options : dict
        Examine the SciPy documentation for valid options to be passed to
        'odeint' as kwargs.

    """
    pass

