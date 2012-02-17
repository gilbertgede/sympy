__all__ = ['cross',
           'dot',
           'express',
           'outer',
           'inertia',
           'mechanics_printing',
           'mprint',
           'mpprint',
           'mlatex',
           'kinematic_equations',
           'code_output']

from sympy.physics.mechanics.essential import (Vector, Dyadic, ReferenceFrame,
                                               MechanicsStrPrinter,
                                               MechanicsPrettyPrinter,
                                               MechanicsLatexPrinter,
                                               dynamicsymbols)
from sympy.physics.mechanics.kane import Kane
from sympy import sympify, diff, sin, cos, Matrix

def cross(vec1, vec2):
    """Cross product convenience wrapper for Vector.cross(): \n"""
    if not isinstance(vec1, (Vector, Dyadic)):
        raise TypeError('Cross product is between two vectors')
    return vec1 ^ vec2
cross.__doc__ += Vector.cross.__doc__

def dot(vec1, vec2):
    """Dot product convenience wrapper for Vector.dot(): \n"""
    if not isinstance(vec1, (Vector, Dyadic)):
        raise TypeError('Dot product is between two vectors')
    return vec1 & vec2
dot.__doc__ += Vector.dot.__doc__

def express(vec, frame, frame2=None):
    """Express convenience wrapper for Vector.express(): \n"""
    if not isinstance(vec, (Vector, Dyadic)):
        raise TypeError('Can only express Vectors')
    if isinstance(vec, Vector):
        return vec.express(frame)
    else:
        return vec.express(frame, frame2)

express.__doc__ += Vector.express.__doc__

def outer(vec1, vec2):
    """Outer prodcut convenience wrapper for Vector.outer():\n"""
    if not isinstance(vec1, Vector):
        raise TypeError('Outer product is between two Vectors')
    return vec1 | vec2
outer.__doc__ += Vector.express.__doc__

def inertia(frame, ixx, iyy, izz, ixy=0, iyz=0, izx=0):
    """Simple way to create inertia Dyadic object.

    If you don't know what a Dyadic is, just treat this like the inertia
    tensor.  Then, do the easy thing and define it in a body-fixed frame.

    Parameters
    ==========

    frame : ReferenceFrame
        The frame the inertia is defined in
    ixx : Sympifyable
        the xx element in the inertia dyadic
    iyy : Sympifyable
        the yy element in the inertia dyadic
    izz : Sympifyable
        the zz element in the inertia dyadic
    ixy : Sympifyable
        the xy element in the inertia dyadic
    iyz : Sympifyable
        the yz element in the inertia dyadic
    izx : Sympifyable
        the zx element in the inertia dyadic

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame, inertia
    >>> N = ReferenceFrame('N')
    >>> inertia(N, 1, 2, 3)
    (N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z)

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Need to define the inertia in a frame')
    ol  = sympify(ixx) * (frame.x | frame.x)
    ol += sympify(ixy) * (frame.x | frame.y)
    ol += sympify(izx) * (frame.x | frame.z)
    ol += sympify(ixy) * (frame.y | frame.x)
    ol += sympify(iyy) * (frame.y | frame.y)
    ol += sympify(iyz) * (frame.y | frame.z)
    ol += sympify(izx) * (frame.z | frame.x)
    ol += sympify(iyz) * (frame.z | frame.y)
    ol += sympify(izz) * (frame.z | frame.z)
    return ol

def mechanics_printing():
    """Sets up interactive printing for mechanics' derivatives.

    The main benefit of this is for printing of time derivatives;
    instead of displaying as Derivative(f(t),t), it will display f'
    This is only actually needed for when derivatives are present and are not
    in a physics.mechanics object.

    Examples
    ========

    >>> # 2 lines below are for tests to function properly
    >>> import sys
    >>> sys.displayhook = sys.__displayhook__
    >>> from sympy import Function, Symbol, diff
    >>> from sympy.physics.mechanics import mechanics_printing
    >>> f = Function('f')
    >>> t = Symbol('t')
    >>> x = Symbol('x')
    >>> diff(f(t), t)
    Derivative(f(t), t)
    >>> mechanics_printing()
    >>> diff(f(t), t)
    f'
    >>> diff(f(x), x)
    Derivative(f(x), x)
    >>> # 2 lines below are for tests to function properly
    >>> import sys
    >>> sys.displayhook = sys.__displayhook__

    """

    import sys
    sys.displayhook = mprint

def mprint(expr, **settings):
    r"""Function for printing of expressions generated in mechanics.

    Extends SymPy's StrPrinter; mprint is equivalent to:
    print sstr()
    mprint takes the same options as sstr.

    Parameters
    ==========

    expr : valid sympy object
        SymPy expression to print
    settings : args
        Same as print for SymPy

    Examples
    ========

    >>> from sympy.physics.mechanics import mprint, dynamicsymbols
    >>> u1 = dynamicsymbols('u1')
    >>> print(u1)
    u1(t)
    >>> mprint(u1)
    u1

    """

    pr = MechanicsStrPrinter(settings)
    outstr = pr.doprint(expr)

    import __builtin__
    if (outstr != 'None'):
        __builtin__._ = outstr
        print(outstr)

def mpprint(expr, **settings):
    r"""Function for pretty printing of expressions generated in mechanics.

    Mainly used for expressions not inside a vector; the output of running
    scripts and generating equations of motion. Takes the same options as
    SymPy's pretty_print(); see that function for more information.

    Parameters
    ==========

    expr : valid sympy object
        SymPy expression to pretty print
    settings : args
        Same as pretty print

    Examples
    ========

    Use in the same way as pprint

    """

    mp = MechanicsPrettyPrinter(settings)
    print(mp.doprint(expr))

def mlatex(expr, **settings):
    r"""Function for printing latex representation of mechanics objects.

    For latex representation of Vectors, Dyadics, and dynamicsymbols. Takes the
    same options as SymPy's latex(); see that function for more information;

    Parameters
    ==========

    expr : valid sympy object
        SymPy expression to represent in LaTeX form
    settings : args
        Same as latex()

    Examples
    ========

    >>> from sympy.physics.mechanics import mlatex, ReferenceFrame
    >>> N = ReferenceFrame('N')
    >>> mlatex(N.x + N.y)
    '\\mathbf{\\hat{n}_x} + \\mathbf{\\hat{n}_y}'

    """

    return MechanicsLatexPrinter(settings).doprint(expr)

def kinematic_equations(speeds, coords, rot_type, rot_order=''):
    """Gives equations relating the qdot's to u's for a rotation type.

    Supply rotation type and order as in orient. Speeds are assumed to be
    body-fixed; if we are defining the orientation of B in A using by rot_type,
    the angular velocity of B in A is assumed to be in the form: speed[0]*B.x +
    speed[1]*B.y + speed[2]*B.z

    Parameters
    ==========

    speeds : list of length 3
        The body fixed angular velocity measure numbers.
    coords : list of length 3 or 4
        The coordinates used to define the orientation of the two frames.
    rot_type : str
        The type of rotation used to create the equations. Body, Space, or
        Quaternion only
    rot_order : str
        If applicable, the order of a series of rotations.

    Examples
    ========

    >>> from sympy.physics.mechanics import dynamicsymbols
    >>> from sympy.physics.mechanics import kinematic_equations, mprint
    >>> u1, u2, u3 = dynamicsymbols('u1 u2 u3')
    >>> q1, q2, q3 = dynamicsymbols('q1 q2 q3')
    >>> mprint(kinematic_equations([u1,u2,u3], [q1,q2,q3], 'body', '313'),
    ...     order=None)
    [-(u1*sin(q3) + u2*cos(q3))/sin(q2) + q1', -u1*cos(q3) + u2*sin(q3) + q2', (u1*sin(q3) + u2*cos(q3))*cos(q2)/sin(q2) - u3 + q3']

    """

    # Code below is checking and sanitizing input
    approved_orders = ('123', '231', '312', '132', '213', '321', '121', '131',
                       '212', '232', '313', '323', '1', '2', '3', '')
    rot_order = str(rot_order).upper() # Now we need to make sure XYZ = 123
    rot_type  = rot_type.upper()
    rot_order = [i.replace('X', '1') for i in rot_order]
    rot_order = [i.replace('Y', '2') for i in rot_order]
    rot_order = [i.replace('Z', '3') for i in rot_order]
    rot_order = ''.join(rot_order)

    if not isinstance(speeds,(list, tuple)):
        raise TypeError('Need to supply speeds in a list')
    if len(speeds) != 3:
        raise TypeError('Need to supply 3 body-fixed speeds')
    if not isinstance(coords,(list, tuple)):
        raise TypeError('Need to supply coordinates in a list')
    if rot_type.lower() in ['body', 'space']:
        if rot_order not in approved_orders:
            raise ValueError('Not an acceptable rotation order')
        if len(coords) != 3:
            raise ValueError('Need 3 coordinates for body or space')
        # Actual hard-coded kinematic differential equations
        q1, q2, q3 = coords
        q1d, q2d, q3d = [diff(i, dynamicsymbols._t) for i in coords]
        w1, w2, w3 = speeds
        s1, s2, s3 = [sin(q1), sin(q2), sin(q3)]
        c1, c2, c3 = [cos(q1), cos(q2), cos(q3)]
        if rot_type.lower() == 'body':
            if rot_order == '123':
                return [q1d - (w1 * c3 - w2 * s3) / c2, q2d - w1 * s3 - w2 *
                        c3, q3d - (-w1 * c3 + w2 * s3) * s2 / c2 - w3]
            if rot_order == '231':
                return [q1d - (w2 * c3 - w3 * s3) / c2, q2d - w2 * s3 - w3 *
                        c3, q3d - w1 - (- w2 * c3 + w3 * s3) * s2 / c2]
            if rot_order == '312':
                return [q1d - (-w1 * s3 + w3 * c3) / c2, q2d - w1 * c3 - w3 *
                        s3, q3d - (w1 * s3 - w3 * c3) * s2 / c2 - w2]
            if rot_order == '132':
                return [q1d - (w1 * c3 + w3 * s3) / c2, q2d + w1 * s3 - w3 *
                        c3, q3d - (w1 * c3 + w3 * s3) * s2 / c2 - w2]
            if rot_order == '213':
                return [q1d - (w1 * s3 + w2 * c3) / c2, q2d - w1 * c3 + w2 *
                        s3, q3d - (w1 * s3 + w2 * c3) * s2 / c2 - w3]
            if rot_order == '321':
                return [q1d - (w2 * s3 + w3 * c3) / c2, q2d - w2 * c3 + w3 *
                        s3, q3d - w1 - (w2 * s3 + w3 * c3) * s2 / c2]
            if rot_order == '121':
                return [q1d - (w2 * s3 + w3 * c3) / s2, q2d - w2 * c3 + w3 *
                        s3, q3d - w1 + (w2 * s3 + w3 * c3) * c2 / s2]
            if rot_order == '131':
                return [q1d - (-w2 * c3 + w3 * s3) / s2, q2d - w2 * s3 - w3 *
                        c3, q3d - w1 - (w2 * c3 - w3 * s3) * c2 / s2]
            if rot_order == '212':
                return [q1d - (w1 * s3 - w3 * c3) / s2, q2d - w1 * c3 - w3 *
                        s3, q3d - (-w1 * s3 + w3 * c3) * c2 / s2 - w2]
            if rot_order == '232':
                return [q1d - (w1 * c3 + w3 * s3) / s2, q2d + w1 * s3 - w3 *
                        c3, q3d + (w1 * c3 + w3 * s3) * c2 / s2 - w2]
            if rot_order == '313':
                return [q1d - (w1 * s3 + w2 * c3) / s2, q2d - w1 * c3 + w2 *
                        s3, q3d + (w1 * s3 + w2 * c3) * c2 / s2 - w3]
            if rot_order == '323':
                return [q1d - (-w1 * c3 + w2 * s3) / s2, q2d - w1 * s3 - w2 *
                        c3, q3d - (w1 * c3 - w2 * s3) * c2 / s2 - w3]
        if rot_type.lower() == 'space':
            if rot_order == '123':
                return [q1d - w1 - (w2 * s1 + w3 * c1) * s2 / c2, q2d - w2 *
                        c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / c2]
            if rot_order == '231':
                return [q1d - (w1 * c1 + w3 * s1) * s2 / c2 - w2, q2d + w1 *
                        s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / c2]
            if rot_order == '312':
                return [q1d - (w1 * s1 + w2 * c1) * s2 / c2 - w3, q2d - w1 *
                        c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / c2]
            if rot_order == '132':
                return [q1d - w1 - (-w2 * c1 + w3 * s1) * s2 / c2, q2d - w2 *
                        s1 - w3 * c1, q3d - (w2 * c1 - w3 * s1) / c2]
            if rot_order == '213':
                return [q1d - (w1 * s1 - w3 * c1) * s2 / c2 - w2, q2d - w1 *
                        c1 - w3 * s1, q3d - (-w1 * s1 + w3 * c1) / c2]
            if rot_order == '321':
                return [q1d - (-w1 * c1 + w2 * s1) * s2 / c2 - w3, q2d - w1 *
                        s1 - w2 * c1, q3d - (w1 * c1 - w2 * s1) / c2]
            if rot_order == '121':
                return [q1d - w1 + (w2 * s1 + w3 * c1) * c2 / s2, q2d - w2 *
                        c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / s2]
            if rot_order == '131':
                return [q1d - w1 - (w2 * c1 - w3 * s1) * c2 / s2, q2d - w2 *
                        s1 - w3 * c1, q3d - (-w2 * c1 + w3 * s1) / s2]
            if rot_order == '212':
                return [q1d - (-w1 * s1 + w3 * c1) * c2 / s2 - w2, q2d - w1 *
                        c1 - w3 * s1, q3d - (w1 * s1 - w3 * c1) / s2]
            if rot_order == '232':
                return [q1d + (w1 * c1 + w3 * s1) * c2 / s2 - w2, q2d + w1 *
                        s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / s2]
            if rot_order == '313':
                return [q1d + (w1 * s1 + w2 * c1) * c2 / s2 - w3, q2d - w1 *
                        c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / s2]
            if rot_order == '323':
                return [q1d - (w1 * c1 - w2 * s1) * c2 / s2 - w3, q2d - w1 *
                        s1 - w2 * c1, q3d - (-w1 * c1 + w2 * s1) / s2]
    elif rot_type.lower() == 'quaternion':
        if rot_order != '':
            raise ValueError('Cannot have rotation order for quaternion')
        if len(coords) != 4:
            raise ValueError('Need 4 coordinates for quaternion')
        # Actual hard-coded kinematic differential equations
        e0, e1, e2, e3 = coords
        w = Matrix(speeds + [0])
        E = Matrix([[e0, -e3, e2, e1], [e3, e0, -e1, e2], [-e2, e1, e0, e3],
            [-e1, -e2, -e3, e0]])
        edots = Matrix([diff(i, dynamicsymbols._t) for i in [e1, e2, e3, e0]])
        return list(edots.T - 0.5 * w.T * E.T)
    else:
        raise ValueError('Not an approved rotation type for this function')


class code_output:
    """
    Class designed to generate code for numerical integration from a
    sympy.physics.mechanics class which has generated equations of motion.

    This is created by taking in an mechanics class which has generated
    equations of motion (currently just the Kane class), and, with some options,
    produces a file with simulation settings and another file which provides
    a function which returns the time derivatives of the coordinates and speeds
    (the right hand side function). The file containing the r.h.s. function can
    be either run directly, or have its rhs function imported to another Python
    script.

    Code output options currently include:
        NumPy/SciPy/Matplotlib

    Parameters
    ==========

    eom_class : Kane (currently only supported class)
        A class which has generated equations of motion already.
    fname : string
        The name of the file which will hold the r.h.s. function; the settings
        name file will also be named off of this.

    """

    def __init__(self, eom_class, fname):
        if not isinstance(eom_class, Kane):
            raise TypeError('Need to supply a valid EoM generating object')

        self._fname = fname
        self._settings_fname = None

        # States
        self._q = eom_class._q + eom_class._qdep
        self._u = eom_class._u + eom_class._udep
        self._qu = self._q + self._u
        # set initial conditions to 0
        self._init_conds = [0] * len(self._qu)

        # The actual equations
        self._mm = eom_class.mass_matrix_full.subs(eom_class.kindiffdict())
        self._f = eom_class.forcing_full.subs(eom_class.kindiffdict())

        # Format these nicely to put into _find functions
        mm_list = self._mm.tolist()
        mm_list = [item for sublist in mm_list for item in sublist]
        f_list = self._f.tolist()
        f_list = [item for sublist in f_list for item in sublist]
        templ = mm_list + f_list
        # parameters which need numerical values (basically, sympy symbols)
        self._params = eom_class._find_othersymbols(templ)
        self._param_vals = [0] * len(self._params)
        # specified quantities (things like forces, time-varying masses, etc.)
        self._specif = eom_class._find_dynamicsymbols(templ, self._q + self._u)

        # default integration options
        self._odeint_time = [0, 10, 100]
        self._odeint_opts = {}

    def get_parameters(self):
        """
        This returns the parameters which need to have values supplied.

        If nothing is supplied, all parameters will be set to 0 in the settings
        file.
        """
        return self._params

    def give_parameters(self, param_list):
        if len(self._params) != len(param_list):
            raise ValueError('Supply either all or none of the parameters')
        self._param_vals = param_list

    def get_initialconditions(self):
        """
        Returns a list of states which need to have initial values supplied.
        """
        return self._q + self._u

    def give_initialconditions(self, ic_list):
        """
        Supply the initial conditions for simulation.

        If nothing is supplied, all initial conditions will be set to 0 in the
        settings file. Also, only all or none of the I.C.'s can be specified.

        Parameters
        ==========

        ic_list : list
            A list of the values from the initial conditions, set in the order
            specified by .get_initialconditions().

        """

        if len(self._init_conds) != len(ic_list):
            raise ValueError('Supply either all or none of the initial'
                             ' conditions')
        self._init_conds = ic_list

    def give_time_int(self, time_list, ode_opts=None):
        """
        Supply a list of timesteps and a dictionary with integration settings.

        Parameters
        ==========

        time_list : list-like
            A list which has t_initial, t_final, n_timesteps
        ode_opts : dictionary
            A dictionary which takes in the integration options for
            scipy.integrate.odeint. Supply as {'kwarg': val}

        """

        self._odeint_time = time_list
        self._odeint_opts = ode_opts

    def write_settings_file(self, fname=None):
        """
        Write the simulation settings to a file.

        Parameters
        ==========

        fname : string
            If no filename for the settings file is supplied, a default one
            will be created. If a filename is supplied, than it will be the
            default one for the right hand side file.

        """

        mstr = MechanicsStrPrinter().doprint
        old_str = dynamicsymbols._str
        dynamicsymbols._str = 'p'
        if fname is None:
            fname = self._fname + '_settings'
        self._settings_fname = fname + '.csv'
        f = open("%s.csv" % fname, 'w')
        # Now write settings to file
        qu = self._q + self._u
        params = self._params
        for i, v in enumerate(qu):
            val = self._init_conds[i]
            f.write("%s; %f\n" % (mstr(v), val))
        for i, v in enumerate(params):
            val = self._param_vals[i]
            f.write("%s; %f\n" % (mstr(v), val))
        f.write("%s; %s\n" % ('time_list', self._odeint_time))
        f.write("%s; %s\n" % ('odeint_opts', self._odeint_opts))
        f.close()
        dynamicsymbols._str = old_str

    def write_rhs_file(self, code_type=None):
        """
        Write the file to run the integration.

        Choose what language you'll use for numerical integration. Right now,
        only SciPy is supported.

        Parameters
        ==========

        code_type : string
            The type of code output you want.

        """

        if code_type == None:
            raise ValueError('A type of code output needs to be supplied')

        if code_type.lower() == 'scipy':
            self._scipy()
        else:
            raise ValueError('Unsupported Output Type')

    def _scipy(self):
        if self._settings_fname == None:
            self.write_settings_file()
        # Printer needs different settings for file output
        mstr = MechanicsStrPrinter().doprint
        old_str = dynamicsymbols._str
        dynamicsymbols._str = 'p'
        f = open("%s.py" % self._fname, 'w')
        f.write('#!/usr/bin/env python\n\n')
        f.write('import sys\n')
        f.write('from numpy import *\n')
        f.write('acos = arccos\n')
        f.write('asin = arcsin\n')
        f.write('atan = arctan\n\n')
        # Writing the rhs function
        f.write('def rhs(y, t, params):\n')
        qu = self._q + self._u
        params = self._params
        l1 = '    ' + mstr(qu) + ' = y\n'
        f.write(l1)
        l2 = '    ' + mstr(params) + ' = params\n'
        f.write(l2)

        f.write('    M = zeros(%s)\n' % str(self._mm.shape))
        f.write('    f = zeros(%s)\n' % str(self._f.shape))

        ma, mb = self._mm.shape
        fa = len(self._f)
        for i in range(ma):
            for j in range(mb):
                f.write('    M[%s][%s] = %s\n' %
                        (i, j, mstr(self._mm.extract([i],[j]))[1:-1]))
        for i in range(fa):
            f.write('    f[%s] = %s\n' % (i, mstr(self._f[i])))
        f.write('    return linalg.solve(M, f).flatten().tolist()\n\n')

        # Writing the main function, which does some time integration
        temp = ('def main(argv=None):\n'
                '    from scipy.integrate import odeint\n'
                '    import matplotlib.pyplot as plt\n\n'
                '    if argv is None:\n'
                '        argv = sys.argv\n'
                "    # load settings file\n"
                '    if len(argv) == 2:\n'
                '        settings_fname = argv[1]\n'
                '    elif len(argv) == 1:\n'
                "        settings_fname = '%s'\n" % self._settings_fname +
                "    else:\n"
                "        raise ValueError('Maximum of 1 settings file')\n\n"
                "    f = open(settings_fname, 'r')\n"
                "    # rearrange settings file (odeint opts, ic's, params)\n\n")
        temp += (
                "    param_names = %s\n" % [mstr(i) for i in params] +
                "    state_names = %s\n" % [mstr(i) for i in qu])
        temp += (
                "    params_dict = {}\n"
                "    state_dict = {}\n"
                "    odeint_opts = {}\n"
                "    for line in f:\n"
                "        line = line.split('#')[0].strip()\n"
                "        if not line:\n"
                "            continue\n\n"
                "        name, value = line.split(';')\n"
                "        if name in param_names:\n"
                "            params_dict.update({name: float(value)})\n"
                "        elif name in state_names:\n"
                "            state_dict.update({name: float(value)})\n"
                "        elif name == 'time_list':\n"
                "            time_list = value\n"
                "        elif name == 'odeint_opts':\n"
                "            val = value.split(',')\n"
                "            odeint_opts = dict(zip(val[::2], val[1::2]))\n\n"
                "    times = time_list.strip()[1:-1].split(',')\n"
                "    t0, tf, nsteps = [float(i) for i in times]\n"
                "    params = [params_dict[i] for i in param_names]\n"
                "    y0 = [state_dict[i] for i in state_names]\n"
                "    t = linspace(t0, tf, nsteps)\n"
                "\n\n\n    # perform integration\n\n"
                '    n = len(y0)\n'
                '    y = odeint(rhs, y0, t, args=(params,), **odeint_opts)\n\n'
                '    plt.plot(t, y)\n'
                "    plt.xlabel('Time (s)')\n"
                "    plt.ylabel('States')\n"
                "    plt.legend(%s)\n" % str([mstr(i) for i in qu]) +
                '    plt.show()\n\n')
        f.write(temp)

        # Writing the main function call, which may or may not be called
        f.write('if __name__ == "__main__":\n    sys.exit(main())')

        # Return printer to default settings
        dynamicsymbols._str = old_str
