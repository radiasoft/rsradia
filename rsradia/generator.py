import jinja2
import os
from pykern import pkio

from rsradia import _PARALLEL_RADIA_TEMPLATE, _PARALLEL_RADIA_SCRIPT, _TEMPLATE_PATH
from rsradia import _SCRIPT_DUMP_PATH, _NOTEBOOK_DUMP_PATH

# TODO: Could be more flexible about how solver arguments may be provided

_SOLVERS = {'RlxAuto': ['prec', 'maxiter', 'meth', 'ZeroM'],
            'RlxMan': [],
            'Solve': ['prec', 'maxiter', 'meth']}

def _check_solver_setup(solver, args):
    assert solver in _SOLVERS.keys(), f"solver name {solver} not understood.\nAvailable options are: " + \
                                      ', '.join(_SOLVERS.keys())
    assert type(args) == list, "Solver arguments must be given as a list. " \
                               "See Radia documentation for arguments and ordering"
    assert len(args) == len(_SOLVERS[solver]), f"Incorrect number arguments for {solver}. " + \
                                               "Expected {exp} but received {rec}.".format(exp=len(_SOLVERS[solver]),
                                                                                           rec=len(args))
    return_args = args.copy()
    for i in range(len(args)):
        if type(args[i]) == str:
            return_args[i] = "\"{}\"".format(args[i])
            
    return return_args


def generate_script(solver, args):
    """
    Create Radia parallel solve script from a template.
    :param solver: (str) One of the Radia relaxation commands. Available options are: 'RlxAuto', 'RlxMan', or 'Solve'
    :param args: (list) List of arguments to configure the solver.
    :return: None
    """
    args = _check_solver_setup(solver, args)

    template_loader = jinja2.FileSystemLoader(searchpath=_TEMPLATE_PATH)
    template_env = jinja2.Environment(loader=template_loader)
    template = template_env.get_template(_PARALLEL_RADIA_TEMPLATE)


    notebook_dump_dir, notebook_dump_file = os.path.split(_NOTEBOOK_DUMP_PATH)

    output_template = template.render(_RELAXATION_COMMAND=solver,
                                      _RELAXATION_ARGS=args,
                                      _INPUT_DUMP_FILE=notebook_dump_file,
                                      _OUTPUT_DUMP_FILE=_SCRIPT_DUMP_PATH)

    file_path = os.path.join(notebook_dump_dir, _PARALLEL_RADIA_SCRIPT)
    pkio.write_text(file_path, output_template)
