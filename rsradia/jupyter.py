import subprocess
import radia
from jupyter_radia import generator
from rsradia import _NOTEBOOK_DUMP_PATH, _SCRIPT_DUMP_PATH, _SOLVE_SCRIPT, _MPI


def _create_dump(path, radia_obj):
    with open(path, 'wb') as ff:
        mag_dump = radia.UtiDmp(radia_obj, 'bin')
        ff.write(mag_dump)

def _subprocess_solve(mpi_command, script, load_dump):
    """Run solve through MPI and return output from the shell"""
    # TODO: This command can probably be adjusted to make the source portion unnecessary
    command = f'source ~/.bashrc && {mpi_command} -n 4 python {script} {load_dump}'
    cmd = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = cmd.communicate()

    return out, err

def _check_solve_output(out, err):
    """
    Verify from shell output that the Radia solve completed satisfactorily.
    Will print any output from the shell.
    :param out: (bytes) STDOUT piped from shell
    :param err: (bytes) STDERR piped from shell
    :return: (bool) True if STDERR is empty. False if STDERR is non-empty.
    """
    if err:
        print("Error occurred during Radia solve procedure:")
        print(out.decode())
        print(err.decode())
        return False
    else:
        print("Solve finished")
        print(out.decode())
        return True

def _load_dump(path):
    """
    Load binary Radia dump object
    :param path: (str) Path to the Radia dump file
    :return: (Radia object)
    """
    with open(path, 'rb') as ff:
        result = ff.read()
        solved_magnet = radia.UtiDmpPrs(result)
        return solved_magnet

def mpi_solve(radia_obj, solver_args, solver='RlxAuto'):
    """
    Wrapper to solve magnetostatic problems in parallel from a Jupyter notebook.

    Produces a dump file and run script that will be executed with MPI to build the interaction
    matrix and run the relaxation procedure. Relaxed object is dumped by the run script and read back into the notebook.

    :param radia_obj: (Radia object): Object holding magnet(s) to perform
    :param solver_args: (list) List of arguments that will be passed to the solver command at run time.
                               All arguments must be set except for the input object. For `RlxAuto` and `RlxMan` the
                               `intrc` entry should not be included on the list, or for `Solve` the value for `obj`
                               should not be included. These entries will automatically be filled by `mpi_solve`.
    :param solver: (str) Change the Radia relaxation command to use. Options are: 'RlxAuto', 'RlxMan', or 'Solve'.
                            The default is to use 'RlxAuto'
    :return: (Radia object): Returns the magnet after performing relaxation.
    """
    if solver_args is None:
        solver_args = [1e-4, 1000, 4, 'ZeroM->False']
    generator.generate_script(solver, solver_args)

    _create_dump(_NOTEBOOK_DUMP_PATH, radia_obj)

    out, err = _subprocess_solve(_MPI, _SOLVE_SCRIPT, _NOTEBOOK_DUMP_PATH)

    solve_complete = _check_solve_output(out, err)

    if solve_complete:
        solved_magnet = _load_dump(_SCRIPT_DUMP_PATH)
        return solved_magnet
    else:
        return None
