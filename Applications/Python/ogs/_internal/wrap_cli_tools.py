import subprocess
import os
from .provide_ogs_cli_tools_via_wheel import binaries_list

OGS_BIN_DIR = os.path.join(os.path.join(os.path.dirname(__file__), "..", "bin"))


class CLI:
    def __init__(self):
        for b in binaries_list:
            if b == "ogs":
                continue  # command provided separately
            setattr(self, b, CLI._get_run_cmd(b))

    def ogs(self, *args, **kwargs):
        """
        This function wraps the commandline tool ogs for easy use from Python.

        It returns the return code of the commandline tool.

        The entries of args are passed as is to the commandline tool.
        The entries of kwargs are transformed: one-letter keys get a single
        dash as a prefix, multi-letter keys are prefixed with two dashes,
        underscores are replaced with dashes.

        Thereby, commandline tools can be used in a "natural" way from Python, e.g.:

        >>> cli = CLI()
        >>> cli.ogs("--help") # prints a help text
        ...
        >>> cli.ogs(help=None) # special, does the same
        ...

        A more useful example. The following will create a line mesh:

        >>> outfile = "line.vtu"
        >>> cli.generateStructuredMesh(e="line", lx=1, nx=10, o=outfile)
        """

        import ogs.simulator as sim

        cmdline = CLI._get_cmdline("ogs", *args, **kwargs)
        return_code = sim.initialize(cmdline)

        # map mangled TCLAP status to usual exit status
        if return_code == 3:  # EXIT_ARGPARSE_FAILURE
            sim.finalize()
            return 1  # EXIT_FAILURE
        elif return_code == 2:  # EXIT_ARGPARSE_EXIT_OK
            sim.finalize()
            return 0  # EXIT_SUCCESS

        if return_code != 0:
            sim.finalize()
            return return_code

        return_code = sim.executeSimulation()
        sim.finalize()
        return return_code

    @staticmethod
    def _format_kv(kwargs):
        for k, v in kwargs.items():
            key = k.replace("_", "-")

            if len(k) == 1:
                yield f"-{key}"
            else:
                yield f"--{key}"

            if v is not None:
                yield f"{v}"

    @staticmethod
    def _get_cmdline(cmd, *args, **kwargs):
        str_kwargs = [arg for arg in CLI._format_kv(kwargs)]
        return [cmd] + str_kwargs + list(args)

    @staticmethod
    def _get_run_cmd(attr):
        def run_cmd(*args, **kwargs):
            # TODO provide override via os.environ?
            cmd = os.path.join(OGS_BIN_DIR, attr)
            cmdline = CLI._get_cmdline(cmd, *args, **kwargs)
            return subprocess.call(cmdline)

        run_cmd.__doc__ = f"""
            This function wraps the commandline tool {attr} for easy use from Python.

            It returns the return code of the commandline tool.

            The entries of args are passed as is to the commandline tool.
            The entries of kwargs are transformed: one-letter keys get a single
            dash as a prefix, multi-letter keys are prefixed with two dashes,
            underscores are replaced with dashes.

            Thereby, commandline tools can be used in a "natural" way from Python, e.g.:

            >>> cli = CLI()
            >>> cli.{attr}("--help") # prints a help text
            ...
            >>> cli.{attr}(help=None) # special, does the same
            ...

            A more useful example. The following will create a line mesh:

            >>> outfile = "line.vtu"
            >>> cli.generateStructuredMesh(e="line", lx=1, nx=10, o=outfile)
            """

        return run_cmd


cli = CLI()
