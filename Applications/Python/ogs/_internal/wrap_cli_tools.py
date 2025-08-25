import subprocess

from . import OGS_USE_PATH
from .binaries_list import binaries_list
from .get_bin_dir import get_bin_dir
from .provide_ogs_cli_tools_via_wheel import ogs_with_args


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
        >>> cli.ogs(help=True) # flags without any additional value can be set to True
        ...

        A more useful example. The following will create a line mesh:

        >>> outfile = "line.vtu"
        >>> cli.generateStructuredMesh(e="line", lx=1, nx=10, o=outfile)
        """

        cmdline = CLI._get_cmdline("ogs", *args, **kwargs)

        if OGS_USE_PATH:
            print("OGS_USE_PATH is true: ogs from $PATH is used!")
            return subprocess.call(cmdline)

        return ogs_with_args(cmdline)

    @staticmethod
    def _format_kv(kwargs):
        for key, v in kwargs.items():
            # Convert None to True
            if v is None:
                # TODO: Remove after 2025/08
                print(
                    f"Deprecation warning: Setting {v} to `None` is deprecated, set to `True` instead!"
                )
                v = True  # noqa: PLW2901

            # If value is False then ignore
            if isinstance(v, bool) and not v:
                continue

            if len(key) == 1:
                yield f"-{key}"
            else:
                yield f"--{key}"

            # Pass value if not bool
            if not isinstance(v, bool):
                yield f"{v}"

    @staticmethod
    def _get_cmdline(cmd, *args, **kwargs):
        str_kwargs = list(CLI._format_kv(kwargs))
        return [cmd] + str_kwargs + list(args)

    @staticmethod
    def _get_run_cmd(attr):
        def run_cmd(*args, **kwargs):
            cmd = get_bin_dir() / attr
            if OGS_USE_PATH:
                cmd = attr
            cmdline = CLI._get_cmdline(cmd, *args, **kwargs)
            return subprocess.call(cmdline)

        # TODO: Only arguments with underscores work. Arguments with dashes not.
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
            >>> cli.ogs(help=True) # flags without any additional value can be set to True
            ...

            A more useful example. The following will create a line mesh:

            >>> outfile = "line.vtu"
            >>> cli.generateStructuredMesh(e="line", lx=1, nx=10, o=outfile)
            """

        return run_cmd


cli = CLI()
