import importlib.util
import json
import subprocess
import sysconfig
from importlib.metadata import Distribution, PackageNotFoundError
from pathlib import Path

from . import OGS_USE_PATH
from .provide_ogs_cli_tools_via_wheel import binaries_list, ogs_with_args

pkg_is_editable = False
if importlib.util.find_spec("ogs") is not None:
    try:
        dist = Distribution.from_name("ogs")
        direct_url = dist.read_text("direct_url.json")
        if direct_url:
            pkg_is_editable = (
                json.loads(direct_url).get("dir_info", {}).get("editable", False)
            )
    except PackageNotFoundError:
        # ogs is importable but not an installed distribution (e.g. source tree, dev mode)
        pkg_is_editable = False
    except FileNotFoundError:
        # distribution is found but no direct_url.json
        pass

if pkg_is_editable:
    site_packages_path = sysconfig.get_paths()["purelib"]
    OGS_BIN_DIR = Path(site_packages_path) / "bin"
else:
    # Here, we assume that this script is installed, e.g., in a virtual environment
    # alongside a "bin" directory.
    OGS_BIN_DIR = Path(__file__).parent.parent.parent / "bin"  # installed wheel

if not OGS_BIN_DIR.exists():
    OGS_BIN_DIR = OGS_BIN_DIR.parent  # build directory


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
            cmd = OGS_BIN_DIR / attr
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
