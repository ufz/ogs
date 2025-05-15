import re
from pathlib import Path


def file_from_template(param, template, file):
    """
    Replaces strings in files masked by "%%".

    Example:
    content of template file: <parameter>%%param_name%%</parameter>
    param: {param_name: 100}
    content of resulting file: <parameter>100</parameter>

    Args:
        param (dict): key (str) is string in template without %%, value (str) replacement str
        template (str): path to template (input) file
        file (str): path to output file
    """
    param = add_prefix_suffix(param, r"%%")
    sed(param, template, file)


def add_prefix_suffix(d, prefix_suffix):
    """
    Adds a prefix and suffix to the keys in the dictionary.

    Args:
        d (dict): Input dictionary
        prefix_suffix (str): The prefix and suffix to add to each key

    Returns:
        dict: New dictionary with modified keys
    """
    new_dict = {}
    for key, value in d.items():
        new_key = f"{prefix_suffix}{key}{prefix_suffix}"
        new_dict[new_key] = str(value)
    return new_dict


def sed(replace, source, output):
    """Replaces strings in source file and writes the changes to the output file.

    In each line, replaces pattern with replace.

    Args:
        replace (dict)
            key (str): pattern to match (can be re.pattern)
            value (str): replacement str
        source  (str): input filename
        output  (str): output filename
    """

    try:
        with Path.open(source) as fin:
            try:
                with Path.open(output, "w") as fout:
                    num_replaced = 0

                    for line in fin:
                        out = line

                        for k, v in replace.items():
                            out = re.sub(k, v, out)

                        fout.write(out)

                        if out != line:
                            num_replaced += 1

                    print(f"{num_replaced} replacements made.")
            except OSError as e:
                print(f"Error writing the file '{output}': {e}")
    except OSError as e:
        print(f"Error reading the file '{source}': {e}")
