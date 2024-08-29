"""Small utility to obtain an absolute path with respect to the caller's file for a given input relative path."""

import inspect, os, pathlib


def this_dir(path: str = "", manual_dir: str = "") -> str:
    """Small utility to obtain an absolute path with respect to the caller's file for a given input relative path.

    Parameters
    ----------
    ``path`` :
        The relative path for which to obtain an absolute path. \
            The path must be relative to the caller's file's location.

    Returns
    -------
        The absolute path
    """
    if not manual_dir:
        where = inspect.stack()[1].filename
    else:
        where = f"{manual_dir}/dummy.placeholder"
    where_dir = pathlib.Path(where).absolute().parent
    return os.path.join(where_dir, path)
