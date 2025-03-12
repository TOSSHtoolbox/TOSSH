import matlab.engine
import numpy as np
import pandas as pd
from datetime import datetime, timedelta


def run_tossh_function(data_path, function_name, eng, *args, **kwargs):
    """
    Calls a MATLAB function with specified arguments and options.

    Parameters:
    - data_path: Path to MATLAB toolbox or folder containing the function.
    - function_name: Name of the MATLAB function to call.
    - eng: MATLAB engine instance.
    - *args: Positional arguments to pass to the MATLAB function.
    - **kwargs: Additional options:
        - 'nargout': Number of output arguments to retrieve from the MATLAB function.
        - Other key-value pairs are passed as varargin to the MATLAB function.

    Returns:
    - Result(s) from the MATLAB function call.
    """
    # Add the specified path to MATLAB's search path
    eng.addpath(eng.genpath(data_path), nargout=0)

    # Convert positional arguments to MATLAB-compatible format
    args = [convert_to_matlab(arg) for arg in args]

    # Extract 'nargout' from kwargs (default is 1)
    nargout = kwargs.pop('nargout', 1)

    # Convert remaining kwargs into a list of key-value pairs for varargin
    varargin = []
    for key, value in kwargs.items():
        varargin.append(key)
        varargin.append(convert_to_matlab(value))

    # Get the MATLAB function handle
    func = getattr(eng, function_name)

    # Call the MATLAB function with args and varargin
    return func(*args, *varargin, nargout=nargout)


def convert_to_matlab(obj):
    """
    Converts Python objects into MATLAB-compatible types.

    Parameters:
    - obj: Python object (e.g., DataFrame, NumPy array, list, bool).

    Returns:
    - MATLAB-compatible object.
    """
    if isinstance(obj, bool):
        return matlab.logical([obj])

    elif isinstance(obj, pd.DataFrame):
        return {col: convert_to_matlab(obj[col].values) for col in obj}

    elif isinstance(obj, np.ndarray):
        if np.issubdtype(obj.dtype, np.datetime64):
            # Convert numpy datetime64 to Python datetime
            python_dates = obj.astype('M8[ms]').astype(datetime)

            # Convert Python datetime to MATLAB datenum
            matlab_datenum = [
                (d - datetime(1, 1, 1)).days + (d - datetime(1, 1, 1)).seconds / 86400 + 367
                for d in python_dates
            ]
            return matlab.double(matlab_datenum, size=(len(matlab_datenum), 1))

        elif np.issubdtype(obj.dtype, np.bool_):
            return matlab.logical(obj.tolist())

        return matlab.double(obj.flatten().tolist(), size=(len(obj), 1))

    elif isinstance(obj, (list, tuple)):
        if all(isinstance(x, bool) for x in obj):
            return matlab.logical(obj)
        return matlab.double(obj, size=(len(obj), 1))

    elif isinstance(obj, (int, float)):
        return matlab.double([obj])

    elif isinstance(obj, str):
        return obj

    return obj
