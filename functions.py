# This Python file uses the following encoding: utf-8

from ase import Atoms
from ase.data import atomic_numbers
import pandas as pd
import numpy as np
import os


def stemcl_format_xyz(
    atoms,
    fname_out=None,
    occupancy=1.0,
    debye_waller=0.0,
    float_format="%.6f",
    **kwargs,
):
    """

    Write .xyz atomic format from ase.Atoms or .xyz file for stemcl simulation.

    Parameters
    ----------
    atoms: str or ase.Atoms
        .xyz file name or collection of Atoms.
    fname_out: None or str
        Output file name, file will be created.
        If None then will use fname appended with '_stemcl'.
        Must be defined if atoms is ase.Atoms.
    occupancy, debye_waller: float or array-like
        If float then the value will apply to all atoms.
        If array-like then must be the same length as number of atoms.
    float_format: str
        String formatting code for written floats.
    kwargs:
        Passed to pandas.read_csv.
        By default values for skiprows, delimiter, and names are defined.

    """
    columns = ("Symbol", "x", "y", "z")

    if isinstance(atoms, Atoms):
        # make sure save file name provided
        assert fname_out is not None, "fname_out must be defined if atoms is ase.Atoms."
        df = pd.DataFrame(
            data=np.column_stack((atoms.get_atomic_numbers(), atoms.get_positions())),
            columns=columns,
        )
        # force atomic numbers to be integers
        df[columns[0]] = df[columns[0]].astype(np.int8)

    elif isinstance(atoms, str):
        kwargs.setdefault("skiprows", 2)
        kwargs.setdefault("delimiter", " ")
        kwargs.setdefault("names", columns)

        # read file and get atomic data
        df = pd.read_csv(atoms, **kwargs)

        # change symbol for atomic number
        df["Symbol"] = [atomic_numbers[i] for i in df["Symbol"]]

        # sort out save file name
        if fname_out is None:
            basename, ext = os.path.splitext(atoms)
            fname_out = "".join((basename + "_stemcl", ext))
    else:
        raise TypeError("atoms must be either ase.Atoms or path to .xyz file.")

    # add occupancy column and debye-waller factor
    if isinstance(occupancy, (int, float)):
        occupancy = [occupancy] * len(df)
    if isinstance(debye_waller, (int, float)):
        debye_waller = [debye_waller] * len(df)

    df["Occupancy"] = occupancy
    df["Debye-Waller"] = debye_waller

    # write file
    with open(fname_out, "w") as f:
        f.write(f"{len(df)}\n")

        for _min, _max in zip(df[["x", "y", "z"]].min(), df[["x", "y", "z"]].max()):
            f.write(" ".join((float_format, float_format)) % (_min, _max) + "\n")

        df.to_csv(f, sep=" ", float_format=float_format, header=False, index=False)
