#!/usr/bin/env python3
"""
Read AMReX plotfiles and make Matplotlib plots.

This script supports a directory tree like

    outputs/
      FDTD_solver/
        plt00000/
        plt00001/
        ...
      pml/
        plt00000/
        plt00001/
        ...

Each pltNNNNN directory is assumed to be a standard AMReX/BoxLib plotfile
containing a Header file and one or more Level_* directories.

The script uses yt to read the plotfiles and Matplotlib to plot a selected
field/component.

Examples
--------
Plot Ez from one plotfile:

    python plot_amrex_data.py outputs/FDTD_solver/plt00010 --field Ez --axis z

Plot Ez from index 10 inside a run directory:

    python plot_amrex_data.py outputs/FDTD_solver --field Ez --axis z --index 10

Plot the last available snapshot:

    python plot_amrex_data.py outputs/FDTD_solver --field Ez --axis z --last

Make one image per plotfile:

    python plot_amrex_data.py outputs/FDTD_solver --field Ez --axis z --all --output-dir frames

Compare FDTD_solver and pml at the same index:

    python plot_amrex_data.py outputs/FDTD_solver --field Ez --axis z --index 10
    python plot_amrex_data.py outputs/pml         --field Ez --axis z --index 10

List available fields:

    python plot_amrex_data.py outputs/FDTD_solver/plt00000 --list-fields
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import yt
from matplotlib.colors import LogNorm, SymLogNorm


PLT_RE = re.compile(r"^plt(\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read AMReX plotfiles with yt and plot 2D slices using Matplotlib."
    )

    parser.add_argument(
        "path",
        type=str,
        help=(
            "Path to either a single AMReX plotfile, e.g. outputs/FDTD_solver/plt00010, "
            "or a run directory containing pltNNNNN folders, e.g. outputs/FDTD_solver."
        ),
    )

    parser.add_argument(
        "--field",
        type=str,
        default="rho",
        help="Field name to plot, e.g. rho, Ex, Ey, Ez, Hx, Hy, Hz.",
    )

    parser.add_argument(
        "--axis",
        type=str,
        choices=("x", "y", "z"),
        default="z",
        help="Slice normal direction. Default: z.",
    )

    parser.add_argument(
        "--coord",
        type=float,
        default=None,
        help="Slice coordinate in physical units. Default: domain center along selected axis.",
    )

    parser.add_argument(
        "--level",
        type=int,
        default=None,
        help="AMR level used for the covering grid. Default: finest available level.",
    )

    parser.add_argument(
        "--index",
        type=int,
        default=None,
        help="Snapshot index to load when path is a run directory. Example: --index 10 loads plt00010.",
    )

    parser.add_argument(
        "--last",
        action="store_true",
        help="When path is a run directory, load the largest available pltNNNNN index.",
    )

    parser.add_argument(
        "--all",
        action="store_true",
        help="When path is a run directory, generate one image for every pltNNNNN plotfile.",
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output image filename for a single plot. Default: <run>_<plt>_<field>_<axis>slice.png.",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default="frames",
        help="Output directory when using --all. Default: frames.",
    )

    parser.add_argument(
        "--cmap",
        type=str,
        default="viridis",
        help="Matplotlib colormap. Default: viridis.",
    )

    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Lower colorbar limit.",
    )

    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Upper colorbar limit.",
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help="Use logarithmic color normalization. Field values must be positive.",
    )

    parser.add_argument(
        "--symlog",
        action="store_true",
        help="Use symmetric logarithmic normalization for signed fields.",
    )

    parser.add_argument(
        "--linthresh",
        type=float,
        default=None,
        help="Linear threshold for --symlog. Default: 1e-3 times max(abs(field)).",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Output image resolution. Default: 200.",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the plot interactively. Usually only useful for a single plot.",
    )

    parser.add_argument(
        "--list-fields",
        action="store_true",
        help="List available fields and exit.",
    )

    return parser.parse_args()


def is_plotfile(path: Path) -> bool:
    return path.is_dir() and (path / "Header").is_file()


def plotfile_index(path: Path) -> Optional[int]:
    m = PLT_RE.match(path.name)
    if m is None:
        return None
    return int(m.group(1))


def find_plotfiles(path: Path) -> list[Path]:
    """Return sorted pltNNNNN plotfiles inside path."""
    if is_plotfile(path):
        return [path]

    if not path.is_dir():
        raise FileNotFoundError(f"Not a directory: {path}")

    plotfiles = []
    for child in path.iterdir():
        if is_plotfile(child) and plotfile_index(child) is not None:
            plotfiles.append(child)

    plotfiles.sort(key=lambda p: plotfile_index(p))

    if not plotfiles:
        raise FileNotFoundError(f"No AMReX plotfiles of the form pltNNNNN found in {path}")

    return plotfiles


def select_plotfiles(path: Path, args: argparse.Namespace) -> list[Path]:
    plotfiles = find_plotfiles(path)

    if is_plotfile(path):
        return plotfiles

    selectors = sum([args.index is not None, args.last, args.all])
    if selectors > 1:
        raise ValueError("Use only one of --index, --last, or --all.")

    if args.all:
        return plotfiles

    if args.last:
        return [plotfiles[-1]]

    if args.index is not None:
        target = f"plt{args.index:05d}"
        for p in plotfiles:
            if p.name == target:
                return [p]
        available = ", ".join(p.name for p in plotfiles)
        raise FileNotFoundError(f"Could not find {target} in {path}. Available: {available}")

    # Default for a run directory: use the last snapshot.
    return [plotfiles[-1]]


def find_field(ds: yt.data_objects.static_output.Dataset, name: str):
    """Return the yt field tuple matching a user field name."""
    all_fields = list(ds.field_list) + list(ds.derived_field_list)

    # First try exact field-name match, ignoring field type.
    matches = [f for f in all_fields if f[-1] == name]
    if matches:
        return matches[0]

    # Then try direct tuple input style, in case the user passed boxlib,rho.
    if "," in name:
        parts = tuple(s.strip() for s in name.split(","))
        if parts in all_fields:
            return parts

    available = sorted({f[-1] for f in all_fields})
    raise ValueError(
        f"Could not find field '{name}'. Available field names include:\n"
        + "\n".join(available)
    )


def get_domain_info(ds, level: Optional[int]):
    """Return AMR level, left edge, right edge, and covering-grid dimensions."""
    if level is None:
        level = int(ds.max_level)

    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge

    # ds.domain_dimensions is the base-level resolution.
    # For standard AMReX refinement by 2 per level, this gives the finest uniform grid.
    dims = ds.domain_dimensions * (2**level)

    return level, left_edge, right_edge, dims


def extract_slice(
    ds,
    field,
    axis: str,
    coord: Optional[float],
    level: Optional[int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, str, str, float]:
    """
    Extract a 2D slice from a uniform covering grid.

    Returns
    -------
    x2d, y2d, data2d, xlabel, ylabel, slice_coord
    """
    level, left_edge, right_edge, dims = get_domain_info(ds, level)

    cg = ds.covering_grid(
        level=level,
        left_edge=left_edge,
        dims=dims,
        num_ghost_zones=0,
    )

    data = np.asarray(cg[field].to_ndarray())

    left = np.asarray(left_edge.to_value())
    right = np.asarray(right_edge.to_value())
    dims_np = np.asarray(dims, dtype=int)

    x = np.linspace(left[0], right[0], dims_np[0], endpoint=False)
    y = np.linspace(left[1], right[1], dims_np[1], endpoint=False)
    z = np.linspace(left[2], right[2], dims_np[2], endpoint=False)

    # Move from cell-left coordinates to approximate cell-center coordinates.
    dx = (right - left) / dims_np
    x += 0.5 * dx[0]
    y += 0.5 * dx[1]
    z += 0.5 * dx[2]

    if axis == "x":
        if coord is None:
            coord = 0.5 * (left[0] + right[0])
        idx = int(np.argmin(np.abs(x - coord)))
        data2d = data[idx, :, :].T
        xx, yy = np.meshgrid(y, z, indexing="xy")
        xlabel, ylabel = "y", "z"
        slice_coord = x[idx]

    elif axis == "y":
        if coord is None:
            coord = 0.5 * (left[1] + right[1])
        idx = int(np.argmin(np.abs(y - coord)))
        data2d = data[:, idx, :].T
        xx, yy = np.meshgrid(x, z, indexing="xy")
        xlabel, ylabel = "x", "z"
        slice_coord = y[idx]

    elif axis == "z":
        if coord is None:
            coord = 0.5 * (left[2] + right[2])
        idx = int(np.argmin(np.abs(z - coord)))
        data2d = data[:, :, idx].T
        xx, yy = np.meshgrid(x, y, indexing="xy")
        xlabel, ylabel = "x", "y"
        slice_coord = z[idx]

    else:
        raise ValueError(f"Invalid axis: {axis}")

    return xx, yy, data2d, xlabel, ylabel, float(slice_coord)


def make_norm(data: np.ndarray, args: argparse.Namespace):
    """Construct an optional Matplotlib normalization."""
    if args.log and args.symlog:
        raise ValueError("Use either --log or --symlog, not both.")

    vmin = args.vmin
    vmax = args.vmax

    if args.log:
        positive = data[data > 0]
        if positive.size == 0:
            raise ValueError("--log requested, but the selected field has no positive values.")
        if vmin is None:
            vmin = float(np.nanmin(positive))
        if vmax is None:
            vmax = float(np.nanmax(positive))
        return LogNorm(vmin=vmin, vmax=vmax)

    if args.symlog:
        absmax = float(np.nanmax(np.abs(data)))
        linthresh = args.linthresh if args.linthresh is not None else max(1e-30, 1e-3 * absmax)
        if vmin is None:
            vmin = -absmax
        if vmax is None:
            vmax = absmax
        return SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax)

    return None


def output_name_for_plotfile(plotfile: Path, args: argparse.Namespace) -> Path:
    if args.output is not None:
        return Path(args.output)

    run_name = plotfile.parent.name
    filename = f"{run_name}_{plotfile.name}_{args.field}_{args.axis}slice.png"

    if args.all:
        outdir = Path(args.output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        return outdir / filename

    return Path(filename)


def plot_slice(
    xx,
    yy,
    data2d,
    field_name: str,
    axis: str,
    slice_coord: float,
    xlabel: str,
    ylabel: str,
    title_prefix: str,
    args,
):
    fig, ax = plt.subplots(figsize=(7.0, 5.5), constrained_layout=True)

    norm = make_norm(data2d, args)

    im = ax.pcolormesh(
        xx,
        yy,
        data2d,
        shading="auto",
        cmap=args.cmap,
        norm=norm,
        vmin=None if norm is not None else args.vmin,
        vmax=None if norm is not None else args.vmax,
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(field_name)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{title_prefix}: {field_name}, {axis} = {slice_coord:.6g}")

    return fig, ax


def process_one_plotfile(plotfile: Path, args: argparse.Namespace) -> None:
    print(f"Reading {plotfile}")
    ds = yt.load(str(plotfile))

    if args.list_fields:
        print("Available fields:")
        for f in sorted(ds.field_list):
            print(f"  {f}")
        return

    field = find_field(ds, args.field)

    xx, yy, data2d, xlabel, ylabel, slice_coord = extract_slice(
        ds=ds,
        field=field,
        axis=args.axis,
        coord=args.coord,
        level=args.level,
    )

    fig, _ = plot_slice(
        xx=xx,
        yy=yy,
        data2d=data2d,
        field_name=args.field,
        axis=args.axis,
        slice_coord=slice_coord,
        xlabel=xlabel,
        ylabel=ylabel,
        title_prefix=f"{plotfile.parent.name}/{plotfile.name}",
        args=args,
    )

    output = output_name_for_plotfile(plotfile, args)
    fig.savefig(output, dpi=args.dpi)
    print(f"Saved {output}")

    if args.show:
        plt.show()

    plt.close(fig)


def main() -> None:
    args = parse_args()

    path = Path(args.path)
    plotfiles = select_plotfiles(path, args)

    for plotfile in plotfiles:
        process_one_plotfile(plotfile, args)


if __name__ == "__main__":
    main()
