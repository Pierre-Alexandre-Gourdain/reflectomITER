#!/usr/bin/env python3

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import yt


# ============================================================
# Physical constants
# ============================================================

qe   = 1.602176634e-19
me0  = 9.1093837015e-31
eps0 = 8.8541878128e-12
c0   = 299792458.0


# ============================================================
# Argument parsing
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
Generate PNG frames from an AMReX / yt plotfile sequence listed in
FDTD_solver.visit. The output PNGs are suitable for LaTeX inclusion.

Example:

    python make_reflectometry_frames.py \\
        --directory data_40 \\
        --frequency 40e9 \\
        --mode O

Modes:
    O     : O-mode cutoff
    X_R   : right-hand X-mode cutoff
    X_L   : left-hand X-mode cutoff
    UH    : upper-hybrid resonance
""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--directory",
        "-d",
        type=str,
        required=False,
        help="Directory containing FDTD_solver.visit, for example data_40.",
    )

    parser.add_argument(
        "--frequency",
        "-f",
        type=float,
        required=False,
        help="Probe frequency in Hz, for example 40e9.",
    )

    parser.add_argument(
        "--mode",
        "-m",
        type=str,
        required=False,
        choices=["O", "X_R", "X_L", "UH"],
        help="Cutoff mode: O, X_R, X_L, or UH.",
    )

    parser.add_argument(
        "--axis",
        type=str,
        default="y",
        choices=["x", "y", "z"],
        help="Slice axis. Default: y.",
    )

    parser.add_argument(
        "--coord",
        type=float,
        default=None,
        help="Slice coordinate. Default: domain center.",
    )

    parser.add_argument(
        "--level",
        type=int,
        default=0,
        help="AMR level to plot. Default: 0.",
    )

    parser.add_argument(
        "--read-stride",
        type=int,
        default=1,
        help="Stride used while reading/slicing data. Default: 1.",
    )

    parser.add_argument(
        "--plot-stride",
        type=int,
        default=1,
        help="Stride used while plotting data. Default: 1.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=250,
        help="PNG output resolution. Default: 250.",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="movie_frames",
        help="Output folder name. Default: movie_frames.",
    )

    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional maximum number of frames, useful for testing.",
    )

    args = parser.parse_args()

    if args.directory is None or args.frequency is None or args.mode is None:
        parser.print_help()
        sys.exit(0)

    return args


# ============================================================
# Physics helper
# ============================================================

def compute_cutoff_frequency(ne, Te, Babs, mode="O"):
    """
    Return selected angular frequency.

    ne   : electron density in m^-3
    Te   : electron temperature in eV
    Babs : magnetic-field magnitude in Tesla

    mode:
        O   : omega_pe
        X_R : right-hand X-mode cutoff
        X_L : left-hand X-mode cutoff
        UH  : upper-hybrid resonance
    """

    ne_safe = np.maximum(ne, 0.0)
    Te_safe = np.maximum(Te, 0.0)
    B_safe  = np.maximum(Babs, 0.0)

    me_rel = me0 * np.sqrt(
        1.0 + 5.0 * qe * Te_safe / (me0 * c0**2)
    )

    omega_pe = np.sqrt(qe**2 * ne_safe / (eps0 * me_rel))
    Omega_e = qe * B_safe / me_rel

    mode = mode.upper()

    if mode == "O":
        return omega_pe

    if mode == "X_R":
        return 0.5 * Omega_e + np.sqrt(0.25 * Omega_e**2 + omega_pe**2)

    if mode == "X_L":
        return -0.5 * Omega_e + np.sqrt(0.25 * Omega_e**2 + omega_pe**2)

    if mode == "UH":
        return np.sqrt(Omega_e**2 + omega_pe**2)

    raise ValueError("mode must be O, X_R, X_L, or UH")


# ============================================================
# Geometry helper
# ============================================================

def extract_2d_slice(data3d, x, y, z, axis="z", coord=None, stride=1):
    """
    Extract a coordinate-aligned 2D plane from a 3D array.

    axis = z gives x-y plane at fixed z.
    axis = y gives x-z plane at fixed y.
    axis = x gives y-z plane at fixed x.

    data3d is assumed to have shape (nx, ny, nz).
    """

    if axis == "x":
        if coord is None:
            coord = 0.5 * (x.min() + x.max())

        idx = int(np.argmin(np.abs(x - coord)))

        yy = y[::stride]
        zz = z[::stride]

        data2d = data3d[idx, ::stride, ::stride].T
        X, Y = np.meshgrid(yy, zz, indexing="xy")

        xlabel = "y (m)"
        ylabel = "z (m)"
        actual_coord = x[idx]

    elif axis == "y":
        if coord is None:
            coord = 0.5 * (y.min() + y.max())

        idx = int(np.argmin(np.abs(y - coord)))

        xx = x[::stride]
        zz = z[::stride]

        data2d = data3d[::stride, idx, ::stride].T
        X, Y = np.meshgrid(xx, zz, indexing="xy")

        xlabel = "r (m)"
        ylabel = "z (m)"
        actual_coord = y[idx]

    elif axis == "z":
        if coord is None:
            coord = 0.5 * (z.min() + z.max())

        idx = int(np.argmin(np.abs(z - coord)))

        xx = x[::stride]
        yy = y[::stride]

        data2d = data3d[::stride, ::stride, idx].T
        X, Y = np.meshgrid(xx, yy, indexing="xy")

        xlabel = "r (m)"
        ylabel = "y (m)"
        actual_coord = z[idx]

    else:
        raise ValueError("axis must be x, y, or z")

    return X, Y, data2d, xlabel, ylabel, actual_coord


# ============================================================
# Plotting helper
# ============================================================

def add_selected_cutoff_contour(
    ax,
    X,
    Y,
    omega_cutoff,
    omega_probe,
    color="white",
):
    quantity = omega_cutoff - omega_probe

    qmin = np.nanmin(quantity)
    qmax = np.nanmax(quantity)

    if qmin <= 0.0 <= qmax:
        ax.contour(
            X,
            Y,
            quantity,
            levels=[0.0],
            colors=color,
            linewidths=2.5,
            linestyles=":",
        )

def figure_size_from_domain(
    X,
    Y,
    nrows=4,
    ncols=2,
    panel_width=4.2,
    colorbar_extra=0.18,
    title_extra=0.8,
):
    """
    Choose figsize so each subplot has approximately the physical
    aspect ratio of the plotted domain, while keeping columns close.
    """

    x_extent = np.nanmax(X) - np.nanmin(X)
    y_extent = np.nanmax(Y) - np.nanmin(Y)

    if x_extent <= 0 or y_extent <= 0:
        return (12, 9)

    data_aspect = y_extent / x_extent

    fig_width = ncols * panel_width * (1.0 + colorbar_extra)
    fig_height = nrows * panel_width * data_aspect + title_extra

    return (fig_width, fig_height)
    
def make_frame(
    plotfile,
    frame_path,
    frequency_hz,
    mode,
    axis="y",
    coord=None,
    level=0,
    read_stride=1,
    plot_stride=1,
    dpi=250,
):
    field_type = "boxlib"
    omega_probe = 2.0 * np.pi * frequency_hz

    def safe_array(cg, field_name):
        return np.asarray(cg[(field_type, field_name)].to_ndarray())

    ds = yt.load(str(plotfile))

    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    dims = ds.domain_dimensions * (2**level)

    cg = ds.covering_grid(
        level=level,
        left_edge=left_edge,
        dims=dims,
        num_ghost_zones=0,
    )

    left = np.asarray(left_edge.to_value())
    right = np.asarray(right_edge.to_value())
    dims_np = np.asarray(dims, dtype=int)

    x = np.linspace(left[0], right[0], dims_np[0], endpoint=False)
    y = np.linspace(left[1], right[1], dims_np[1], endpoint=False)
    z = np.linspace(left[2], right[2], dims_np[2], endpoint=False)

    dx = (right - left) / dims_np

    x += 0.5 * dx[0]
    y += 0.5 * dx[1]
    z += 0.5 * dx[2]

    X, Y, _, xlabel, ylabel, actual_coord = extract_2d_slice(
        safe_array(cg, "Ex"),
        x,
        y,
        z,
        axis=axis,
        coord=coord,
        stride=read_stride,
    )

    Ex_2d = extract_2d_slice(safe_array(cg, "Ex"), x, y, z, axis, coord, read_stride)[2]
    Ey_2d = extract_2d_slice(safe_array(cg, "Ey"), x, y, z, axis, coord, read_stride)[2]
    Ez_2d = extract_2d_slice(safe_array(cg, "Ez"), x, y, z, axis, coord, read_stride)[2]

    Hx_2d = extract_2d_slice(safe_array(cg, "Hx"), x, y, z, axis, coord, read_stride)[2]
    Hy_2d = extract_2d_slice(safe_array(cg, "Hy"), x, y, z, axis, coord, read_stride)[2]
    Hz_2d = extract_2d_slice(safe_array(cg, "Hz"), x, y, z, axis, coord, read_stride)[2]

    Bx_2d = extract_2d_slice(safe_array(cg, "Bx"), x, y, z, axis, coord, read_stride)[2]
    By_2d = extract_2d_slice(safe_array(cg, "By"), x, y, z, axis, coord, read_stride)[2]
    Bz_2d = extract_2d_slice(safe_array(cg, "Bz"), x, y, z, axis, coord, read_stride)[2]

    Jx_2d = extract_2d_slice(safe_array(cg, "Jx"), x, y, z, axis, coord, read_stride)[2]
    Jy_2d = extract_2d_slice(safe_array(cg, "Jy"), x, y, z, axis, coord, read_stride)[2]
    Jz_2d = extract_2d_slice(safe_array(cg, "Jz"), x, y, z, axis, coord, read_stride)[2]

    ne_2d = extract_2d_slice(safe_array(cg, "ne"), x, y, z, axis, coord, read_stride)[2]
    Te_2d = extract_2d_slice(safe_array(cg, "Te"), x, y, z, axis, coord, read_stride)[2]

    rho_2d = extract_2d_slice(safe_array(cg, "rho"), x, y, z, axis, coord, read_stride)[2]
    Eerr_2d = extract_2d_slice(safe_array(cg, "Eerr"), x, y, z, axis, coord, read_stride)[2]

    Eerr_2d = np.log10(10**Eerr_2d + 1e-18)

    E_norm = np.sqrt(Ex_2d**2 + Ey_2d**2 + Ez_2d**2)
    H_norm = np.sqrt(Hx_2d**2 + Hy_2d**2 + Hz_2d**2)
    B_norm = np.sqrt(Bx_2d**2 + By_2d**2 + Bz_2d**2)
    J_norm = np.sqrt(Jx_2d**2 + Jy_2d**2 + Jz_2d**2)

    omega_cutoff = compute_cutoff_frequency(
        ne_2d,
        Te_2d,
        B_norm,
        mode=mode,
    )

    plot_data = [
        [E_norm, ne_2d],
        [H_norm, Te_2d],
        [J_norm, B_norm],
        [rho_2d, Eerr_2d],
    ]

    plot_cmaps = [
        ["viridis", "jet"],
        ["viridis", "inferno"],
        ["viridis", "nipy_spectral"],
        ["RdBu", "nipy_spectral"],
    ]

    plot_center_zero = [
        [False, False],
        [False, False],
        [False, False],
        [True,  False],
    ]

    plot_titles = [
        [r"$|\mathbf{E}|$", r"$n_e$"],
        [r"$|\mathbf{H}|$", r"$T_e$"],
        [r"$|\mathbf{J}|$", r"$|\mathbf{B}|$"],
        [r"$\rho$", r"$\log_{10}|\nabla\cdot\mathbf{E}|$ error"],
    ]

    plot_scales = [
        [1e-3, 1e19],
        [1e-6, 1e3],
        [1e-3, 1],
        [1e-15, 1],
    ]

    plot_units = [
        [r"$mV/m$", r"$\times 10^{19}\,m^{-3}$"],
        [r"$\mu A/m$", r"$keV$"],
        [r"$mA/m^2$", r"$T$"],
        [r"$fC/m^3$", r""],
    ]

    figsize = figure_size_from_domain(
        X,
        Y,
        nrows=4,
        ncols=2,
        panel_width=4.2,
        colorbar_extra=0.18,
    )    
    
    fig, axes = plt.subplots(
        4,
        2,
        figsize=figsize,
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )

    fig.set_constrained_layout_pads(
        w_pad=0.1,
        h_pad=0.1,
        wspace=0.05,
        hspace=0.05,
    )

    for i in range(4):
        for j in range(2):
            ax = axes[i, j]

            data_plot = plot_data[i][j] / plot_scales[i][j]

            if plot_center_zero[i][j]:
                v = np.nanmax(np.abs(data_plot))

                if v == 0 or not np.isfinite(v):
                    v = 1.0

                norm = TwoSlopeNorm(
                    vmin=-v,
                    vcenter=0.0,
                    vmax=v,
                )

                im = ax.pcolormesh(
                    X[::plot_stride],
                    Y[::plot_stride],
                    data_plot[::plot_stride],
                    shading="auto",
                    cmap=plot_cmaps[i][j],
                    norm=norm,
                    rasterized=True,
                )

            else:
                vmin = np.nanmin(data_plot)
                vmax = np.nanmax(data_plot)

                if not np.isfinite(vmin):
                    vmin = 0.0

                if not np.isfinite(vmax):
                    vmax = 1.0

                if vmin == vmax:
                    vmax = vmin + 1.0

                im = ax.pcolormesh(
                    X[::plot_stride],
                    Y[::plot_stride],
                    data_plot[::plot_stride],
                    shading="auto",
                    cmap=plot_cmaps[i][j],
                    vmin=vmin,
                    vmax=vmax,
                    rasterized=True,
                )

            fig.colorbar(
                im,
                ax=ax,
                shrink=0.75,
                fraction=0.046,
                pad=0.04,
                label=plot_units[i][j],
            )

            if i <= 2:
                add_selected_cutoff_contour(
                    ax,
                    X,
                    Y,
                    omega_cutoff,
                    omega_probe,
                    color="white",
                )

            ax.set_aspect("equal", adjustable="box")
            ax.set_title(plot_titles[i][j])

    for ax in axes[-1, :]:
        ax.set_xlabel(xlabel)

    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)

    fig.suptitle(
        fr"{plotfile.name}: ${mode}$-mode cutoff "
        fr"{frequency_hz / 1e9:.4g} GHz at {axis} = {actual_coord:.4g}",
        fontsize=36,
    )

    fig.savefig(
        frame_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)

    del ds, cg
    gc.collect()


# ============================================================
# Visit-file reader
# ============================================================

def read_visit_file(data_dir):
    visit_file = data_dir / "FDTD_solver.visit"

    if not visit_file.exists():
        raise FileNotFoundError(f"Could not find visit file: {visit_file}")

    plotfiles = []
    seen = set()

    with open(visit_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            p = data_dir / line

            # VisIt/AMReX files sometimes list:
            #     plt00000/Header
            # but yt.load wants:
            #     plt00000
            if p.name == "Header":
                p = p.parent

            if p.exists():
                p = p.resolve()

                if p not in seen:
                    plotfiles.append(p)
                    seen.add(p)
            else:
                print(f"Warning: missing plotfile listed in visit file: {p}")

    return plotfiles

# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    data_dir = Path(args.directory).expanduser().resolve()

    if not data_dir.exists():
        print(f"Error: directory does not exist: {data_dir}")
        sys.exit(1)

    plotfiles = read_visit_file(data_dir)

    if args.max_frames is not None:
        plotfiles = plotfiles[:args.max_frames]

    if len(plotfiles) == 0:
        print("Error: no valid plotfiles found.")
        sys.exit(1)

    frame_dir = data_dir / args.output
    frame_dir.mkdir(parents=True, exist_ok=True)

    print("============================================================")
    print("Reflectometry PNG frame generator")
    print("============================================================")
    print(f"Data directory : {data_dir}")
    print(f"Visit file     : {data_dir / 'FDTD_solver.visit'}")
    print(f"Output folder  : {frame_dir}")
    print(f"Frequency      : {args.frequency:.6e} Hz")
    print(f"Mode           : {args.mode}")
    print(f"Slice axis     : {args.axis}")
    print(f"Slice coord    : {args.coord}")
    print(f"AMR level      : {args.level}")
    print(f"Frames         : {len(plotfiles)}")
    print("============================================================")

    for iframe, plotfile in enumerate(plotfiles):
        frame_path = frame_dir / f"frame_{iframe:05d}.png"

        print(f"[{iframe + 1}/{len(plotfiles)}] {plotfile.name} -> {frame_path.name}")

        make_frame(
            plotfile=plotfile,
            frame_path=frame_path,
            frequency_hz=args.frequency,
            mode=args.mode,
            axis=args.axis,
            coord=args.coord,
            level=args.level,
            read_stride=args.read_stride,
            plot_stride=args.plot_stride,
            dpi=args.dpi,
        )

    print("============================================================")
    print("Done.")
    print(f"Saved PNG frames in: {frame_dir}")
    print("============================================================")


if __name__ == "__main__":
    main()