"""
Meridional field plot for a NON-rotating pair of polly runs: the last model in
a density sweep without a horizon, stacked over the first one with.

Each row is one run - psi on the left, the physical density rho on the right -
and the horizon, where a run has one, is drawn over both panels of its row.
The two rows share a colour scale per column, so they can be read straight
against each other: one step up in rho_0 barely moves either field, and the
thing that changes is that a horizon has appeared.

This is field_plot.py with the rotation taken out. At v_0 = 0, W and v are
zero everywhere, so those two panels would be blank and are dropped, and every
field left is even about both the axis and the equator.

rho is rho-raw.dat, the physical density. With EFFECTIVE_RHO on (source.hpp)
it is rho_eff / psi^5 once v = 0, so it dips where psi peaks.

Only the plotted window is read off disk - see the note in field_plot.py.

Run with the venv that has numpy + matplotlib:
    ./plt-env/bin/python field_plot_static.py
"""

import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke

# ----------------------------------------------------------------------------
# what to plot
# ----------------------------------------------------------------------------

# the two runs, relative to this file: the top row has no horizon, the bottom
# row does.
#
# paper-sph sweeps a uniform ball of radius 1 over rho_0 = 3.00 .. 3.10 in steps
# of 0.02. up to rho3-04 the chaser collapses to r ~ 1e-5 and gives up ("no
# usable horizon") and pollax_surf.dat is zero-filled; from rho3-06 on pollax
# finds one at r ~ 0.71, inside the matter, with a Hawking mass within 0.1% of
# M_ADM. so this pair brackets horizon formation as tightly as the sweep can
RUN_NO_AH = "../batch/paper-sph/rho3-04"
RUN_AH = "../batch/paper-sph/rho3-06"

# the e = 0.99 pancake (paper-pan) brackets it between rho5-00 (rho_0 = 5.0)
# and rho5-10 (rho_0 = 5.1), the finest step that sweep took. pass the pair on the command
# line, as below, to get field_plot_static_paper-pan_rho5-00_rho5-10.
#
# the e = 0.99 baguette (paper-bag) brackets it between rho6-80 (rho_0 = 6.8)
# and rho7-00 (rho_0 = 7.0), again the finest step that sweep took
# NB paper-pan/rho3-50 is mislabelled - its params.dat has rho_0 = 4.5

# which horizon to overlay - see field_plot.py for why pollax is the default
SURF_FILE = "pollax_surf.dat"

# half-width of the plotted window, in code units, or None to size it from the
# runs: whichever is bigger of the matter and either horizon
HALF_WIDTH = None

# padding around whatever the auto window ends up framing
WINDOW_MARGIN = 1.12

OUT_STEM = "field_plot_static"

# label each row down the left with its rho_0. unlike field_plot.py's run
# header this is what tells the two rows apart, so it is on by default; turn it
# off once a caption says which row is which
ROW_LABELS = True

# ----------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))

# two directories on the command line override the pair, no-horizon first:
#   ./plt-env/bin/python field_plot_static.py ../batch/x/run_a ../batch/x/run_b
#
# the sweep directory goes into the name too - every sweep calls its runs
# rhoX-YY, so the run names alone would collide between paper-sph and paper-pan
if len(sys.argv) == 3:
    RUN_NO_AH, RUN_AH = sys.argv[1], sys.argv[2]
    OUT_STEM = ("field_plot_static_"
                + os.path.basename(os.path.dirname(os.path.normpath(RUN_NO_AH)))
                + "_" + os.path.basename(os.path.normpath(RUN_NO_AH)) + "_"
                + os.path.basename(os.path.normpath(RUN_AH)))
elif len(sys.argv) != 1:
    sys.exit("usage: field_plot_static.py [RUN_NO_AH RUN_AH]")

run_dirs = [RUN_NO_AH, RUN_AH]
runs = [r if os.path.isabs(r) else os.path.join(HERE, r) for r in run_dirs]

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 220,
    "savefig.bbox": "tight",
    "font.size": 11,
    "axes.labelsize": 31.5,
    "axes.titlesize": 30,
    "mathtext.fontset": "cm",
})


def read_params(path):
    """params.dat is one 'key value' pair per line, as main_batch writes it."""
    p = {}
    with open(path) as f:
        for line in f:
            bits = line.split()
            if len(bits) == 2:
                p[bits[0]] = float(bits[1])
    return p


def load_window(run, name, nj, ni):
    """Read the first nj rows and ni columns of a quadrant .dat file."""
    path = os.path.join(run, name)
    arr = np.loadtxt(path, delimiter=",", max_rows=nj, usecols=range(ni))
    print(f"  {os.path.basename(run)}/{name:12s} -> {arr.shape[0]} x {arr.shape[1]}")
    return arr


def reflect(arr):
    """
    Mirror a quadrant out to the full meridional plane, evenly in both
    directions - there is no rotation here to make anything odd.
    Returned array pairs with origin="lower".
    """
    upper = np.hstack([np.fliplr(arr), arr])
    return np.vstack([np.flipud(upper), upper])


def horizon_outline(path):
    """
    (theta, f) from a surface file -> a closed (x, y) curve in the meridional
    plane, ordered pole -> equator -> pole so the loop closes without a stray
    chord. Returns (None, None) for a zero-filled file - see field_plot.py.
    """
    data = np.loadtxt(path, comments="#")
    theta, f = data[:, 0], data[:, 1]

    if f.max() <= 0. or np.isnan(data[:, 3]).all():
        return None, None

    x = f * np.sin(theta)
    y = f * np.cos(theta)

    xs = np.concatenate([x, x[::-1], -x, -x[::-1], x[:1]])
    ys = np.concatenate([y, -y[::-1], -y, y[::-1], y[:1]])
    return xs, ys


# ---- params and horizons, both runs ---------------------------------------

pars = [read_params(os.path.join(r, "params.dat")) for r in runs]

for name, par in zip(run_dirs, pars):
    # a rotating run would plot fine here and quietly lose its W and v_phi -
    # that is what field_plot.py is for
    if par["v_0"] != 0.0:
        print(f"WARNING  : {name} has v_0 = {par['v_0']:g}, so it rotates - "
              f"W and v_phi are not shown here, use field_plot.py")

# the two rows share an extent and are read off the same cells, so the grids
# have to match
nx, ny = int(pars[0]["nx"]), int(pars[0]["ny"])
for key in ("nx", "ny", "xL", "xU", "yL", "yU"):
    if pars[1][key] != pars[0][key]:
        sys.exit(f"the two runs are on different grids ({key} = "
                 f"{pars[0][key]:g} vs {pars[1][key]:g}) - pick a matching pair")
dx = (pars[0]["xU"] - pars[0]["xL"]) / nx
dy = (pars[0]["yU"] - pars[0]["yL"]) / ny

outlines = []
for name, run in zip(run_dirs, runs):
    surf_path = os.path.join(run, SURF_FILE)
    hx = hy = None
    if not os.path.exists(surf_path):
        print(f"Horizon  : no {SURF_FILE} in {name}")
    else:
        hx, hy = horizon_outline(surf_path)
        if hx is None:
            print(f"Horizon  : {name}/{SURF_FILE} is zeroed (no horizon found)")
        else:
            print(f"Horizon  : {name}/{SURF_FILE},  "
                  f"r = {np.hypot(hx, hy).min():.4f} .. {np.hypot(hx, hy).max():.4f}")
    outlines.append((hx, hy))

# the whole point of the figure is the top row without and the bottom row with
if outlines[0][0] is not None:
    print(f"WARNING  : the top run {run_dirs[0]} HAS a horizon")
if outlines[1][0] is None:
    print(f"WARNING  : the bottom run {run_dirs[1]} has NO horizon")

# ---- window ----------------------------------------------------------------

half = HALF_WIDTH
if half is None:
    reach = max(max(p["a"], p["b"]) for p in pars)
    for hx, hy in outlines:
        if hx is not None:
            reach = max(reach, np.abs(hx).max(), np.abs(hy).max())
    half = WINDOW_MARGIN * reach

ni = min(nx, int(np.ceil(half / dx)))
nj = min(ny, int(np.ceil(half / dy)))

print(f"Grid     : {nx} x {ny}, dx = {dx:.6g}")
print(f"Window   : |x| <= {ni*dx:.4g}, |y| <= {nj*dy:.4g}  "
      f"({ni} x {nj} cells of the quadrant)"
      f"{'  [auto]' if HALF_WIDTH is None else '  [HALF_WIDTH]'}")
print("Loading  :")

psis = [load_window(r, "psi.dat", nj, ni) for r in runs]
rhos = [load_window(r, "rho-raw.dat", nj, ni) for r in runs]

for name, psi, rho in zip(run_dirs, psis, rhos):
    print(f"\n  {name}")
    print(f"    psi   {psi.min():.6f} .. {psi.max():.6f}")
    print(f"    rho   {rho.min():.6g} .. {rho.max():.6g}")

# ---- panels ----------------------------------------------------------------

extent = [-ni * dx, ni * dx, -nj * dy, nj * dy]

# one scale per column across both runs, or each row would stretch its own
# range over the full map and the pair would look more different than it is
psi_lims = dict(vmin=min(p.min() for p in psis), vmax=max(p.max() for p in psis))
rho_lims = dict(vmin=0.0, vmax=max(r.max() for r in rhos))

columns = [
    (psis, r"$\psi$", "viridis", psi_lims),
    (rhos, r"$\rho_\mathrm{D}$", "magma", rho_lims),
]

fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.8), constrained_layout=True,
                         sharex=True, sharey=True)

for c, (fields, title, cmap, lims) in enumerate(columns):
    # matplotlib >= 3.9 hands back a COPY from get_cmap, so the modified map
    # has to be built and passed to imshow
    cmap = matplotlib.colormaps[cmap].with_extremes(bad="0.85")

    for r in range(2):
        ax = axes[r, c]
        im = ax.imshow(reflect(fields[r]), origin="lower", extent=extent,
                       cmap=cmap, aspect="equal", interpolation="nearest",
                       **lims)

        # white over a dark stroke, so it survives both ends of both maps
        hx, hy = outlines[r]
        if hx is not None:
            ax.plot(hx, hy, color="white", linewidth=1.6,
                    path_effects=[withStroke(linewidth=3.2, foreground="black")],
                    label="horizon")

        # ax.plot can reach past the imshow extent and matplotlib autoscales
        # rather than clipping - pin the limits, as in field_plot.py
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])

        if r == 0:
            # the default 6 pt pad was set for a 20 pt title - at 30 pt the
            # descender of psi reaches down onto the panel
            ax.set_title(title, pad=12)
        if r == 1:
            ax.set_xlabel(r"$r$")
        if c == 0:
            ax.set_ylabel(r"$z$")

        ax.axhline(0.0, color="white", linewidth=0.4, alpha=0.25)
        ax.axvline(0.0, color="white", linewidth=0.4, alpha=0.25)

    # the panels come out wide enough that matplotlib picks 0.25 steps, and
    # nine of those run together along r. the axes are shared, so this one
    # call sets every panel
    axes[0, 0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    axes[0, 0].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))

    # the scale is shared down the column, so one bar spanning both rows. a
    # colorbar keeps a fixed length:width aspect (20 by default), which would
    # stop it at about half the column - stretch it to the full height
    fig.colorbar(im, ax=axes[:, c], fraction=0.046, pad=0.02, aspect=45)

if ROW_LABELS:
    for r in range(2):
        # pinned to the z label rather than to an axes fraction, so it clears
        # the tick labels however wide they come out
        axes[r, 0].annotate(
            rf"$\rho_0 = {pars[r]['rho_0']:.2f}$",
            xy=(0.0, 0.5), xycoords=axes[r, 0].yaxis.label,
            xytext=(-12, 0), textcoords="offset points",
            rotation=90, ha="right", va="center", fontsize=27)

for r in range(2):
    if outlines[r][0] is not None:
        axes[r, 0].legend(loc="upper right", fontsize=13.5, framealpha=0.8)
        break

fig.savefig(os.path.join(HERE, OUT_STEM + ".png"))
fig.savefig(os.path.join(HERE, OUT_STEM + ".pdf"))
print(f"\nSaved {OUT_STEM}.png and {OUT_STEM}.pdf in vis/")

plt.show()
