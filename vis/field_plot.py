"""
Four-panel meridional field plot for one polly run, with a horizon contour.

polly stores only the x >= 0, z >= 0 quadrant, so every field is reflected out
to the full meridional plane here. The reflection is NOT the same for every
field: psi, rho and v are even about both the axis and the equator, but W is
ODD about the rotation axis - the x < 0 half of the picture is the meridional
plane at phi + pi, where the azimuthal shift reverses sign. That is exactly the
parity Grid::fill_all_ghosts_robin imposes on the ghost cells, so plotting W
evenly would contradict the boundary condition the solve was run with.

v_phys is the PHYSICAL rotational velocity of the dust,

    v_phys = | v_raw * x * psi^2 |

which is the quantity Source::sanity_check compares against 1: at 1 the matter
is moving at the speed of light and the raw density changes sign, so the run is
unphysical from there on. It is not a saved field - it is rebuilt here from
v-raw.dat and psi.dat.

It is plotted here on its own scale, as a velocity. For the threshold-style
version - vmax pinned at 1 so the fill reads as a fraction of the light speed
limit, white reserved for "at or past it" - see field_plot_lum.py.

Only the plotted window is read off disk. The .dat files are a 4096^2 CSV each
(~300-400 MB), so reading all of every one of them to show the middle 8% would
cost minutes and gigabytes; max_rows/usecols cut that to a fraction of a second.

Run with the venv that has numpy + matplotlib:
    ./plt-env/bin/python field_plot.py
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

# the run directory, relative to this file - resolved below so the script can
# be run from anywhere, not just from vis/
RUN_DIR = "../batch/twist/test-anti"
# RUN_DIR = "../batch/paper-diffrotpan-HV/v100"
# RUN_DIR = "../batch/paper-diffrotsph/v18"

# NOTE bag-v2 (rho_0 = 5.1, v_0 = 2) has NO horizon - its output.txt says
# "no horizon found, AH rows are zeroed", and chased/pollax_surf.dat are
# zero-filled placeholders rather than surfaces. Point RUN_DIR at a collapsed
# run to see the horizon overlay, e.g. the same baguette at rho_0 = 7, v_0 = 10:
# RUN_DIR = "../batch/paper-bag-v/rho-7_00-v10"

# which horizon to overlay. the three finders in cauchy-tools each write the
# same "# theta  f  psi  div_s  Ricci" format:
#   pollax_surf.dat   the polygrid relaxer   <- the default
#   chased_surf.dat   the chaser
#   shot_surf.dat     the shooter
#
# pollax and the chaser are independent relaxations of the same condition and
# agree to six digits where both converge - on naka999/e999 both land on
# M = 19.9082, r = 9.4895 .. 10.8631 at res ~ 1e-11. the shooter is the one to
# distrust: on that same run it returns M = 22.83 over r = 3.01 .. 8.16, a
# different surface entirely, so it should not be the thing a plot defaults to.
# pollax gets the default over the chaser only because it converges in far
# fewer iterations for the same answer
SURF_FILE = "pollax_surf.dat"

# is W ODD about the equator, as Grid::fill_all_ghosts_antisymmW imposes? then
# the hemispheres counter-rotate, the net angular momentum vanishes, and the far
# field leads with the l = 2 term ~ sin(th)cos(th)/R^3 rather than l = 1
# ~ sin(th)/R^2.
#
# polly stores only the y >= 0 quadrant, so nothing in the saved data records
# this - it is the reflection below that has to carry the right sign into the
# lower half plane, for W and for v_phi, which reverses with it.
ANTISYMM_W = False

# half-width of the plotted window, in code units, or None to size it from the
# run itself. None is almost always what you want: it frames whichever is
# bigger, the horizon or the matter, so a run whose horizon sits out at r ~ 11
# (the e999 baguette) needs no retuning of a number in here. a hardcoded value
# SMALLER than the horizon is what used to read as the plot zooming out into
# white space - see the set_xlim note in the panel loop for why.
HALF_WIDTH = None

# padding around whatever the auto window ends up framing
WINDOW_MARGIN = 1.12

OUT_STEM = "field_plot"

# run header across the top: the directory name, rho_0, v_0, a, b, e, M_ADM and
# which horizon file was drawn. useful when flicking through runs, redundant
# once the figure has a caption carrying the same numbers, so it is off by
# default and the panels get the space back
SHOW_HEADER = False

# ----------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))

# a directory on the command line overrides RUN_DIR, so one run can be checked
# without editing the file:  ./plt-env/bin/python field_plot.py ../batch/.../run
if len(sys.argv) > 1:
    RUN_DIR = sys.argv[1]
    OUT_STEM = "field_plot_" + os.path.basename(os.path.normpath(RUN_DIR))

run = RUN_DIR if os.path.isabs(RUN_DIR) else os.path.join(HERE, RUN_DIR)

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 220,
    "savefig.bbox": "tight",
    "font.size": 11,
    "axes.labelsize": 21,
    "axes.titlesize": 20,
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


def load_window(name, nj, ni):
    """Read the first nj rows and ni columns of a quadrant .dat file."""
    path = os.path.join(run, name)
    arr = np.loadtxt(path, delimiter=",", max_rows=nj, usecols=range(ni))
    print(f"  {name:16s} -> {arr.shape[0]} x {arr.shape[1]}  (of the full grid)")
    return arr


def reflect(arr, parity_x=+1, parity_y=+1):
    """
    Mirror a quadrant out to the full meridional plane.

    parity_x = -1 for a field that is odd about the rotation axis (W),
    parity_y = -1 for one odd about the equator (nothing here is).
    Returned array pairs with origin="lower".
    """
    upper = np.hstack([parity_x * np.fliplr(arr), arr])
    lower = parity_y * np.flipud(upper)
    return np.vstack([lower, upper])


def horizon_outline(path):
    """
    (theta, f) from a surface file -> a closed (x, y) curve in the meridional
    plane. theta runs from the pole to the equator over the first quadrant, so
    x = f sin(theta) is the cylindrical radius and y = f cos(theta) the height.
    The four segments are ordered pole -> equator -> pole so the loop closes
    without a stray chord across the middle.
    """
    data = np.loadtxt(path, comments="#")
    theta, f = data[:, 0], data[:, 1]

    # cauchy-tools still writes the file when the finder fails - every f is 0
    # and div_s/Ricci are nan ("no horizon found, AH rows are zeroed"). drawing
    # that gives a dot at the origin, which reads as a real surface
    if f.max() <= 0. or np.isnan(data[:, 3]).all():
        return None, None

    x = f * np.sin(theta)
    y = f * np.cos(theta)

    xs = np.concatenate([x, x[::-1], -x, -x[::-1], x[:1]])
    ys = np.concatenate([y, -y[::-1], -y, y[::-1], y[:1]])
    return xs, ys


# ---- grid geometry ---------------------------------------------------------

par = read_params(os.path.join(run, "params.dat"))
nx, ny = int(par["nx"]), int(par["ny"])
dx = (par["xU"] - par["xL"]) / nx
dy = (par["yU"] - par["yL"]) / ny

# ---- the horizon, read before the window is sized -------------------------
#
# a run with no horizon, or one where the finder failed, still plots its fields.
# this has to come before the window is chosen, because the outline is part of
# what the window is sized to
surf_path = os.path.join(run, SURF_FILE)
hx = hy = None
if not os.path.exists(surf_path):
    print(f"Horizon  : no {SURF_FILE} in this run - plotting without one")
else:
    hx, hy = horizon_outline(surf_path)
    if hx is None:
        print(f"Horizon  : {SURF_FILE} is zeroed (no horizon found) - none drawn")
    else:
        print(f"Horizon  : {SURF_FILE},  "
              f"r = {np.hypot(hx, hy).min():.4f} .. {np.hypot(hx, hy).max():.4f}")

# ---- window ----------------------------------------------------------------
#
# frame the horizon and the matter together. the matter semi-axes alone are the
# fallback for a run that never collapsed, and the domain is the hard ceiling -
# asking for more than was saved would just pad the picture with nothing
half = HALF_WIDTH
if half is None:
    reach = max(par["a"], par["b"])
    if hx is not None:
        reach = max(reach, np.abs(hx).max(), np.abs(hy).max())
    half = WINDOW_MARGIN * reach

# how much of the quadrant the window actually needs
ni = min(nx, int(np.ceil(half / dx)))
nj = min(ny, int(np.ceil(half / dy)))

print(f"Run      : {RUN_DIR}")
print(f"Grid     : {nx} x {ny}, dx = {dx:.6g}")
print(f"Window   : |x| <= {ni*dx:.4g}, |y| <= {nj*dy:.4g}  "
      f"({ni} x {nj} cells of the quadrant)"
      f"{'  [auto]' if HALF_WIDTH is None else '  [HALF_WIDTH]'}")
print("Loading  :")

psi = load_window("psi.dat", nj, ni)
W = load_window("W.dat", nj, ni)
rho_raw = load_window("rho-raw.dat", nj, ni)
rho_eff = load_window("rho-eff.dat", nj, ni)
v_raw = load_window("v-raw.dat", nj, ni)

# cell-centred cylindrical radius, as Field::get_x computes it
x_centres = (np.arange(ni) + 0.5) * dx + par["xL"]
v_phys = np.abs(v_raw * x_centres[None, :] * psi**2)

# zeroed outside the matter ellipse, which is also exactly the set
# Source::sanity_check counts over. the vacuum carries the same factor on zero
# density, and since v keeps its v_0 floor everywhere while x keeps growing,
# the empty grid runs up to v_phys ~ 25 at the window edge - a number that
# describes nothing, but which would own the colour scale and bury the range
# that actually matters
# the mask is taken from rho_eff, NOT rho_raw, although the two share a support.
# rho_raw carries a factor 1 - v_phys^2 and goes NEGATIVE once rotation drives
# the matter superluminal, so "rho_raw > 0" would quietly drop exactly the cells
# the luminal panel exists to show. rho_eff cannot change sign
matter = rho_eff > 0.0
v_phys = np.where(matter, v_phys, 0.0)

v_max = v_phys.max()
n_bad = int((v_phys > 1.0).sum())
n_neg = int((rho_raw < 0.0).sum())
print(f"\n  matter cells in window = {int(matter.sum())}")
if n_neg:
    print(f"  *** rho (physical) is NEGATIVE in {n_neg} cells - superluminal matter,"
          f" min = {rho_raw.min():.6g}.  the panel clips these to 0 ***")
print(f"  v_phys max = {v_max:.6f}"
      f"   cells above 1 = {n_bad}"
      f"{'   *** SUPERLUMINAL ***' if n_bad else ''}")

# ---- panels ----------------------------------------------------------------

extent = [-ni * dx, ni * dx, -nj * dy, nj * dy]

# parity about the equator for the two fields that carry the rotation
par_y = -1 if ANTISYMM_W else +1

# v_phi is built as a magnitude, so when it is reflected oddly the lower half
# plane comes out negative. that is the physical picture - the hemispheres turn
# opposite ways - but it needs a diverging map about zero rather than a
# sequential one running up from 0, or everything below the equator reads as
# empty black
if ANTISYMM_W:
    v_cmap, v_lims = "RdBu_r", dict(vmin=-v_max, vmax=v_max)
else:
    v_cmap, v_lims = "inferno", dict(vmin=0.0, vmax=v_max)

# sequential for the magnitudes, diverging about zero for the signed field.
# the house scripts reach for "rainbow", which banks contrast at arbitrary
# values and invents edges the data does not have - these are perceptually
# uniform instead, so equal colour steps mean equal field steps
W_lim = np.abs(W).max()

# v_phys is drawn as a magnitude on its own scale, like psi and rho, so the
# structure of the rotation profile is visible rather than only how close it
# comes to 1. that costs comparability between runs - a panel here says nothing
# about the light speed limit on its own - which is the trade field_plot_lum.py
# takes the other way
panels = [
    (reflect(psi), r"$\psi$", "viridis",
     dict(vmin=psi.min(), vmax=psi.max()), "neither"),
    (reflect(W, parity_x=-1, parity_y=par_y), r"$\mathcal{W}$", "RdBu_r",
     dict(vmin=-W_lim, vmax=W_lim), "neither"),
    (reflect(rho_raw), r"$\rho$", "magma",
     dict(vmin=0.0, vmax=rho_raw.max()), "neither"),
    (reflect(v_phys, parity_y=par_y),
     r"$v_\phi$", v_cmap,
     v_lims, "neither"),
]

# matplotlib >= 3.9 hands back a COPY from get_cmap, so set_bad on that is
# silently discarded - the modified map has to be built and passed to imshow
panels = [(f, t,
           matplotlib.colormaps[c].with_extremes(bad="0.85")
           if isinstance(c, str) else c,
           l, e)
          for f, t, c, l, e in panels]

# 2 x 2 rather than 1 x 4. the window is square, so four square panels pack
# into a square figure with far less wasted paper than a 19-inch strip, and
# each one lands about 25% larger on the page
fig, axes = plt.subplots(2, 2, figsize=(11, 9.8), constrained_layout=True,
                         sharex=True, sharey=True)

for k, (ax, (field, title, cmap, lims, extend)) in enumerate(zip(axes.flat, panels)):
    im = ax.imshow(field, origin="lower", extent=extent, cmap=cmap,
                   aspect="equal", interpolation="nearest", **lims)

    # the horizon, on every panel. white over a dark stroke so it survives
    # both the light and the dark end of all four colour maps
    if hx is not None:
        ax.plot(hx, hy, color="white", linewidth=1.6,
                path_effects=[withStroke(linewidth=3.2, foreground="black")],
                label="horizon")

    # the luminal surface itself, where there is one to draw
    if v_max > 1.0:
        ax.contour(np.linspace(extent[0], extent[1], field.shape[1]),
                   np.linspace(extent[2], extent[3], field.shape[0]),
                   reflect(v_phys, parity_y=par_y), levels=[1.0],
                   colors="red", linewidths=1.2)

    # imshow sets the limits to the extent it painted, but ax.plot above can
    # reach past it, and matplotlib answers that by AUTOSCALING rather than
    # clipping. that leaves the field sitting in the middle of a larger, empty
    # axes - which looks like the plot zoomed out and the data went missing,
    # when in fact the white is bare figure background around a window that was
    # simply too small. pin the limits to what was actually loaded, so anything
    # outside is cropped and the picture can never again lie about its extent
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    ax.set_title(title)

    # sharex/sharey already hide the inner TICK labels; the axis labels have to
    # be skipped by hand, or they repeat down the middle of the figure
    if k // 2 == 1:
        ax.set_xlabel(r"$r$")
    if k % 2 == 0:
        ax.set_ylabel(r"$z$")

    ax.axhline(0.0, color="white", linewidth=0.4, alpha=0.25)
    ax.axvline(0.0, color="white", linewidth=0.4, alpha=0.25)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, extend=extend)

if hx is not None:
    axes.flat[0].legend(loc="upper right", fontsize=9, framealpha=0.8)

if SHOW_HEADER:
    fig.suptitle(
        f"{os.path.basename(run)}    "
        rf"$\rho_0 = {par['rho_0']:.2f}$,  $v_0 = {par['v_0']:.1f}$,  "
        rf"$a = {par['a']:.3f}$,  $b = {par['b']:.3f}$,  "
        rf"$e = {par['ecc']:.3f}$,  $M_{{\rm ADM}} = {par['M_ADM']:.4f}$"
        + "\n"
        + (f"horizon: {SURF_FILE}" if hx is not None else "no horizon"),
        fontsize=11,
    )

fig.savefig(os.path.join(HERE, OUT_STEM + ".png"))
fig.savefig(os.path.join(HERE, OUT_STEM + ".pdf"))
print(f"\nSaved {OUT_STEM}.png and {OUT_STEM}.pdf in vis/")

plt.show()
