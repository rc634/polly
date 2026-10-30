"""
Meridional field plot for one ROTATING polly run, laid out and styled to match
field_plot_static.py: the metric down the left, the physical matter down the
right,

    psi     rho_D
    W       psi^2 r v^phi

with the horizon, where the run has one, over all four panels.

The matter panels are the physical dust quantities, not the rescaled ones that
source the solve. In the paper's notation

    rho_N   = psi^5 rho_D / (1 - (v^phi)^2 psi^4 r^2)      rho-eff.dat
    v_N^phi = psi^5 v^phi                                  v-eff.dat

so rho_D is rho-raw.dat and v^phi is v-raw.dat. The paper writes the physical
velocity as plain v^phi, and so does the panel title; the comments and printout
below call it v^phi_D, to sit with rho_D as physical.

v^phi_D itself is an angular velocity, and under differential rotation
(v_N^phi ~ v_0 / r) it runs off like 1/r into the axis, so a plot of it is one
bright column of cells and nothing else. The panel shows the speed
psi^2 r v^phi_D instead - the quantity in rho_N's denominator, which has to
stay below 1 for rho_N to exist, and the one Source::sanity_check tests. It is
finite on the axis, and on a fixed 0 .. 1 scale (V_MAX) it reads directly as a
fraction of the speed of light, the same way in every figure. It has a colour
map of its own, dark to cyan - see speed_cmap for why none of matplotlib's do.

It is cut to the matter, as in field_plot.py: v-raw.dat keeps its v_0 floor out
into the vacuum, where the speed means nothing. The cut is to NaN, drawn black
to match the empty space around rho_D.

W is ODD about the rotation axis - see field_plot.py for why - and the speed is
even about it.

Run with the venv that has numpy + matplotlib:
    ./plt-env/bin/python field_plot_spin.py
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

# the run directory, relative to this file.
#
# paper-sph-v spins the rho_0 = 3.10 ball of paper-sph up through v_0 = 0 .. 4,
# rigidly in the rescaled variable: v-eff.dat is v_0 in every cell. v0 to v3
# all keep the horizon the static ball already had at rho_0 = 3.10; v4 is the
# first to lose it - the chaser gives up and pollax_surf.dat is zero-filled.
#
# the differentially rotating runs go on the command line:
#   ../batch/paper-diffrotsph/v18
#   ../batch/paper-diffrotpan-HV/v100
RUN_DIR = "../batch/paper-sph-v/rho-3_10-v4"

# which horizon to overlay - see field_plot.py for why pollax is the default
SURF_FILE = "pollax_surf.dat"

# is W ODD about the equator, as Grid::fill_all_ghosts_antisymmW imposes? see
# field_plot.py. if so, W and the speed both reverse in the lower half plane
ANTISYMM_W = False

# top of the speed colour scale. 1 pins it at the speed of light, so every
# figure reads on the same scale; None stretches it to the run's own maximum
V_MAX = 1.0

# half-width of the plotted window, in code units, or None to size it from the
# run: whichever is bigger of the matter and the horizon
HALF_WIDTH = None

# padding around whatever the auto window ends up framing
WINDOW_MARGIN = 1.12

OUT_STEM = "field_plot_spin"

# ----------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))

# a directory on the command line overrides RUN_DIR, as in field_plot.py. the
# sweep directory goes into the name, as in field_plot_static.py
if len(sys.argv) == 2:
    RUN_DIR = sys.argv[1]
    OUT_STEM = ("field_plot_spin_"
                + os.path.basename(os.path.dirname(os.path.normpath(RUN_DIR)))
                + "_" + os.path.basename(os.path.normpath(RUN_DIR)))
elif len(sys.argv) != 1:
    sys.exit("usage: field_plot_spin.py [RUN_DIR]")

run = RUN_DIR if os.path.isabs(RUN_DIR) else os.path.join(HERE, RUN_DIR)

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


def load_window(name, nj, ni):
    """Read the first nj rows and ni columns of a quadrant .dat file."""
    path = os.path.join(run, name)
    arr = np.loadtxt(path, delimiter=",", max_rows=nj, usecols=range(ni))
    print(f"  {name:12s} -> {arr.shape[0]} x {arr.shape[1]}")
    return arr


def reflect(arr, parity_x=+1, parity_y=+1):
    """
    Mirror a quadrant out to the full meridional plane.

    parity_x = -1 for a field that is odd about the rotation axis (W),
    parity_y = -1 for one odd about the equator (W and the speed, with
    ANTISYMM_W). Returned array pairs with origin="lower".
    """
    upper = np.hstack([parity_x * np.fliplr(arr), arr])
    lower = parity_y * np.flipud(upper)
    return np.vstack([lower, upper])


def speed_cmap(n=256):
    """
    The speed panel's colour map: near-black at 0 to a bright mint at 1.

    matplotlib's perceptually uniform maps that start dark - magma, inferno,
    plasma, cividis - all finish in yellow or cream, which left the speed panel
    next to rho_D in two near-identical maps; viridis is psi's. This one is
    built the way they were, in CIELAB: the lightness L* climbs linearly, so
    equal steps in speed are equal steps in perceived brightness (and it still
    reads in greyscale), while the hue turns from plum through blue and teal to
    mint green.

    The hue has to keep turning all the way up. A version that sat on cyan for
    the top third was just as uniform on paper, but saturated cyans are hard to
    tell apart by eye, and v18 - which lives at 0.6 .. 0.98 - came out as one
    flat disc with its two fast lobes gone.

    The chroma follows the most vivid colour the sRGB gamut allows at each step,
    smoothed so it has no kinks. Clipping each step to the gamut on its own
    leaves kinks in the chroma, and those read as bright bands that are not in
    the data; scaling one fixed chroma curve down to fit the tightest step
    instead leaves the rest of the map washed out.
    """
    t = np.linspace(0.0, 1.0, n)
    L = 4.0 + (95.0 - 4.0) * t
    h = np.radians(330.0 + (145.0 - 330.0) * t)

    def lab_to_linear_rgb(L, C, h):
        # CIELAB (D65 white) -> XYZ -> linear sRGB, for arrays of n steps
        fy = (L + 16.0) / 116.0
        fx = fy + C * np.cos(h) / 500.0
        fz = fy - C * np.sin(h) / 200.0
        finv = lambda f: np.where(f > 6/29, f**3, 3 * (6/29)**2 * (f - 4/29))
        xyz = np.array([0.95047 * finv(fx), finv(fy), 1.08883 * finv(fz)])
        return np.array([[ 3.2406, -1.5372, -0.4986],
                         [-0.9689,  1.8758,  0.0415],
                         [ 0.0557, -0.2040,  1.0570]]) @ xyz

    # the largest chroma each (L*, h) can take inside the gamut, bisected for
    # all n steps at once
    lo, hi = np.zeros(n), np.full(n, 150.0)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        ok = np.all((lambda c: (c >= 0.0) & (c <= 1.0))(lab_to_linear_rgb(L, mid, h)),
                    axis=0)
        lo, hi = np.where(ok, mid, lo), np.where(ok, hi, mid)
    c_max = lo

    # just inside the gamut, capped so the middle does not go garish, then
    # gaussian-smoothed (sigma 14 of the 256 steps) back under the gamut a few
    # times over, and finally scaled down only if a smoothed bump poked out
    k = np.exp(-0.5 * (np.arange(-42, 43) / 14.0)**2)
    k /= k.sum()
    smooth = lambda y: np.convolve(np.pad(y, 42, mode="edge"), k, mode="valid")
    C = np.minimum(0.92 * c_max, 60.0)
    for _ in range(6):
        C = smooth(np.minimum(C, 0.92 * c_max))
    C *= min(1.0, np.min(0.99 * c_max[C > 0] / C[C > 0]))

    rgb = np.clip(lab_to_linear_rgb(L, C, h), 0.0, 1.0).T
    srgb = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb**(1/2.4) - 0.055)
    return matplotlib.colors.ListedColormap(srgb, name="speed")


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


# ---- grid geometry ---------------------------------------------------------

par = read_params(os.path.join(run, "params.dat"))
nx, ny = int(par["nx"]), int(par["ny"])
dx = (par["xU"] - par["xL"]) / nx
dy = (par["yU"] - par["yL"]) / ny

if par["v_0"] == 0.0:
    print(f"WARNING  : {RUN_DIR} has v_0 = 0, so W and the speed are blank - "
          f"field_plot_static.py is the plot for a static run")

# ---- the horizon, read before the window is sized -------------------------

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

half = HALF_WIDTH
if half is None:
    reach = max(par["a"], par["b"])
    if hx is not None:
        reach = max(reach, np.abs(hx).max(), np.abs(hy).max())
    half = WINDOW_MARGIN * reach

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
rho_D = load_window("rho-raw.dat", nj, ni)
rho_eff = load_window("rho-eff.dat", nj, ni)
v_D = load_window("v-raw.dat", nj, ni)

# the matter mask comes from rho_eff, not rho_D - rho_D goes negative once the
# matter is superluminal, and would drop exactly the cells worth seeing
matter = rho_eff > 0.0

# the speed psi^2 r v^phi_D, with r the cell-centred cylindrical radius as
# Field::get_x computes it. zero outside the matter here, for the luminal
# contour - the panel's NaN copy is made below
x_centres = (np.arange(ni) + 0.5) * dx + par["xL"]
speed = np.where(matter, np.abs(v_D * x_centres[None, :] * psi**2), 0.0)

v_max = speed.max()
n_bad = int((speed > 1.0).sum())
n_neg = int((rho_D < 0.0).sum())

print(f"\n  psi      {psi.min():.6f} .. {psi.max():.6f}")
print(f"  W        {W.min():.6g} .. {W.max():.6g}")
print(f"  rho_D    {rho_D.min():.6g} .. {rho_D.max():.6g}")
print(f"  v^phi_D  {v_D[matter].min():.6g} .. {v_D[matter].max():.6g}"
      f"   (in the matter)")
print(f"  psi^2 r v^phi_D max = {v_max:.6f}"
      f"   cells above 1 = {n_bad}"
      f"{'   *** SUPERLUMINAL ***' if n_bad else ''}")
if n_neg:
    print(f"  *** rho_D is NEGATIVE in {n_neg} cells - superluminal matter,"
          f" min = {rho_D.min():.6g}.  the panel clips these to 0 ***")

# ---- panels ----------------------------------------------------------------

extent = [-ni * dx, ni * dx, -nj * dy, nj * dy]

par_y = -1 if ANTISYMM_W else +1

W_lim = np.abs(W).max()
v_top = v_max if V_MAX is None else V_MAX

# NaN outside the matter, drawn in the map's bad colour - see the top
speed_panel = np.where(matter, speed, np.nan)

# the speed is reflected oddly about the equator with ANTISYMM_W, so it needs a
# diverging map about zero there - see field_plot.py
if ANTISYMM_W:
    v_cmap, v_lims = "RdBu_r", dict(vmin=-v_top, vmax=v_top)
else:
    v_cmap, v_lims = speed_cmap(), dict(vmin=0.0, vmax=v_top)

# matter past the top of a pinned scale - superluminal, at V_MAX = 1 - gets an
# arrow on the colour bar rather than quietly saturating into the top colour
v_extend = ("neither" if v_max <= v_top
            else "both" if ANTISYMM_W else "max")

# metric down the left, matter down the right
panels = [
    [(reflect(psi), r"$\psi$", "viridis",
      dict(vmin=psi.min(), vmax=psi.max()), "neither"),
     (reflect(rho_D), r"$\rho_\mathrm{D}$", "magma",
      dict(vmin=0.0, vmax=rho_D.max()), "neither")],
    [(reflect(W, parity_x=-1, parity_y=par_y), r"$\mathcal{W}$", "RdBu_r",
      dict(vmin=-W_lim, vmax=W_lim), "neither"),
     (reflect(speed_panel, parity_y=par_y), r"$\psi^2 r\, v^\phi$",
      v_cmap, v_lims, v_extend)],
]

fig, axes = plt.subplots(2, 2, figsize=(11.6, 10.4), constrained_layout=True,
                         sharex=True, sharey=True)

for r in range(2):
    for c in range(2):
        field, title, cmap, lims, extend = panels[r][c]
        ax = axes[r, c]

        # matplotlib >= 3.9 hands back a COPY from get_cmap, so the modified
        # map has to be built and passed to imshow. only the speed has NaN
        # cells, its vacuum, and black there matches the rho_D panel. the speed
        # map is built here rather than registered, so it arrives as an object
        cmap = (matplotlib.colormaps[cmap] if isinstance(cmap, str)
                else cmap).with_extremes(bad="black")
        im = ax.imshow(field, origin="lower", extent=extent, cmap=cmap,
                       aspect="equal", interpolation="nearest", **lims)

        # white over a dark stroke, so it survives both ends of all four maps
        if hx is not None:
            ax.plot(hx, hy, color="white", linewidth=1.6,
                    path_effects=[withStroke(linewidth=3.2, foreground="black")],
                    label="horizon")

        # the luminal surface, where there is one to draw
        if n_bad:
            ax.contour(np.linspace(extent[0], extent[1], field.shape[1]),
                       np.linspace(extent[2], extent[3], field.shape[0]),
                       reflect(speed, parity_y=par_y), levels=[1.0],
                       colors="red", linewidths=1.2)

        # ax.plot can reach past the imshow extent and matplotlib autoscales
        # rather than clipping - pin the limits, as in field_plot.py
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])

        # every panel is its own field here, so every panel gets a title - the
        # pad is the one field_plot_static.py needed at this title size
        ax.set_title(title, pad=12)
        if r == 1:
            ax.set_xlabel(r"$r$")
        if c == 0:
            ax.set_ylabel(r"$z$")

        ax.axhline(0.0, color="white", linewidth=0.4, alpha=0.25)
        ax.axvline(0.0, color="white", linewidth=0.4, alpha=0.25)

        # four different fields, so four bars
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, extend=extend)

# shared axes, so this sets every panel - see field_plot_static.py
axes[0, 0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
axes[0, 0].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))

if hx is not None:
    axes[0, 0].legend(loc="upper right", fontsize=13.5, framealpha=0.8)

fig.savefig(os.path.join(HERE, OUT_STEM + ".png"))
fig.savefig(os.path.join(HERE, OUT_STEM + ".pdf"))
print(f"\nSaved {OUT_STEM}.png and {OUT_STEM}.pdf in vis/")

plt.show()
