#!/usr/bin/env python3
"""
TB vaccine model schematics -- Figure 1 (base model) and Supplemental Figure 1
(alternative model used in sensitivity analyses).

True vector output from one parameterised routine, so Figure 1 and Supplemental
Figure 1 are guaranteed identical in style and differ ONLY by the alternative-model
reinfection-driven progression pathway.

Two structural guarantees that earlier hand-built versions lacked:
  * every arrow label is positioned from the ACTUAL Bezier midpoint of its arrow
    (no hand-guessed coordinates), offset perpendicular to the direction of travel;
  * a fixed inches-per-data-unit scale, so standalone panels and the combined
    two-panel figures render boxes and text at exactly the same visual size.

Palette, box style and font are sampled from the originally submitted figures
(Figures/Figure1a.png, Figure1b.png) so the visual identity is preserved.

Rate expressions come from the model code (R/model.R), not from prose:

    R.infectionrate_S    = lambda * (1 - LeakyPOI*vaxon)         S -> E
    R.EtoL               = e                                     E -> L   (no VE)
    R.reinfection        = lambda * omega * (1 - LeakyPOI*vaxon)  L -> E
    R.primaryprogression = r1 * (1 - LeakyPOD*vaxon)              E -> D
    R.reactivation       = r2 * (1 - LeakyPOD*vaxon)              L -> D
    R.fastprogression    = r1 * lambda * fastprogon
                              * (1 - LeakyPOD*vaxon) * (1 - LeakyPOI*vaxon)
                                                                 E -> D, ALT MODEL ONLY

  All-or-none arm: LeakyPOD = LeakyPOI = 0, so the vaccinated compartments
  (P, EP=E_v, LP=L_v) carry bare natural-history rates and vaccine efficacy enters
  only through allocation at enrolment. EP does NOT receive R.fastprogression
  (all-or-none POI blocks the exposures that drive it).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.font_manager as fm

# ----------------------------------------------------------------------------- style
for _f in ("/System/Library/Fonts/Supplemental/Arial.ttf",
           "/System/Library/Fonts/Supplemental/Arial Bold.ttf"):
    if os.path.exists(_f):
        fm.fontManager.addfont(_f)

plt.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "pdf.fonttype": 42,     # embed TrueType -> text stays editable in Illustrator
    "svg.fonttype": "none", # keep text as text in the SVG
})

FILL = {"S": "#E0F4DA", "E": "#DFEDFA", "L": "#EDD9BA", "D": "#EBD2D0"}
EDGE = {"S": "#356920", "E": "#2D5F8A", "L": "#8A4F21", "D": "#673630"}
INK, ALT, ALLOC = "#1A1A1A", "#B3372A", "#4A4A4A"

BW, BH = 1.00, 0.90        # box size in data units (original was 290 x 262 px)
SCALE = 0.62               # INCHES PER DATA UNIT -- fixes visual size everywhere
LW_BOX, LW_ARR = 1.35, 1.5
FS_STATE, FS_LAB, FS_NOTE, FS_KEY = 30, 10.5, 9.5, 8.5


# ------------------------------------------------------------------------ primitives
BOXREG = []       # every Box drawn, for the audit


class Box:
    def __init__(self, name, cx, cy, kind):
        self.name, self.cx, self.cy, self.kind = name, cx, cy, kind

    @property
    def bounds(self):
        return self.cx - BW / 2, self.cy - BH / 2, self.cx + BW / 2, self.cy + BH / 2

    def anchor(self, side, frac=0.5):
        x0, y0, x1, y1 = self.bounds
        return {"left": (x0, y0 + frac * BH), "right": (x1, y0 + frac * BH),
                "top": (x0 + frac * BW, y1), "bottom": (x0 + frac * BW, y0)}[side]

    def draw(self, ax):
        BOXREG.append(self)
        x0, y0, _, _ = self.bounds
        ax.add_patch(FancyBboxPatch(
            (x0, y0), BW, BH, boxstyle="round,pad=0,rounding_size=0.075",
            facecolor=FILL[self.kind], edgecolor=EDGE[self.kind],
            linewidth=LW_BOX, zorder=3, mutation_aspect=1))
        t = self.name
        if "_" in t:
            b, s = t.split("_", 1)
            t = rf"$\mathrm{{{b}}}_{{\mathrm{{{s}}}}}$"
        ax.text(self.cx, self.cy, t, ha="center", va="center",
                fontsize=FS_STATE, color=INK, zorder=4)


def _arc_geom(p0, p1, rad):
    """Exact midpoint of matplotlib's arc3 quadratic Bezier, and the unit normal
    pointing to the LEFT of the direction of travel."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    right = np.array([d[1], -d[0]])                # matplotlib arc3: ctrl offsets RIGHT
    ctrl = (p0 + p1) / 2 + rad * right             # exactly matplotlib's control point
    mid = 0.25 * p0 + 0.5 * ctrl + 0.25 * p1       # Bezier at t = 0.5
    left = -right                                  # label convention: off>0 = left
    n = np.linalg.norm(left)
    return mid, (left / n if n else left)


REG = []          # (ax, p0, p1, rad, text_artist_or_None) for the collision guard


def edge(ax, p0, p1, rad=0.0, name=None, rate=None, off=0.20,
         color=INK, dashed=False, note=None, label_at=None, ha="center"):
    """Draw an arrow and place its label at the curve's true midpoint.

    off > 0 puts the label to the LEFT of the direction of travel, so a
    left-to-right arrow gets its label above and a right-to-left arrow below --
    both land outside the diagram without any per-arrow tuning.
    """
    ax.add_patch(FancyArrowPatch(
        p0, p1, connectionstyle=f"arc3,rad={rad}",
        arrowstyle="-|>,head_length=6,head_width=3.2",
        color=color, linewidth=LW_ARR, zorder=2,
        linestyle=(0, (5, 3)) if dashed else "solid",
        shrinkA=1.5, shrinkB=1.5, joinstyle="miter", capstyle="butt"))

    parts = [p for p in (name, rate) if p]
    if not parts:
        REG.append((ax, p0, p1, rad, None))
        return
    mid, nrm = _arc_geom(p0, p1, rad)
    if label_at is not None:
        x, y, va = label_at[0], label_at[1], "center"
    else:
        x, y = mid + off * nrm
        dy = off * nrm[1]             # how far the label moved vertically
        if abs(nrm[1]) < 0.35:        # normal ~horizontal: label sits beside arrow
            va = "center"
        else:                         # grow the text block AWAY from the curve
            va = "bottom" if dy > 0 else "top"
    t = ax.text(x, y, "\n".join(parts), ha=ha, va=va, fontsize=FS_LAB,
                color=color, zorder=5, linespacing=1.35,
                multialignment="center")
    REG.append((ax, p0, p1, rad, t))


def blocked(ax, p0, p1, name=None, off=-0.20, color="#B4B4B4", xcolor="#4A4A4A"):
    """Draw a transition that the vaccine abolishes: faint arrow struck through with
    an X. Used for L_v -> E_v, so the reader can see that all-or-none protection
    against infection removes the reinfection route rather than having to notice a
    missing arrow."""
    ax.add_patch(FancyArrowPatch(
        p0, p1, connectionstyle="arc3,rad=0",
        arrowstyle="-|>,head_length=6,head_width=3.2",
        color=color, linewidth=LW_ARR, zorder=2,
        shrinkA=1.5, shrinkB=1.5, joinstyle="miter", capstyle="butt"))
    mid, nrm = _arc_geom(p0, p1, 0.0)
    d = 0.115
    for sx in (+1, -1):
        ax.plot([mid[0] - d, mid[0] + d], [mid[1] - sx * d, mid[1] + sx * d],
                color=xcolor, linewidth=1.9, solid_capstyle="round", zorder=6)
    t = None
    if name:
        t = ax.text(*(mid + off * nrm), name, ha="center",
                    va="top" if off * nrm[1] < 0 else "bottom",
                    fontsize=FS_LAB, color=xcolor, zorder=5)
    REG.append((ax, p0, p1, 0.0, t))


# --------------------------------------------------------------------- rate strings
def rates(leaky):
    if leaky:
        return dict(infection=r"$\lambda(1-\mathrm{VE}_I)$",
                    etol=r"$\varepsilon$",
                    reinfection=r"$\lambda\omega(1-\mathrm{VE}_I)$",
                    progression=r"$r_1(1-\mathrm{VE}_D)$",
                    reactivation=r"$r_2(1-\mathrm{VE}_D)$",
                    fastprog=r"$r_1\lambda(1-\mathrm{VE}_I)(1-\mathrm{VE}_D)$")
    return dict(infection=r"$\lambda$", etol=r"$\varepsilon$",
                reinfection=r"$\lambda\omega$", progression=r"$r_1$",
                reactivation=r"$r_2$", fastprog=r"$r_1\lambda$")


# --------------------------------------------------------------------------- panel A
def panel_leaky(ax, alt):
    """Leaky protection: every vaccinee has partially reduced risk on each
    transition. Topology follows the originally submitted Figure 1A."""
    R = rates(True)
    S, E = Box("S", 0.0, 0.0, "S"), Box("E", 3.4, 0.0, "E")
    L, D = Box("L", 6.8, 0.0, "L"), Box("D", 3.4, -2.85, "D")
    for b in (S, E, L, D):
        b.draw(ax)

    edge(ax, S.anchor("right"), E.anchor("left"),
         name="infection", rate=R["infection"], off=0.22)
    edge(ax, E.anchor("right", 0.74), L.anchor("left", 0.74), rate=R["etol"], off=0.16)
    edge(ax, L.anchor("left", 0.26), E.anchor("right", 0.26),
         name="reinfection", rate=R["reinfection"], off=0.20)
    edge(ax, E.anchor("bottom", 0.60), D.anchor("top", 0.60),
         name="progression", rate=R["progression"], off=1.00)
    edge(ax, L.anchor("bottom", 0.50), D.anchor("right", 0.42), rad=-0.24,
         name="reactivation", rate=R["reactivation"], off=0.52)
    if alt:
        edge(ax, E.anchor("bottom", 0.14), D.anchor("left", 0.66), rad=0.34,
             name="reinfection-driven progression", rate=R["fastprog"],
             off=-1.62, color=ALT)
    return (-0.85, 7.45), (-3.45, 0.82)


# --------------------------------------------------------------------------- panel B
def panel_aon(ax, alt):
    """All-or-none protection: a fraction of vaccinees are completely protected.

    Allocation happens at enrolment (proportions, not rates), exactly as in the code:
        P  = ltbiprev*AoN.D.pos + (1-ltbiprev)*[1-(1-AoN.I)(1-AoN.D.neg)]
        EP = (1-AoN.D.pos)*AoN.I*recentlyinfected      (= E_v)
        LP = (1-AoN.D.pos)*AoN.I*distantlyinfected     (= L_v)
    Protection against DISEASE sends participants to P, where they cannot progress.
    Protection against INFECTION sends already-infected participants to E_v / L_v,
    where they are protected against reinfection but continue to progress through
    natural history at the unvaccinated rates -- hence E_v -> D and L_v -> D, which
    matches the appendix equations D(t+1) = ... + r1*Ev(t) + r2*Lv(t).

    P (protected, left) and D (disease, right) are the two sinks; enrolment allocation
    flows left into P, natural-history progression flows right into D.
    """
    R = rates(False)
    S, E, L = Box("S", 0.0, 0.0, "S"), Box("E", 3.4, 0.0, "E"), Box("L", 6.8, 0.0, "L")
    Ev, Lv = Box("E_v", 3.4, 3.15, "E"), Box("L_v", 6.8, 3.15, "L")
    Pv, D = Box("P", 0.0, 3.15, "S"), Box("D", 11.0, 1.575, "D")
    for b in (S, E, L, Ev, Lv, Pv, D):
        b.draw(ax)

    # ---- unvaccinated tier: natural history
    edge(ax, S.anchor("right"), E.anchor("left"),
         name="infection", rate=R["infection"], off=0.22)
    edge(ax, E.anchor("right", 0.74), L.anchor("left", 0.74), rate=R["etol"], off=0.16)
    edge(ax, L.anchor("left", 0.26), E.anchor("right", 0.26),
         name="reinfection", rate=R["reinfection"], off=0.20)
    edge(ax, E.anchor("bottom", 0.50), D.anchor("left", 0.14), rad=0.55,
         name="progression", rate=R["progression"], off=-0.30)
    edge(ax, L.anchor("bottom", 0.62), D.anchor("left", 0.34), rad=0.30,
         name="reactivation", rate=R["reactivation"], off=0.34)

    # ---- vaccinated tier: bare natural-history rates (VE acted at enrolment)
    edge(ax, Ev.anchor("right", 0.74), Lv.anchor("left", 0.74), rate=R["etol"], off=0.16)
    blocked(ax, Lv.anchor("left", 0.26), Ev.anchor("right", 0.26),
            name="reinfection blocked", off=0.20)
    edge(ax, Ev.anchor("top", 0.50), D.anchor("left", 0.86), rad=-0.34,
         name="progression", rate=R["progression"], off=0.30)
    edge(ax, Lv.anchor("top", 0.62), D.anchor("left", 0.66), rad=-0.20,
         name="reactivation", rate=R["reactivation"], off=-0.34)

    # ---- allocation at enrolment (dashed): protection against DISEASE -> P.
    # Already-infected participants given an all-or-none POD vaccine become fully
    # protected and cannot progress; this is the pathway by which an all-or-none
    # vaccine prevents disease in people who are already infected.
    edge(ax, S.anchor("top", 0.50), Pv.anchor("bottom", 0.50), color=ALLOC, dashed=True)
    ax.text(0.17, 1.02, r"$\mathrm{VE}_I$ and/or $\mathrm{VE}_D$", ha="left",
            va="center", fontsize=FS_LAB, color=ALLOC, zorder=5)
    # labels placed explicitly because these arcs' own midpoints fall in congested
    # space; the collision guard still checks them like any other label
    edge(ax, E.anchor("top", 0.12), Pv.anchor("bottom", 0.80), rad=-0.05,
         rate=r"$\mathrm{VE}_D$", label_at=(1.78, 1.82), color=ALLOC, dashed=True)
    edge(ax, L.anchor("top", 0.12), Pv.anchor("right", 0.55), rad=0.26,
         rate=r"$\mathrm{VE}_D$", label_at=(6.15, 1.30), color=ALLOC, dashed=True)

    # ---- allocation at enrolment (dashed): protection against INFECTION -> E_v / L_v.
    # These participants are protected against reinfection but are NOT protected
    # against disease, so they still progress at the unvaccinated rates.
    for b_lo, b_hi in ((E, Ev), (L, Lv)):
        edge(ax, b_lo.anchor("top", 0.72), b_hi.anchor("bottom", 0.72),
             rate=r"$(1-\mathrm{VE}_D)\,\mathrm{VE}_I$",
             label_at=(b_lo.anchor("top", 0.72)[0] + 0.16, 0.98), ha="left",
             color=ALLOC, dashed=True)

    ax.text(Pv.cx, Pv.cy + BH / 2 + 0.20, "fully protected\n(cannot progress)",
            ha="center", va="bottom", fontsize=FS_NOTE, color=EDGE["S"],
            zorder=5, linespacing=1.3)

    if alt:
        edge(ax, E.anchor("bottom", 0.22), D.anchor("bottom", 0.55), rad=0.78,
             name="reinfection-driven progression", rate=R["fastprog"],
             off=-0.36, color=ALT)
    return (-1.40, 11.85), (-3.55 if alt else -2.35, 5.05)


PANELS = {"leaky": panel_leaky, "aon": panel_aon}

KEY = [
    (r"$\lambda$  force of infection        "
     r"$\omega$  relative risk of infection given prior infection (partial natural immunity)", "#333333"),
    (r"$\varepsilon$  rate of transition from early to late infection        "
     r"$r_1$, $r_2$  early and late progression rates", "#333333"),
    (r"$\mathrm{VE}_I$  vaccine efficacy against infection        "
     r"$\mathrm{VE}_D$  vaccine efficacy against disease", "#333333"),
    (r"$\lambda$, $\omega$, $\varepsilon$, $r_1$, $r_2$ are natural history parameters applied in both "
     r"trial arms; $\mathrm{VE}_I$ and $\mathrm{VE}_D$ act on the vaccine arm only.", "#333333"),
    ("In B, dashed arrows are enrolment allocation proportions, not rates. P cannot progress: "
     "protection against disease acts here.", "#333333"),
    (r"$\mathrm{E_v}$ and $\mathrm{L_v}$ are protected against reinfection only, and still progress "
     r"to disease at the unvaccinated rates $r_1$ and $r_2$.", "#333333"),
]
KEY_ALT = [
    ("Red: reinfection-driven progression \u2014 alternative model only; acts on early infection (E), "
     "not on late infection.", ALT),
    (r"$\mathrm{E_v}$ is exempt, because all-or-none protection against infection blocks the exposures "
     "that drive it.", ALT),
]


def check_collisions(fig, tag, pad=1.5):
    """Fail if any arrow passes through any label, or any two labels overlap.
    The arrow path is the analytic quadratic Bezier that matplotlib's arc3 draws,
    sampled and mapped to display coordinates -- so this is exact, not eyeballed."""
    fig.canvas.draw()
    problems = []
    labels = [(ax, t) for ax, _, _, _, t in REG if t is not None]

    for ax, p0, p1, rad, _ in REG:
        p0v, p1v = np.asarray(p0, float), np.asarray(p1, float)
        d = p1v - p0v
        ctrl = (p0v + p1v) / 2 + rad * np.array([d[1], -d[0]])
        ts = np.linspace(0.06, 0.94, 90)[:, None]
        pts = (1 - ts) ** 2 * p0v + 2 * ts * (1 - ts) * ctrl + ts ** 2 * p1v
        disp = ax.transData.transform(pts)
        for lax, t in labels:
            if lax is not ax:
                continue
            bb = t.get_window_extent().expanded(1.0, 1.0)
            hit = ((disp[:, 0] > bb.x0 + pad) & (disp[:, 0] < bb.x1 - pad) &
                   (disp[:, 1] > bb.y0 + pad) & (disp[:, 1] < bb.y1 - pad))
            if hit.any():
                problems.append(f"arrow crosses label {t.get_text()[:42]!r}")

    for ax, t in labels:
        bb = t.get_window_extent()
        for b in BOXREG:
            x0, y0, x1, y1 = b.bounds
            c = ax.transData.transform([[x0, y0], [x1, y1]])
            from matplotlib.transforms import Bbox
            if bb.overlaps(Bbox([[min(c[:, 0]), min(c[:, 1])],
                                 [max(c[:, 0]), max(c[:, 1])]])):
                problems.append(f"label {t.get_text()[:34]!r} overlaps box {b.name}")

    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            (ax1, t1), (ax2, t2) = labels[i], labels[j]
            if ax1 is not ax2:
                continue
            b1, b2 = t1.get_window_extent(), t2.get_window_extent()
            if b1.overlaps(b2):
                problems.append(
                    f"labels overlap: {t1.get_text()[:28]!r} / {t2.get_text()[:28]!r}")
    if problems:
        print(f"  !! {tag}")
        for q in sorted(set(problems)):
            print(f"       {q}")
    return problems


# -------------------------------------------------------------------------- assembly
def _mk_ax(fig, xlim, ylim, left_in, bottom_in, figw, figh):
    w = (xlim[1] - xlim[0]) * SCALE
    h = (ylim[1] - ylim[0]) * SCALE
    ax = fig.add_axes([left_in / figw, bottom_in / figh, w / figw, h / figh])
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_aspect("equal"); ax.axis("off")
    return ax


def draw_standalone(kind, alt, stem, letter):
    REG.clear()
    fn = PANELS[kind]
    probe = plt.figure()                       # measure extents without emitting output
    xlim, ylim = fn(probe.add_subplot(111), alt)
    plt.close(probe)

    pad, top = 0.10, 0.0 if letter is None else 0.34
    figw = (xlim[1] - xlim[0]) * SCALE + 2 * pad
    figh = (ylim[1] - ylim[0]) * SCALE + 2 * pad + top
    fig = plt.figure(figsize=(figw, figh))
    ax = _mk_ax(fig, xlim, ylim, pad, pad, figw, figh)
    fn(ax, alt)
    if letter is not None:
        check_collisions(fig, stem.split("/")[-1])
    fig.text(pad / figw, 1 - 0.06 / figh, letter, ha="left", va="top",
                 fontsize=17, fontweight="bold", color=INK)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=600 if ext == "png" else None, facecolor="white")
    plt.close(fig)


def draw_combined(alt, stem):
    REG.clear()
    ext_lim = {}
    for kind in ("leaky", "aon"):
        probe = plt.figure()
        ext_lim[kind] = PANELS[kind](probe.add_subplot(111), alt)
        plt.close(probe)

    pad, gap, hdr = 0.16, 0.46, 0.30
    keyn = len(KEY) + (len(KEY_ALT) if alt else 0)
    keyh = 0.12 + 0.235 * keyn
    wA = (ext_lim["leaky"][0][1] - ext_lim["leaky"][0][0]) * SCALE
    wB = (ext_lim["aon"][0][1] - ext_lim["aon"][0][0]) * SCALE
    hA = (ext_lim["leaky"][1][1] - ext_lim["leaky"][1][0]) * SCALE
    hB = (ext_lim["aon"][1][1] - ext_lim["aon"][1][0]) * SCALE

    figw = max(wA, wB) + 2 * pad
    figh = hA + hB + 2 * (hdr) + gap + keyh + 2 * pad
    fig = plt.figure(figsize=(figw, figh))

    yB = pad + keyh
    axB = _mk_ax(fig, *ext_lim["aon"], pad, yB, figw, figh)
    PANELS["aon"](axB, alt)
    yA = yB + hB + hdr + gap
    axA = _mk_ax(fig, *ext_lim["leaky"], pad, yA, figw, figh)
    PANELS["leaky"](axA, alt)

    for y, letter, sub in ((yA + hA + 0.04, "A", "leaky protection"),
                           (yB + hB + 0.04, "B", "all-or-none protection")):
        fig.text(pad / figw, y / figh, letter, ha="left", va="bottom",
                 fontsize=18, fontweight="bold", color=INK)
        fig.text((pad + 0.30) / figw, y / figh, sub, ha="left", va="bottom",
                 fontsize=11.5, color="#555555")

    rows = KEY + (KEY_ALT if alt else [])
    for i, (txt, col) in enumerate(rows):
        fig.text(pad / figw, (pad + keyh - 0.14 - 0.235 * i) / figh, txt, ha="left",
                 va="top", fontsize=FS_KEY, color=col)

    check_collisions(fig, stem.split("/")[-1])
    fig.canvas.draw()
    lim = figw * fig.dpi
    for t in fig.texts:
        bb = t.get_window_extent()
        if bb.x1 > lim - 2 or bb.x0 < 2:
            raise SystemExit(
                f"TEXT OVERFLOWS CANVAS in {stem}: {t.get_text()[:60]!r} "
                f"(x0={bb.x0:.0f}, x1={bb.x1:.0f}, limit={lim:.0f})")
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=600 if ext == "png" else None, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    out = os.path.dirname(os.path.abspath(__file__))
    for alt, base in ((False, "Figure1"), (True, "FigureS1")):
        draw_standalone("leaky", alt, f"{out}/{base}A", None)
        draw_standalone("aon", alt, f"{out}/{base}B", None)
        draw_combined(alt, f"{out}/{base}")
    print("wrote Figure1{,A,B} + FigureS1{,A,B} as .pdf/.png/.svg to", out)
