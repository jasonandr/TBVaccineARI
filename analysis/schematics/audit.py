"""Audit: does the figure draw exactly the transitions the model implements?"""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import make_schematics as M

def which(pt):
    """Nearest compartment to a point (arrows start/end on box perimeters)."""
    best, bd = None, 1e9
    for b in M.BOXREG:
        x0,y0,x1,y1 = b.bounds
        dx = max(x0-pt[0], 0, pt[0]-x1); dy = max(y0-pt[1], 0, pt[1]-y1)
        d = (dx*dx+dy*dy)**0.5
        if d < bd: best, bd = b.name, d
    return best

def inventory(kind, alt):
    M.REG.clear(); M.BOXREG.clear()
    fig = plt.figure(); ax = fig.add_subplot(111)
    xlim, ylim = M.PANELS[kind](ax, alt)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    out = []
    for _, p0, p1, rad, t in M.REG:
        lab = t.get_text().replace("\n", " ") if t is not None else ""
        out.append((which(p0), which(p1), lab))
    plt.close(fig)
    return sorted(out)

# ---- transitions implemented in R/model.R (lines 127-162), by hand, with source
MODEL = {
 ("leaky", False): {
   ("S","E"),("E","L"),("E","D"),("L","E"),("L","D")},
 ("leaky", True): {
   ("S","E"),("E","L"),("E","D"),("L","E"),("L","D"),
   ("E","D")},                                    # + R.fastprogression (same pair)
 ("aon", False): {
   ("S","E"),("E","L"),("E","D"),("L","E"),("L","D"),
   ("E_v","L_v"),("E_v","D"),("L_v","D"),         # EP/LP dynamics
   ("S","P"),("E","P"),("L","P"),                 # enrolment allocation -> P
   ("E","E_v"),("L","L_v"),                       # enrolment allocation -> EP/LP
   ("L_v","E_v")},                                # drawn struck-through = rate 0
 ("aon", True): None,                             # = aon False + ("E","D") fastprog
}
MODEL[("aon", True)] = MODEL[("aon", False)] | {("E","D")}

for kind in ("leaky","aon"):
    for alt in (False, True):
        inv = inventory(kind, alt)
        pairs = [(a,b) for a,b,_ in inv]
        want = MODEL[(kind,alt)]
        extra  = sorted(set(pairs) - want)
        missing= sorted(want - set(pairs))
        # ("E","D") legitimately appears twice in alt models (progression + fastprog)
        ndup = sum(1 for x in pairs if x == ("E","D"))
        exp_dup = 2 if alt else 1
        tag = f'{kind:5s} alt={str(alt):5s}'
        ok = not extra and not missing and ndup == exp_dup
        print(f'{tag}  {"PASS" if ok else "FAIL"}   arrows={len(pairs)}  E->D count={ndup} (expect {exp_dup})')
        if extra:   print('     EXTRA (in figure, not in model):', extra)
        if missing: print('     MISSING (in model, not in figure):', missing)
        for a,b,l in inv:
            print(f'        {a:4s} -> {b:4s}   {l}')
        print()
