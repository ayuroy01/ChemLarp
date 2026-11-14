#!/usr/bin/env python3
# Usage: python -I cif2sdf.py <file.cif> <out.sdf> [charge]   (needs RDKit; run on CIFs in a scratch folder)
# Asymmetric unit of a small-molecule CIF -> SDF (bond orders from RDKit's DetermineBondOrders on the
# metal-free part; each Fe/Ru bonded to the ring carbons within 2.4 Å).
import sys, math, shlex
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D

def parse(path):
    toks = []
    for line in open(path):
        if line.startswith(';'):
            continue
        line = line.split('#')[0] if not line.lstrip().startswith("'") else line
        try:
            toks += shlex.split(line, posix=True)
        except ValueError:
            toks += line.split()
    cell = {}
    atoms = []
    i = 0
    while i < len(toks):
        t = toks[i]
        if t.startswith('_cell_') and i + 1 < len(toks):
            cell[t] = float(toks[i + 1].split('(')[0]); i += 2; continue
        if t == 'loop_':
            j = i + 1; cols = []
            while j < len(toks) and toks[j].startswith('_'):
                cols.append(toks[j]); j += 1
            vals = []
            while j < len(toks) and not toks[j].startswith('_') and toks[j] != 'loop_' and not toks[j].startswith('data_'):
                vals.append(toks[j]); j += 1
            if '_atom_site_fract_x' in cols:
                n = len(cols)
                for k in range(0, len(vals) - n + 1, n):
                    row = dict(zip(cols, vals[k:k + n]))
                    if row.get('_atom_site_disorder_group', '.') not in ('.', '1'):
                        continue
                    atoms.append(row)
            i = j; continue
        i += 1
    return cell, atoms

def cart(cell):
    a, b, c = cell['_cell_length_a'], cell['_cell_length_b'], cell['_cell_length_c']
    al, be, ga = (math.radians(cell[k]) for k in ('_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'))
    v = math.sqrt(1 - math.cos(al)**2 - math.cos(be)**2 - math.cos(ga)**2 + 2*math.cos(al)*math.cos(be)*math.cos(ga))
    M = [[a, b*math.cos(ga), c*math.cos(be)],
         [0, b*math.sin(ga), c*(math.cos(al) - math.cos(be)*math.cos(ga))/math.sin(ga)],
         [0, 0, c*v/math.sin(ga)]]
    return lambda f: [sum(M[r][k]*f[k] for k in range(3)) for r in range(3)]

cell, rows = parse(sys.argv[1])
to = cart(cell)
metals = {'Fe', 'Ru', 'Cr', 'Co', 'Ni', 'Os', 'Mn'}
els, xyz = [], []
for r in rows:
    el = r.get('_atom_site_type_symbol') or ''.join(ch for ch in r['_atom_site_label'] if ch.isalpha())[:2]
    el = el.rstrip('+-0123456789'); el = el[0].upper() + el[1:].lower()
    f = [float(r[k].split('(')[0]) for k in ('_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z')]
    els.append(el); xyz.append(to(f))
org = [i for i, e in enumerate(els) if e not in metals]
met = [i for i, e in enumerate(els) if e in metals]
# Cp rings: carbons within 2.4 A of a metal. Each ring gets a capping H on one CH carbon (pointing away from
# the metal) so the metal-free part is neutral cyclopentadiene; bond orders are guessed with charge 0, then
# the cap is removed and that carbon becomes the Cp- carbanion.
caps = []
pts = [xyz[i] for i in org]
pels = [els[i] for i in org]
hapto = {}
for i in met:
    for k, j in enumerate(org):
        if els[j] == 'C' and math.dist(xyz[i], xyz[j]) < 2.4: hapto[k] = i
def neighbours(k):
    return [q for q in range(len(pts)) if q != k and math.dist(pts[k], pts[q]) < (1.25 if 'H' in (pels[k], pels[q]) else 1.75)]
done = set()
for k in sorted(hapto):
    if k in done: continue
    ring, stack = set(), [k]
    while stack:
        u = stack.pop()
        if u in ring: continue
        ring.add(u)
        stack += [q for q in neighbours(u) if q in hapto and hapto[q] == hapto[k]]
    done |= ring
    if len(ring) != 5: continue
    pick = next((q for q in sorted(ring) if any(pels[h] == 'H' for h in neighbours(q))), sorted(ring)[0])
    mpos = xyz[hapto[k]]
    d = [pts[pick][t] - mpos[t] for t in range(3)]; L = math.sqrt(sum(x*x for x in d))
    caps.append((pick, len(pts)))
    pts.append([pts[pick][t] + 1.09 * d[t] / L for t in range(3)]); pels.append('H')
block = f"{len(pts)}\n\n" + "".join(f"{pels[i]} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n" for i, p in enumerate(pts))
m = Chem.MolFromXYZBlock(block)
rdDetermineBonds.DetermineConnectivity(m)
rdDetermineBonds.DetermineBondOrders(m, charge=int(sys.argv[3]) if len(sys.argv) > 3 else 0)
rw = Chem.RWMol(m)
for c, h in sorted(caps, key=lambda t: -t[1]):
    rw.RemoveAtom(h)
    rw.GetAtomWithIdx(c).SetFormalCharge(-1)
    rw.GetAtomWithIdx(c).SetNoImplicit(True)
for i in met:
    idx = rw.AddAtom(Chem.Atom(els[i]))
    rw.GetAtomWithIdx(idx).SetFormalCharge(sum(1 for c, h in caps if hapto[c] == i))
    conf = rw.GetConformer()
    conf.SetAtomPosition(idx, Point3D(*xyz[i]))
    for k, j in enumerate(org):
        if els[j] == 'C' and math.dist(xyz[i], xyz[j]) < 2.4:
            rw.AddBond(idx, k, Chem.BondType.SINGLE)
mol = rw.GetMol()
mol.SetProp('_Name', sys.argv[1].split('/')[-1])
w = Chem.SDWriter(sys.argv[2]); w.SetKekulize(True); 
try:
    w.write(mol)
except Exception:
    mol.UpdatePropertyCache(strict=False); w.write(mol)
w.close()
print(sys.argv[2], len(org), 'organic atoms,', len(met), 'metal')
