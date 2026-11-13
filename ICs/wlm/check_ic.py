#!/usr/bin/env python3
"""Check a DiskModel Gadget-1 IC (single or multi-file): masses, IDs, structure, Q, v_c."""
import sys, glob, struct
import numpy as np

G = 4.30091e-6  # kpc (km/s)^2 / Msun
UM = 1e10       # code mass unit [Msun]


def blk(f, dtype, n, ncomp=1):
    (l,) = struct.unpack('i', f.read(4))
    a = np.frombuffer(f.read(l), dtype=dtype)
    (l2,) = struct.unpack('i', f.read(4))
    assert l == l2 and a.size == n * ncomp, (l, l2, a.size, n * ncomp)
    return a.reshape(n, ncomp) if ncomp > 1 else a


def read(base):
    files = sorted(glob.glob(base + '.*[0-9]'), key=lambda s: int(s.rsplit('.', 1)[1])) \
        if not glob.glob(base) else [base]
    out = {t: {'pos': [], 'vel': [], 'id': [], 'm': [], 'u': []} for t in range(6)}
    for fn in files:
        with open(fn, 'rb') as f:
            (l,) = struct.unpack('i', f.read(4)); h = f.read(256); f.read(4)
            npart = struct.unpack('6i', h[:24]); mass = struct.unpack('6d', h[24:72])
            n = sum(npart)
            pos = blk(f, 'f4', n, 3); vel = blk(f, 'f4', n, 3); ids = blk(f, 'i4', n)
            nm = sum(npart[t] for t in range(6) if mass[t] == 0)
            m = blk(f, 'f4', nm) if nm else None
            u = blk(f, 'f4', npart[0]) if npart[0] else None
            o = om = 0
            for t in range(6):
                k = npart[t]
                if k == 0: continue
                out[t]['pos'].append(pos[o:o + k]); out[t]['vel'].append(vel[o:o + k]); out[t]['id'].append(ids[o:o + k])
                if mass[t] == 0: out[t]['m'].append(m[om:om + k]); om += k
                else: out[t]['m'].append(np.full(k, mass[t], 'f4'))
                if t == 0: out[t]['u'].append(u)
                o += k
    P = {}
    for t, d in out.items():
        if d['pos']:
            P[t] = {k: np.concatenate(v) for k, v in d.items() if v}
    return P, files


def expfit(R, m, rmin, rmax):
    edges = np.linspace(rmin, rmax, 15); c = 0.5 * (edges[1:] + edges[:-1])
    S = np.histogram(R, edges, weights=m)[0] / (np.pi * np.diff(edges ** 2))
    ok = S > 0
    return -1.0 / np.polyfit(c[ok], np.log(S[ok]), 1)[0]


def main(base, vcfile=None):
    P, files = read(base)
    names = {0: 'gas', 1: 'DM', 2: 'disc stars', 3: 'bulge', 5: 'BH'}
    print(f'files: {len(files)}')
    allid = np.concatenate([P[t]['id'] for t in P]).astype(np.int64)
    print(f'IDs: N={allid.size} unique={np.unique(allid).size == allid.size} min={allid.min()} max={allid.max()}')
    mtot = 0; pc = np.zeros(3); mom = np.zeros(3)
    for t in P:
        m = P[t]['m'].astype(np.float64) * UM; mtot += m.sum()
        pc += (P[t]['pos'] * m[:, None]).sum(0); mom += (P[t]['vel'] * m[:, None]).sum(0)
        r = np.linalg.norm(P[t]['pos'], axis=1)
        print(f'{names[t]:>10}: N={m.size:>11d}  m_part={m[0]:.4g} (min {m.min():.4g} max {m.max():.4g})  M={m.sum():.5g} Msun  r_max={r.max():.1f} kpc  r_99={np.percentile(r, 99):.1f}')
    print(f'total M={mtot:.5g} Msun  COM={pc / mtot} kpc  Vcom={mom / mtot} km/s')

    # spherical enclosed mass -> v_c
    rall = np.concatenate([np.linalg.norm(P[t]['pos'], axis=1) for t in P])
    mall = np.concatenate([P[t]['m'].astype(np.float64) * UM for t in P])
    o = np.argsort(rall); rs = rall[o]; ms = np.cumsum(mall[o])
    def vc(R): return np.sqrt(G * np.interp(R, rs, ms) / R)
    print('v_c (spherical M(<r)) [km/s]: ' + '  '.join(f'{R}:{vc(R):.1f}' for R in (0.5, 1, 2, 3, 5, 10, 20, 45.4)))
    r200 = 45.4; print(f'M(<45.4 kpc) = {np.interp(r200, rs, ms):.4g} Msun; DM(<45.4)={P[1]["m"][np.linalg.norm(P[1]["pos"],axis=1)<r200].astype(np.float64).sum()*UM:.4g}')

    g = P[0]; s = P.get(2)
    Rg = np.hypot(g['pos'][:, 0], g['pos'][:, 1]); zg = g['pos'][:, 2]
    print(f'gas: exp scale length (0.5-5 kpc) = {expfit(Rg, g["m"], 0.5, 5):.3f} kpc; median|z| (R<3) = {np.median(np.abs(zg[Rg < 3])):.3f} kpc')
    u = g['u'].astype(np.float64)  # (km/s)^2
    for mu in (0.59, 1.22):
        T = (2 / 3) * u * 1e10 * mu * 1.6726e-24 / 1.380649e-16
        print(f'gas T (mu={mu}): median {np.median(T):.4g} K  min {T.min():.4g} max {T.max():.4g}')
    if s is not None:
        Rs = np.hypot(s['pos'][:, 0], s['pos'][:, 1]); zs = s['pos'][:, 2]
        print(f'stars: exp scale length (0.5-5 kpc) = {expfit(Rs, s["m"], 0.5, 5):.3f} kpc; median|z| = {np.median(np.abs(zs)):.3f} kpc (Gaussian sigma 0.7 -> 0.472)')

    # Toomre Q from particles
    edges = np.arange(0.25, 6.01, 0.5); cR = 0.5 * (edges[1:] + edges[:-1])
    vcc = vc(cR); dv = np.gradient(vcc, cR); kap = np.sqrt(2 * vcc / cR * (dv + vcc / cR))  # km/s/kpc
    area = np.pi * np.diff(edges ** 2)
    Sg = np.histogram(Rg, edges, weights=g['m'].astype(np.float64) * UM)[0] / area
    cs = np.sqrt((2 / 3) * np.array([np.median(u[(Rg >= a) & (Rg < b)]) if np.any((Rg >= a) & (Rg < b)) else np.nan for a, b in zip(edges[:-1], edges[1:])]))
    Qg = cs * kap / (np.pi * G * Sg)
    phi = np.arctan2(g['pos'][:, 1], g['pos'][:, 0])
    vphi_g = -g['vel'][:, 0] * np.sin(phi) + g['vel'][:, 1] * np.cos(phi)
    vrot = np.array([np.mean(vphi_g[(Rg >= a) & (Rg < b)]) for a, b in zip(edges[:-1], edges[1:])])
    rows = []
    if s is not None:
        ps = np.arctan2(s['pos'][:, 1], s['pos'][:, 0])
        vR = s['vel'][:, 0] * np.cos(ps) + s['vel'][:, 1] * np.sin(ps)
        vps = -s['vel'][:, 0] * np.sin(ps) + s['vel'][:, 1] * np.cos(ps)
        Ss = np.histogram(Rs, edges, weights=s['m'].astype(np.float64) * UM)[0] / area
        sR = np.array([np.std(vR[(Rs >= a) & (Rs < b)]) for a, b in zip(edges[:-1], edges[1:])])
        sz = np.array([np.std(s['vel'][(Rs >= a) & (Rs < b), 2]) for a, b in zip(edges[:-1], edges[1:])])
        vrs = np.array([np.mean(vps[(Rs >= a) & (Rs < b)]) for a, b in zip(edges[:-1], edges[1:])])
        Qs = sR * kap / (3.36 * G * Ss)
    print(' R[kpc]  vc   vrot_gas  vrot_st  sigR_st sigz_st  Sig_g[Msun/pc2] Sig_s   Q_gas   Q_st   Q_tot(WS94)')
    for i, R in enumerate(cR):
        if s is not None:
            Qt = 1 / (1 / Qg[i] + 1 / Qs[i])
            print(f'{R:6.2f} {vcc[i]:6.1f} {vrot[i]:8.1f} {vrs[i]:8.1f} {sR[i]:7.1f} {sz[i]:7.1f} {Sg[i]/1e6:12.3f} {Ss[i]/1e6:8.3f} {Qg[i]:7.2f} {Qs[i]:6.2f} {Qt:8.2f}')
    if vcfile:
        txt = open(vcfile).read()
        if 'Toomre Q:' in txt:
            q = np.array([float(x) for x in txt.split('Toomre Q:')[1].split()])
            q = q[np.isfinite(q) & (q > 0)]
            print(f'generator .vc Toomre Q: min {q.min():.3g}  median {np.median(q):.3g}  ({q.size} radii)')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
