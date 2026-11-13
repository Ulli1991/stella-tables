#!/usr/bin/env python3
"""Convert a single-file DiskModel Gadget-1 IC to a GIZMO HDF5 IC (uint64 IDs, explicit masses)."""
import sys, struct
import numpy as np
import h5py


def main(src, dst):
    raw = np.memmap(src, dtype='u1', mode='r')
    off = 0

    def rec():
        nonlocal off
        (l,) = struct.unpack('i', raw[off:off + 4].tobytes())
        start = off + 4
        (l2,) = struct.unpack('i', raw[start + l:start + l + 4].tobytes())
        assert l == l2
        off = start + l + 4
        return start, l

    s, l = rec(); h = raw[s:s + 256].tobytes()
    npart = np.array(struct.unpack('6i', h[:24]), dtype=np.int64)
    mass = struct.unpack('6d', h[24:72])
    assert struct.unpack('i', h[124:128])[0] == 1, 'multi-file input not supported'
    n = int(npart.sum())
    sp, lp = rec(); sv, lv = rec(); si, li = rec()
    assert lp == 12 * n and lv == 12 * n and li == 4 * n
    nm = int(sum(npart[t] for t in range(6) if mass[t] == 0))
    assert nm == 0, 'per-particle mass block not handled'
    su, lu = rec() if npart[0] else (None, 0)
    assert off == raw.size, 'trailing bytes'
    pos = np.ndarray((n, 3), np.float32, raw, sp); vel = np.ndarray((n, 3), np.float32, raw, sv)
    ids = np.ndarray(n, np.int32, raw, si)
    u = np.ndarray(int(npart[0]), np.float32, raw, su) if npart[0] else None
    first = np.concatenate([[0], np.cumsum(npart)])
    with h5py.File(dst, 'w') as f:
        hd = f.create_group('Header').attrs
        hd['NumPart_ThisFile'] = npart.astype(np.uint32)
        hd['NumPart_Total'] = (npart & 0xffffffff).astype(np.uint32)
        hd['NumPart_Total_HighWord'] = (npart >> 32).astype(np.uint32)
        hd['MassTable'] = np.zeros(6)
        hd['Time'] = 0.0; hd['Redshift'] = 0.0; hd['BoxSize'] = 0.0
        hd['NumFilesPerSnapshot'] = 1
        hd['Omega0'] = 0.0; hd['OmegaLambda'] = 0.0; hd['HubbleParam'] = 1.0
        for k in ('Flag_Sfr', 'Flag_Cooling', 'Flag_StellarAge', 'Flag_Metals', 'Flag_Feedback', 'Flag_DoublePrecision'):
            hd[k] = 0
        hd['Flag_IC_Info'] = 0
        hd['UnitLength_In_CGS'] = 3.085678e21; hd['UnitMass_In_CGS'] = 1.989e43; hd['UnitVelocity_In_CGS'] = 1e5
        for t in range(6):
            if npart[t] == 0: continue
            a, b = first[t], first[t + 1]
            g = f.create_group(f'PartType{t}')
            g.create_dataset('Coordinates', data=pos[a:b])
            g.create_dataset('Velocities', data=vel[a:b])
            g.create_dataset('ParticleIDs', data=ids[a:b].astype(np.uint64))
            g.create_dataset('Masses', data=np.full(b - a, mass[t], np.float32))
            if t == 0:
                g.create_dataset('InternalEnergy', data=u)
    print(f'wrote {dst}: npart={npart.tolist()} masses={mass}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
