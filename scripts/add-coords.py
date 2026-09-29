#!/usr/bin/env python3
"""Make a self-contained restart file from a checkpoint without coordinates.

  add-coords.py <checkpoint> <file with coordinates> <output>

nekRS writes the mesh coordinates only into a job's first checkpoint, but an
interpolated restart onto another mesh (startFrom = <file>+int) needs them.
This copies X from <file with coordinates> (the first checkpoint of the same
job, so the element order matches) into <checkpoint>. Both files must use the
Nek5000 #std format with FP32 data.
"""
import sys

import numpy as np


def layout(path):
    with open(path, "rb") as f:
        hdr = f.read(132)
    w = hdr.decode().split()
    nx, ny, nz, nel = int(w[2]), int(w[3]), int(w[4]), int(w[5])
    code = w[11]
    ncomp = (3 if "X" in code else 0) + (3 if "U" in code else 0) + (1 if "P" in code else 0)
    if "S" in code:
        ncomp += int(code[code.index("S") + 1:code.index("S") + 3])
    return hdr, nel, nx*ny*nz, code, ncomp


def main(ckpt, with_x, out):
    h1, nel, nxyz, code1, ncomp1 = layout(ckpt)
    h0, nel0, nxyz0, code0, ncomp0 = layout(with_x)
    assert code0 == "X" + code1 and (nel, nxyz) == (nel0, nxyz0), (code0, code1)
    head = 132 + 4 + 4*nel
    block = 4*nel*nxyz
    d1 = np.memmap(ckpt, dtype=np.uint8, mode="r")
    d0 = np.memmap(with_x, dtype=np.uint8, mode="r")
    assert len(d1) == head + ncomp1*(block + 8*nel)
    assert len(d0) == head + ncomp0*(block + 8*nel)
    assert np.array_equal(d0[132:head], d1[132:head]), "element order differs"

    time0, time1 = h0.decode().split()[7], h1.decode().split()[7]
    hdr = h0.decode().replace(time0, time1, 1).encode()
    assert len(hdr) == 132
    meta0, meta1 = head + ncomp0*block, head + ncomp1*block
    with open(out, "wb") as f:
        f.write(hdr)
        f.write(d1[132:head])                    # endian marker and element ids
        f.write(d0[head:head + 3*block])         # X
        f.write(d1[head:meta1])                  # U, P, S
        f.write(d0[meta0:meta0 + 3*8*nel])       # X metadata (min/max per element)
        f.write(d1[meta1:])                      # U, P, S metadata

    d2 = np.memmap(out, dtype=np.uint8, mode="r")
    assert len(d2) == len(d0)
    assert np.array_equal(d2[head:head + 3*block], d0[head:head + 3*block])
    assert np.array_equal(d2[head + 3*block:head + ncomp0*block], d1[head:meta1])
    assert np.array_equal(d2[head + ncomp0*block + 24*nel:], d1[meta1:])
    print("{}: {} with coordinates from {} (t = {})".format(out, ckpt, with_x, float(time1)))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
