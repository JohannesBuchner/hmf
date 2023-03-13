"""
Halo mass function module

"Computes" the halo mass function.

"""

from math import ceil

import numpy

"""

Copyright (c) 2017-2018, Johannes Buchner
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



cosmosim queries to get the data:

select (snapnum*10000 + floor(log10(Mvir)*10)) as bucket, snapnum, floor(log10(Mvir)*10) as M, count(*)
from MDPL2.Rockstar
group by bucket

select snapnum, zred from SMDPL.Redshifts

"""


def get_hmf(logMmin=10, logMmax=16, z=0, cumulative=False, box="small"):
    """Compute halo mass function.

    halo definition: Rockstar subhalos in MultiDark simulations.

    Units
    ------
    Masses are all in Msun h^-1. h=0.6777. Planck cosmology; see above links.
    Units are number/cGpc^3.

	Parameters
	-----------
	logMmin: int
        lowest halo mass. Units are number/Gpc^3/dlog Msun.
	logMmax: int
		highest halo mass. see logMmax.
	z: float
		redshift
	cumulative: bool
		if True, returns the density above a given mass (cumulative mass function).
	box: str
		which simulation to load.

		 * 'small' means SMDPL, see https://www.cosmosim.org/cms/simulations/smdpl/
		 * 'mid' means MDPL2, see https://www.cosmosim.org/cms/simulations/mdpl2/
		 * 'big' means BigMDPL, see https://www.cosmosim.org/cms/simulations/bigmdpl/

		 Please give appropriate credit.

	Returns
	----------
	logM: np.array
		logarithmic mass axis between logMmin and logMmax
	density: np.array
		number density
"""
    boxsize = dict(small=0.4, mid=1.0, big=2.5)[box]

    simvolume_Gpc3 = (boxsize * 0.6777) ** 3

    logMmin = int(logMmin * 10)
    logMmax = int(ceil(logMmax * 10))
    if box == "big":
        logMmin = max(115, logMmin)

    Ms = numpy.arange(logMmin, logMmax)
    hist = {}
    data = numpy.loadtxt("%s_hist2_z.csv" % box, delimiter=",", skiprows=1)
    for _, _, i, M, N in data:
        hist[(i, M)] = N
    dM = 0.1
    snapshot_idxs, zs = numpy.loadtxt("%s_z.csv" % box).transpose()
    if cumulative:
        # hist should be N and everything above
        # so propagate downward
        for i in snapshot_idxs:
            for M in Ms[::-1]:
                hist[(i, M)] = hist.get((i, M + 1), 0) + hist.get((i, M), 0)
        dM = 1

    snapshot_idx = snapshot_idxs[((zs - z) ** 2).argmin()]
    values = numpy.array(
        [hist.get((snapshot_idx, M), numpy.nan) / simvolume_Gpc3 / dM for M in Ms]
    )
    return Ms / 10.0, values


if __name__ == "__main__":
    zlist = [0, 0.5, 1, 2, 3, 4.5, 6.2, 8, 10, 15]
    import matplotlib.pyplot as plt

    for cumulative in False, True:
        for z in zlist:
            Ms, rho = get_hmf(z=z, cumulative=cumulative, box="small")
            (l1,) = plt.plot(Ms, rho, label="z=%.1f" % z)
            Ms, rho = get_hmf(z=z, cumulative=cumulative, box="mid")
            (l2,) = plt.plot(Ms, rho, "--", color=l1.get_color())
            Ms, rho = get_hmf(z=z, cumulative=cumulative, box="big")
            (l2,) = plt.plot(Ms, rho, ":", color=l1.get_color())

        plt.xlabel("Halo Mass [log $M_\odot h^{-1}$]")
        plt.yscale("log")
        plt.ylim(0.1, 1e9)
        plt.legend(loc="best")
        if not cumulative:
            plt.ylabel("Number Density $dN/dM_h$ [$\mathrm{Gpc}^{-3}\log M_\odot$]")
        else:
            plt.ylabel("Number Density N(>$M_h$) [$\mathrm{Gpc}^{-3}$]")

        plt.savefig(
            "hmf%s_M.png" % ("_cumulative" if cumulative else ""), bbox_inches="tight"
        )
        plt.close()

        Ms, N = get_hmf(z=z, cumulative=cumulative, box="small")

        for M in 10, 11, 12, 13:
            rhos = [
                get_hmf(
                    logMmin=M, logMmax=M + 0.1, z=z, cumulative=cumulative, box="small"
                )[1]
                for z in zlist
            ]
            plt.plot(zlist, rhos, label="$\log M_h=%.0f h^{-1}$" % (M))

        plt.xlabel("Redshift")
        if not cumulative:
            plt.ylabel("Number Density $dN/dM_h$ [$\mathrm{Gpc}^{-3}\log M_\odot$]")
        else:
            plt.ylabel("Number Density N(>$M_h$) [$\mathrm{Gpc}^{-3}$]")
        plt.ylim(100, 3e9)
        plt.yscale("log")
        plt.legend(loc="best")
        plt.savefig(
            "hmf%s_z.png" % ("_cumulative" if cumulative else ""), bbox_inches="tight"
        )
        plt.close()
