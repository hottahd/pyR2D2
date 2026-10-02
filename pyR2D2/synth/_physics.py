"""合成が使う物理定数とプランク関数（ODF-radiation の odfgen.constants・odfgen.means と同じ値・同じ式）。

R2D2plus-input の opacity/（odfgen）と結果をビット単位で揃えるため、値も式も写したまま変えない。
"""
import numpy as np

c_light = 2.99792458e10  # 光速 [cm/s]
h_planck = 6.62607015e-27  # プランク定数 [erg s]
k_boltz = 1.380649e-16  # ボルツマン定数 [erg/K]


def planck_lambda(lam_A, T):
    """B_lambda(T) [erg/s/cm^2/sr/cm]。lam は Angstrom。"""
    lam = np.asarray(lam_A, dtype=float) * 1.0e-8  # cm
    x = h_planck * c_light / (lam * k_boltz * T)
    with np.errstate(over="ignore"):
        return (2.0 * h_planck * c_light**2 / lam**5) / np.expm1(x)
