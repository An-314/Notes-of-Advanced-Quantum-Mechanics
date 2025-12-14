#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第9次作业],
  author: "AnZrew",
  time: "2025年12月",
)

#exercise(subname: [Sakurai 6.4])[
  Consider a potential
  $
    V = 0 "for" r > R, V = V_0 = "constant" "for" r < R,
  $
  where $V_0$ may be positive or negative. Using the method of partial waves, show that for $abs(V_0) << E = (hbar^2 k^2)/(2m)$ and $k R << 1$, the differential cross section is isotropic and that the total cross section is given by
  $
    sigma_"tot" = (16pi)/9 (m^2 V_0^2 R^6)/hbar^4
  $
  Suppose the energy is raised slightly. Show that the angular distribution can then be written as
  $
    dv(sigma, Omega) = vb(A) + vb(B)cos(theta)
  $
  Obtain an approximate expression for $B/A$.
]

#exercise(subname: [Steck 16.3])[
  The basic delta-scattering potential in 3D is
  $
    V = beta delta^3(vb(r))
  $
  without any regularization for $1/r$ divergent states. Starting with the T-matrix relation (16.43),
  $
    T^+(E) = V + V G_0^+ (E) T^+(E);
  $
  show that the T-matrix element in the momentum basis may be written
  $
    braket(vb(k), hat(T)^+ (E), vb(k)) = 1/(2 pi)^3 (1/beta - I)^(-1)
  $
  in terms of the integral
  $
    I := 1/(2 pi)^3 integral dd(vb(k)) 1/(E - (hbar^2 k^2)/(2m) + i 0^+)
  $
  Then argue that the integral diverges, and thus the scattering amplitude and cross section vanish for the (unregularized) 3D delta scatterer.
  This is quite a result, considering the first Born approximation for the scattering amplitude is finite, and the second term in the expansion is divergent. The problem is that the delta potential is an unphysical
  idealization because it ignores important structure of the scattering potential at short length scales (large $k$). By cutting off the integral $I$ at large $k$ one can obtain a finite (but cutoff-dependent) result for the scattering amplitude.#footnote[
    Indrajit Mitra, Ananda DasGupta, and Binayak Dutta-Roy, ''Regularization and renormalization in scattering from Dirac delta potentials,'' American Journal of Physics 66, 1101 (1998) (#link("http://dx.doi.org/10.1119/1.19051")[doi: 10.1119/1.19051]) regularize the divergence and give what appears to be a finite result for the T -matrix element and scattering cross section in terms of the binding energy of a bound state. However, the binding energies tend to infinity for the 3D delta well, so their results are consistent with the vanishing cross section.
  ]
]
