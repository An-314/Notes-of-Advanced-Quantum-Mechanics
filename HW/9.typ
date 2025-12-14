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
    dv(sigma, Omega) = A + B cos(theta)
  $
  Obtain an approximate expression for $B/A$.
]

#solution[
  分波法的散射振幅
  $
    f(theta) & = 1/k sum_(l=0)^oo (2 l + 1) e^(i delta_l) sin delta_l P_l (cos theta) \
  $
  在$k R << 1$时，相移的行为是
  $
    delta_l prop k^(2 l + 1)
  $
  所以$f$的主要贡献来自$l = 0$项，即
  $
    f(theta) & = 1/k e^(i delta_0) sin delta_0 \
  $
  这就意味着散射截面
  $
    dv(sigma, Omega) & = |f(theta)|^2 \
  $
  在$k R << 1$时，散射截面是各向同性的，即与$theta$无关。

  进一步地，在$abs(V_0) << E$时，在$r<R$区域
  $
    k' = sqrt((2m (E - V_0))/hbar^2) = k sqrt(1 - V_0/E)
  $
  只考虑$l=0$项，径向解取
  $
    u_0 (r) = cases(
      C sin(k' r) "  " r<R,
      sin(k r + delta_0) "  " r>R,
    )
  $
  考虑到在$r=R$处的连续性条件$u'/u$有
  $
    k' cot(k' R) = k cot(k R + delta_0) <=> tan(k R + delta_0) = k/k' tan(k' R)
  $
  令$k R = x, k' R = k'/k x$，对$x$展开到三阶，有
  $
    tan(x + delta_0) & = x + delta_0 + (x^3)/3 + O(x^5) \
    k/k' tan(k'/k x) & = k/k' (k'/k x + (k'/k)^3 (x^3)/3) + O(x^5) \
  $
  从而有
  $
    delta_0 = ((k'/k)^2 - 1) x^3/3 = - V_0/E (k R)^3/3 = - (2 m V_0 k R^3)/(3 hbar^2)
  $
  从而有
  $
    sigma_"tot" & = integral dd(Omega) dv(sigma, Omega) = 4pi abs(f(0))^2 = 4pi abs(delta_0/k)^2 \
                & = 4pi ((2 m V_0 R^3)/(3 hbar^2))^2 = (16pi)/9 (m^2 V_0^2 R^6)/hbar^4
  $
  #newpara()
  能量略升高：加入$delta_1$项，此时
  $
    f(theta) & = 1/k (e^(i delta_0) sin delta_0 + 3 e^(i delta_1) sin delta_1 cos(theta)) \
  $
  于是
  $
    dv(sigma, Omega) & = abs(f(theta))^2 = abs(delta_0/k)^2 + 9 abs(delta_1/k)^2 cos^2 theta + 6 Re((delta_0 delta_1)/k^2) cos(theta) \
    & approx abs(delta_0/k)^2 + 6 (delta_0 delta_1)/k^2 cos(theta)\
    &= A + B cos(theta)
  $
  从而
  $
    B/A & = 6 delta_1/delta_0
  $
  下面我们计算$delta_1$，对于任意$l$，在无穷远处
  $
    u_l (r) = j_l (k r) cos delta_l - n_l (k r) sin delta_l
  $
  $r<R$内的解为
  $
    u_l (r) = A j_l (k' r)
  $
  考虑边界条件$u'/u$连续有
  $
    A j_l (k' R) = j_l (k R) cos delta_l - n_l (k R) sin delta_l\
    A k' j'_l (k' R) = k j'_l (k R) cos delta_l - k n'_l (k R) sin delta_l
  $
  以及递推关系
  $
    j'_l (k' R) = l (j_l (k' R)) /(k' R) - j_(l + 1) (k' R)
  $
  有
  $
    tan delta_l & = ((k' R)^2 j_l (k R) - (2l + 3) j_(l + 1) (k R)) /((k' R)^2 n_l (k R) - (2l + 3) n_(l + 1) (k R)) \
                & = (-2 m V_0 R^2)/(hbar^2 (2l+3)) ((2^l l!)/(2l + 1)!)^2 (k R)^(2l+1)
  $
  $l=1$有
  $
    delta_1 = - (2 m V_0 R^2)/(5 hbar^2) (2/6)^2 (k R)^3 = - V_0/E (k R)^5/45
  $
  从而
  $
    B/A = 6 delta_1/delta_0 = 2/5 (k R)^2
  $
]

#exercise(subname: [Steck 16.3])[
  The basic delta-scattering potential in 3D is
  $
    V = beta delta^3(vb(r))
  $
  without any regularization for $1/r$ divergent states. Starting with the T-matrix relation
  $
    T^+(E) = V + V G_0^+ (E) T^+(E);
  $
  show that the T-matrix element in the momentum basis may be written
  $
    braket(vb(k), hat(T)^+ (E), vb(k)') = 1/(2 pi)^3 (1/beta - I)^(-1)
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

#proof[
  考虑到
  $
    braket(vb(r), vb(k)) = 1/(2pi)^(3/2) e^(i vb(k) dot vb(r))
  $
  则
  $
    braket(vb(k), hat(V) (E), vb(k)') &= 1/(2 pi)^3 integral dd(vb(r)) integral dd(vb(r)') braket(vb(k), vb(r)) braket(vb(r), hat(V) (E), vb(r)') braket(vb(r)', vb(k)') \
    & = beta braket(vb(k), 0) braket(0, vb(k)') = beta/(2 pi)^3
  $
  考虑到
  $
    hat(T)^+ = hat(V) + hat(V) hat(G)_0^+ hat(T)^+
  $
  就有
  $
    braket(vb(k), hat(T)^+, vb(k)') & = braket(vb(k), hat(V), vb(k)') + braket(vb(k), hat(V) hat(G)_0^+ hat(T)^+, vb(k)') \
    &= braket(vb(k), hat(V), vb(k)') + integral dd(vb(k)'') integral dd(vb(k)''') braket(vb(k), hat(V), vb(k)'') braket(vb(k)'', hat(G)_0^+, vb(k)''') braket(vb(k)''', hat(T)^+, vb(k)') \
    &= beta/(2 pi)^3 + integral dd(vb(k)''') beta/(2 pi)^3 1/(E - (hbar^2 k'''^2)/(2m) + i 0^+) braket(vb(k)''', hat(T)^+, vb(k)')\
  $
  该方程右侧对$vb(k)$的依赖仅通过$braket(vb(k)''', hat(T)^+, vb(k)')$存在，即不依赖于$vb(k)$，由于$hat(T)^+$是幺正算符，因此也不依赖于$vb(k)''$。因此，我们可以将$braket(vb(k)''', hat(T)^+, vb(k)')$写为$tau$，即
  $
    tau = beta/(2 pi)^3 + beta/(2 pi)^3 tau integral dd(vb(k)) 1/(E - (hbar^2 k^2)/(2m) + i 0^+)
  $
  从而
  $
    tau = beta/(2 pi)^3 (1 - integral dd(vb(k)) 1/(E - (hbar^2 k^2)/(2m) + i 0^+))^(-1)
  $
  #newpara()
  进一步我们看$I$的敛散性。考虑$E = (hbar^2 k^2)/(2m)$，则
  $
    I & = 1/(2 pi)^3 integral dd(vb(q)) 1/(hbar^2/(2m) (k^2 - q^2) + i 0^+) \
      & = (2m)/hbar^2 1/(2 pi)^3 integral dd(vb(q)) 1/(k^2 - q^2 + i 0^+) \
      & = (2m)/hbar^2 (4pi)/(2 pi)^3 integral_0^oo dd(q) q^2/(k^2 - q^2 + i 0^+) \
  $
  很显然最后一个式子的积分是发散的。于是得到
  $
    hat(T)^+ (E) = 1/(2 pi)^3 (1/beta - I)^(-1) -> 0
  $
  散射振幅$f$和$T$的矩阵元成正比，从而
  $
    f(theta) -> 0 => dv(sigma, Omega)= abs(f(theta))^2 -> 0
  $
]
