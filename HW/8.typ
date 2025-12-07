#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第8次作业],
  author: "AnZrew",
  time: "2025年12月",
)

#exercise(subname: [Sakurai 3.24])[
  We are to add angular momenta $j_1 = 1$ and $j_2 = 1$ to form $j = 2, 1$, and $0$ states. Using either the ladder operator method or the recursion relation, express all (nine) ${j,m}$ eigenkets in terms of $ket(j_1 j_2\; m_1 m_2)$. Write your answer as
  $
    ket(j = 1 \, m = 1) = 1/sqrt(2) ket(+ \, 0) - 1/sqrt(2) ket(0 \, +), ...
  $
  where $+$ and $0$ stand for $m_(1,2) = 1 , 0$, respectively.
]

#solution[
  考虑到
  $
    hat(J)_plus.minus = hat(J)_(1 plus.minus) + hat(J)_(2 plus.minus)
  $
  自旋 1 的单粒子作用为
  $
    hat(J)_plus.minus ket(m) = hbar sqrt(2 - m(m plus.minus 1)) ket(m plus.minus 1)
  $
  - $j=2$

    最高权态
    $
      ket(j=2\, m=2) = ket(+\, +)
    $
    应用降算符得到
    $
      2 hbar ket(2\, 1) & = hat(J)_minus ket(2\, 2) = (hat(J)_(1 minus) + hat(J)_(2 minus)) ket(+\, +) = hbar sqrt(2) ket(0\, +) + hbar sqrt(2) ket(+\, 0) \
      => ket(2\, 1) &= 1/sqrt(2) ket(0\, +) + 1/sqrt(2) ket(+\, 0)
    $
    再次应用降算符得到
    $
      sqrt(6) hbar ket(2\, 0) & = hat(J)_minus ket(2\, 1) = (hat(J)_(1 minus) + hat(J)_(2 minus)) (1/sqrt(2) ket(0\, +) + 1/sqrt(2) ket(+\, 0)) \
      & = 1/sqrt(2) hbar sqrt(2) ket(-\, +) + 1/sqrt(2) hbar sqrt(2) ket(0\, 0) + 1/sqrt(2) hbar sqrt(2) ket(0\, 0) + 1/sqrt(2) hbar sqrt(2) ket(+\, -) \
      => ket(2\, 0) & = 1/sqrt(6) ket(-\, +) + 2/sqrt(6) ket(0\, 0) + 1/sqrt(6) ket(+\, -)
    $
    再次应用降算符得到
    $
      2 hbar ket(2\, -1) & = hat(J)_minus ket(2\, 0) = (hat(J)_(1 minus) + hat(J)_(2 minus)) (1/sqrt(6) ket(-\, +) + 2/sqrt(6) ket(0\, 0) + 1/sqrt(6) ket(+\, -)) \
      & = 1/sqrt(6) hbar sqrt(2) ket(-\, 0) + 2/sqrt(6) hbar sqrt(2) ket(-\, 0) + 2/sqrt(6) hbar sqrt(2) ket(0\, -) + 1/sqrt(6) hbar sqrt(2) ket(0\, -) \
      => ket(2\, -1) & = 1/sqrt(2) ket(-\, 0) + 1/sqrt(2) ket(0\, -)
    $
    最低权态
    $
      ket(2\, -2) = ket(-\, -)
    $
  - $j=1$

    在$m=1$的子空间${ket(+\, 0), ket(0\, +)}$中，与$ket(j=2\, m=1)$正交的态为
    $
      ket(1\, 1) = 1/sqrt(2) ket(0\, +) - 1/sqrt(2) ket(+\, 0)
    $
    应用降算符得到
    $
      sqrt(2) hbar ket(1\, 0) & = hat(J)_minus ket(1\, 1) = (hat(J)_(1 minus) + hat(J)_(2 minus)) (1/sqrt(2) ket(0\, +) - 1/sqrt(2) ket(+\, 0)) \
      & = 1/sqrt(2) hbar sqrt(2) ket(-\, +) - 1/sqrt(2) hbar sqrt(2) ket(0\, 0) + 1/sqrt(2) hbar sqrt(2) ket(0\, 0) - 1/sqrt(2) hbar sqrt(2) ket(+\, -) \
      => ket(1\, 0) & = 1/sqrt(2) ket(-\, +) - 1/sqrt(2) ket(+\, -)
    $
    再次应用降算符得到
    $
      sqrt(2) hbar ket(1\, -1) & = hat(J)_minus ket(1\, 0) = (hat(J)_(1 minus) + hat(J)_(2 minus)) (1/sqrt(2) ket(-\, +) - 1/sqrt(2) ket(+\, -)) \
      & = 1/sqrt(2) hbar sqrt(2) ket(-\, 0) - 1/sqrt(2) hbar sqrt(2) ket(0\, -) \
      => ket(1\, -1) & = 1/sqrt(2) ket(-\, 0) - 1/sqrt(2) ket(0\, -)
    $

  - $j=0$

    在$m=0$的子空间${ket(-\, +), ket(0\, 0), ket(+\, -)}$中，与$ket(j=2\, m=0)$和$ket(j=1\, m=0)$正交的态为
    $
      ket(0\, 0) = 1/sqrt(3) ket(-\, +) - 1/sqrt(3) ket(0\, 0) + 1/sqrt(3) ket(+\, -)
    $

  综上可以得所有九个态为
  $
     ket(2\, 2) & = ket(+\, +) \
     ket(2\, 1) & = 1/sqrt(2) ket(0\, +) + 1/sqrt(2) ket(+\, 0) \
     ket(2\, 0) & = 1/sqrt(6) ket(-\, +) + 2/sqrt(6) ket(0\, 0) + 1/sqrt(6) ket(+\, -) \
    ket(2\, -1) & = 1/sqrt(2) ket(-\, 0) + 1/sqrt(2) ket(0\, -) \
    ket(2\, -2) & = ket(-\, -) \
     ket(1\, 1) & = 1/sqrt(2) ket(0\, +) - 1/sqrt(2) ket(+\, 0) \
     ket(1\, 0) & = 1/sqrt(2) ket(-\, +) - 1/sqrt(2) ket(+\, -) \
    ket(1\, -1) & = 1/sqrt(2) ket(-\, 0) - 1/sqrt(2) ket(0\, -) \
     ket(0\, 0) & = 1/sqrt(3) ket(-\, +) - 1/sqrt(3) ket(0\, 0) + 1/sqrt(3) ket(+\, -) \
  $
]

#exercise(subname: [Sakurai 6.2(a)])[
  Prove
  $
    sigma_"tot" tilde.eq m^2/(pi hbar^4) integral dd(vb(x), 3) integral dd(vb(x'), 3) V(r) V(r') (sin^2 k abs(vb(x) - vb(x)'))/(k^2 abs(vb(x) - vb(x'))^2)
  $
  by integrating the differential cross section computed using the first-order Born approximation.
]

#proof[
  微分散射截面为
  $
    dv(sigma, Omega) = abs(f(theta, phi))^2
  $
  Born近似给出
  $
    f^((1)) (vb(k)', vb(k)) & = - (m)/(2 pi hbar^2) integral dd(vb(r)', 3) e^(i (vb(k) - vb(k)') dot vb(r)') V(vb(r)') \
  $
  其中$abs(vb(k)) = abs(vb(k)') = k$。因此
  $
    dv(sigma, Omega) & = (m^2)/(4 pi^2 hbar^4) abs(integral dd(vb(r)', 3) e^(i (vb(k) - vb(k)') dot vb(r)') V(vb(r)'))^2 \
    & = (m^2)/(4 pi^2 hbar^4) integral dd(vb(x), 3) integral dd(vb(x'), 3) e^(i (vb(k) - vb(k)') dot (vb(x) - vb(x'))) V(vb(x)) V(vb(x'))
  $
  对立体角积分
  $
    sigma_"tot" & = integral dd(Omega) dv(sigma, Omega)\
    &= (m^2)/(4 pi^2 hbar^4) integral dd(vb(x), 3) V(vb(x)) V(vb(x')) integral dd(vb(x'), 3) integral dd(Omega) e^(i (vb(k) - vb(k)') dot (vb(x) - vb(x')))
  $
  记$vb(R) = vb(x) - vb(x')$，则
  $
    integral dd(Omega) e^(i (vb(k) - vb(k)') dot vb(R)) & =e^(i vb(k) dot vb(R)) integral_0^(2 pi) dd(phi) integral_0^pi dd(theta) sin theta e^(- i k R cos theta) \
    & = 2 pi e^(i vb(k) dot vb(R)) integral_(-1)^(1) dd(cos theta) e^(- i k R cos theta) \
    & = 2 pi e^(i vb(k) dot vb(R)) (2 sin(k R))/(k R) \
    & = 4 pi e^(i vb(k) dot vb(R)) (sin(k abs(vb(x) - vb(x'))))/(k abs(vb(x) - vb(x')))
  $
  代回到总截面表达式中
  $
    sigma_"tot" & = (m^2)/(4 pi^2 hbar^4) integral dd(vb(x), 3) V(vb(x)) integral dd(vb(x'), 3) V(vb(x')) 4 pi e^(i vb(k) dot (vb(x) - vb(x'))) (sin(k abs(vb(x) - vb(x'))))/(k abs(vb(x) - vb(x'))) \
  $
  由于势$V(r)$是球对称的，故可以对$vb(k)$做一次角向平均
  $
    e^(i vb(k) dot (vb(x) - vb(x')))_"avg" & = 1/(4 pi) integral dd(Omega_k) e^(i vb(k) dot (vb(x) - vb(x'))) \
    & = 1/(4 pi) integral_0^(2 pi) dd(phi_k) integral_0^pi dd(theta_k) sin theta_k e^(i k abs(vb(x) - vb(x')) cos theta_k) \
    & = 1/2 integral_(-1)^(1) dd(cos theta_k) e^(i k abs(vb(x) - vb(x')) cos theta_k) \
    & = (sin(k abs(vb(x) - vb(x'))))/(k abs(vb(x) - vb(x')))
  $
  因此最终得到
  $
    sigma_"tot" & = (m^2)/(pi hbar^4) integral dd(vb(x), 3) V(vb(x)) integral dd(vb(x'), 3) V(vb(x')) (sin^2(k abs(vb(x) - vb(x'))))/(k^2 abs(vb(x) - vb(x'))^2)
  $
]

#exercise(subname: [Steck 7.4])[
  Remember that the trace is a sum over all diagonal matrix elements, which is true in any basis; in the context here the matrix elements run over a basis corresponding to a particular $j$, so that for example
  $
    Tr[hat(A)] = sum_m braket(j m, hat(A), j m)
  $
  for any operator $hat(A)$. Let $hat(J)_alpha$ denote the Cartesian components of the angular momentum vector $vb(J)$.
  + Show that $Tr[hat(J)_alpha] = 0$.
  + Show that
    $
      Tr[hat(J)_alpha hat(J)_beta] = (hbar^2 j(j+1)(2j+1))/3 delta_(alpha beta)
    $
  + Show that
    $
      Tr[hat(J)_alpha hat(J)_beta hat(J)_gamma] = (i hbar^3 j(j+1)(2j+1))/6 epsilon_(alpha beta gamma)
    $
  _Hint: use the commutation/anticommutation relations of the angular momentum operators and remember that the irreducible symmetric rank-2 tensor is proportional to the Kronecker delta $δ_(i j)$._
]

#proof[
  + 任取单位矢量$vu(n)$，对应无穷小旋转$hat(U) (theta) = e^(- i/hbar theta vb(hat(J)) dot vu(n))$，则有
    $
      Tr[hat(J)_alpha] & = Tr[hat(U) hat(J)_alpha hat(U)^dagger ] \
    $
    #proof[
      下面证明迹在酉变换下不变：
      $
        Tr[hat(U) hat(A) hat(U)^dagger ] & = sum_m braket(j m, hat(U) hat(A) hat(U)^dagger, j m) \
        & = sum_m sum_(m', m'') braket(j m, hat(U), j m') braket(j m', hat(A), j m'') braket(j m'', hat(U)^dagger, j m) \
        & = sum_(m', m'') delta_(m' m'') braket(j m', hat(A), j m'') \
        & = sum_m braket(j m, hat(A), j m) \
        & = Tr[hat(A)]
      $
    ]
    对$theta$取一阶导数，令$theta = 0$，得到
    $
      0 = eval(dv(, theta) Tr[hat(U) hat(J)_alpha hat(U)^dagger ])_(theta=0) = i/hbar Tr([vb(hat(J)) dot vu(n), hat(J)_alpha]) \
      = i/hbar Tr(i hbar epsilon_(beta alpha gamma) hat(J)_gamma vu(n)_beta) = - epsilon_(beta alpha gamma) vu(n)_beta Tr(hat(J)_gamma)
    $
    由于$vu(n)$任意，故$Tr(hat(J)_gamma) = 0$。
  + 考虑二阶张量
    $
      T_(alpha beta) = Tr[hat(J)_alpha hat(J)_beta]
    $
    在固定$j$的不可约表示中，$T_(alpha beta)$为对称张量，可以写成
    $
      T_(alpha beta) = A delta_(alpha beta)
    $
    取迹得到
    $
      Tr(T) = sum_alpha T_(alpha alpha) = 3 A = Tr(vb(hat(J)) dot vb(hat(J))) = (2j+1) hbar^2 j(j+1)
    $
    因此
    $
      Tr[hat(J)_alpha hat(J)_beta] = (hbar^2 j(j+1)(2j+1))/3 delta_(alpha beta)
    $
  + 考虑三阶张量
    $
      S_(alpha beta gamma) = Tr[hat(J)_alpha hat(J)_beta hat(J)_gamma]
    $
    在固定$j$的不可约表示中，$S_(alpha beta gamma)$的反对称部分可以写成
    $
      S_(alpha beta gamma) = B epsilon_(alpha beta gamma)
    $
    用对易/反对易分解：
    $
      hat(J)_alpha hat(J)_beta = 1/2 {hat(J)_alpha, hat(J)_beta} + 1/2 [hat(J)_alpha, hat(J)_beta] \
    $
    因此
    $
      S_(alpha beta gamma) & = 1/2 Tr[{hat(J)_alpha, hat(J)_beta} hat(J)_gamma] + 1/2 Tr([[hat(J)_alpha, hat(J)_beta] hat(J)_gamma]) \
    $
    注意到对称二阶张量只有Kronecker delta，因此
    $
      {hat(J)_alpha, hat(J)_beta} = C delta_(alpha beta)
    $
    从而
    $
      Tr[{hat(J)_alpha, hat(J)_beta} hat(J)_gamma] = C delta_(alpha beta) Tr(hat(J)_gamma) = 0
    $
    而
    $
      Tr([[hat(J)_alpha, hat(J)_beta] hat(J)_gamma]) & = i hbar epsilon_(alpha beta delta) Tr(hat(J)_delta hat(J)_gamma) \
      & = i hbar epsilon_(alpha beta delta) (hbar^2 j(j+1)(2j+1))/3 delta_(delta gamma) \
      & = i hbar^3 j(j+1)(2j+1)/3 epsilon_(alpha beta gamma)
    $
    因此
    $
      S_(alpha beta gamma) = (i hbar^3 j(j+1)(2j+1))/6 epsilon_(alpha beta gamma)
    $
]
