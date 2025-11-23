#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第6次作业],
  author: "AnZrew",
  time: "2025年11月",
)

#exercise(subname: [Sakurai 4.7])[
  + Let $psi(x, t)$ be the wave function of a spinless particle corresponding to a plane wave in three dimensions. Show that $psi^*(x, -t)$ is the wave function for the plane wave with the momentum direction reversed.
  + Let $chi(vu(n))$ be the two-component eigenspinor of $vb(sigma) dot vu(n)$ with eigenvalue $+1$. _Using that the explicit characterize form_ of $chi(vu(n))$ (in terms of the polar and azimuthal angles $beta$ and $gamma$ that characterize $vu(n)$), verify that $- i sigma_2 chi^*(vu(n))$ is the two-component eigenspinor with the spin direction reversed.
]

#proof[

  + 三维平面波可写为
    $
      psi(vb(x), t) prop e^(i/hbar (vb(p) dot vb(x) - E t))
    $
    则
    $
      psi^*(vb(x), -t) prop e^(- i/hbar (vb(p) dot vb(x) + E t)) = e^(i/hbar ((- vb(p)) dot vb(x) - E t))
    $
    即动量方向反转的平面波。

  + 考虑
    $
      vu(u) = mat(sin beta cos gamma, sin beta sin gamma, cos beta)
    $
    $vb(sigma) dot vu(n)$的$+1$本征自旋子可取$sigma_z$基下的形式为
    $
      chi_+ (vu(n)) = mat(cos beta/2; e^(i gamma) sin beta/2)\
    $
    则
    $
      chi_+^* (vu(n)) = mat(cos beta/2; e^(- i gamma) sin beta/2)\
      i sigma_2 = mat(0, -1; 1, 0)\
    $
    直接作用可得
    $
      - i sigma_2 chi_+^* (vu(n)) = mat(- e^(- i gamma) sin beta/2; cos beta/2) = chi_- (vu(n))
    $
    即动量方向反转的本征自旋子。
]

#exercise(subname: [Sakurai 4.11])[
  Suppose a spinless particle is bound to a fixed center by a potential $V(vb(x))$ so asymmetrical that no energy level is degenerate. Using time-reversal invariance, prove
  $
    expval(vb(L)) = 0
  $
  for any energy eigenstate. (This is known as quenching of orbital angular momentum.) If the wave function of such a nondegenerate eigenstate is expanded as
  $
    sum_l sum_m F_(l m)(vb(r)) Y_l^m (theta, phi)
  $
  what kind of phase restrictions do we obtain on $F_(l m)(vb(r))$?
]

#proof[
  对无自旋粒子，时间反演算符$hat(P)$是反幺正的，且在坐标表象下等同于复共轭
  $
    hat(P) psi(vb(x)) = psi^*(vb(x))
  $
  它满足
  $
    hat(P) i hat(P)^dagger = - i\
    hat(P) hat(vb(r)) hat(P)^dagger = hat(vb(r))\
    hat(P) hat(vb(p)) hat(P)^dagger = - hat(vb(p))
  $
  因而
  $
    hat(P) hat(vb(L)) hat(P)^dagger = hat(P) (hat(vb(r)) times hat(vb(p))) hat(P)^dagger = hat(vb(r)) times (- hat(vb(p))) = - hat(vb(L))
  $
  若势能实（时间反演不变），则
  $
    [hat(H), hat(P)] = 0
  $
  对非简并能量本征态$ket(n)$有
  $
    hat(P) ket(n) = e^(i alpha) ket(n)
  $
  用反幺正变换的矩阵元恒等式
  $
    braket(phi, hat(A), psi) = braket(hat(P) psi, hat(P) hat(A) hat(P)^dagger, hat(P) phi)^*
  $
  有
  $
    braket(n, hat(vb(L)), n) = braket(hat(P) n, hat(P) hat(vb(L)) hat(P)^dagger, hat(P) n)^* = braket(e^(i alpha) n, - hat(vb(L)), e^(i alpha) n)^* = - braket(n, hat(vb(L)), n)^*
  $
  因为矩阵元为实数，故
  $
    expval(vb(L)) = braket(n, hat(vb(L)), n) = 0
  $
  #newpara()

  时间反演（对无自旋）在坐标表象等同于复共轭。由的非简并性，可选取相位使本征态满足
  $
    hat(P) ket(n) = ket(n) => psi_n^*(vb(r)) = psi_n (vb(r))
  $
  即波函数可整体选为实函数。球谐函数满足
  $
    Y_l^m (theta, phi)^* = (-1)^m Y_l^(- m) (theta, phi)
  $
  对于
  $
    psi(vb(r)) = sum_l sum_m F_(l m)(vb(r)) Y_l^m (theta, phi)
  $
  取共轭并令$psi^*(vb(r)) = psi(vb(r))$，得到
  $
    sum_l sum_m F_(l m)^* (vb(r)) (-1)^m Y_l^(- m) (theta, phi) = sum_l sum_m F_(l m)(vb(r)) Y_l^m (theta, phi)
  $
  把求和指标$m->-m$后逐项比较，得到系数的相位约束：
  $
    F_(l - m)^* (vb(r)) = (-1)^m F_(l m)(vb(r))
  $
]

#exercise(subname: [Sakurai 4.12])[
  The Hamiltonian for a spin 1 system is given by
  $
    H = A S_z^2 + B (S_x^2 - S_y^2)
  $
  Solve this problem exactly to find the normalized energy eigenstates and eigenvalues. (A spin-dependent Hamiltonian of this kind actually appears in crystal physics.) Is this Hamiltonian invariant under time reversal? How do the normalized eigenstates you obtained transform under time reversal?
]

#solution[
  利用梯算符
  $
    hat(S)_plus.minus = hat(S)_x plus.minus hat(S)_y
  $
  有
  $
    hat(S)_x^2 - hat(S)_y^2 = 1/2 (hat(S)_plus^2 + hat(S)_minus^2)
  $
  在$hat(S)_z$本征态基底$ket(1)$，$ket(0)$，$ket(-1)$下
  $
    hat(S)_z ket(m) = m hbar ket(m)\
  $
  $hat(S)_plus.minus^2$的作用为
  $
    hat(S)_minus^2 ket(1) = 2 hbar^2 ket(-1)\
    hat(S)_plus^2 ket(-1) = 2 hbar^2 ket(1)\
    hat(S)_plus.minus^2 ket(0) = 0
  $
  从而Hamiltonian矩阵表示为
  $
    hat(H) = mat(A hbar^2, 0, B hbar^2; 0, 0, 0; B hbar^2, 0, A hbar^2)
  $
  求解特征值方程
  $
    det(hat(H) - E hat(I)) = 0
  $
  得本征值
  $
    E_+ & = A hbar^2 + B hbar^2, \
    E_0 & = 0, \
    E_- & = A hbar^2 - B hbar^2
  $
  对应的归一化本征态为
  $
    ket(psi_+) & = 1/sqrt(2) (ket(1) + ket(-1)) \
    ket(psi_0) & = ket(0) \
    ket(psi_-) & = 1/sqrt(2) (ket(1) - ket(-1))
  $
  #newpara()
  下面说明时间反演不变性。自旋体系的时间反演算符$hat(P)$满足
  $
    hat(P) vb(S) hat(P)^dagger = - vb(S), hat(P) i hat(P)^dagger = - i
  $
  因此
  $
    hat(P) hat(S)_x^2 hat(P)^dagger = hat(S)_x^2\
    hat(P) hat(S)_y^2 hat(P)^dagger = hat(S)_y^2\
    hat(P) hat(S)_z^2 hat(P)^dagger = hat(S)_z^2
  $
  进而
  $
    hat(P) hat(H) hat(P)^dagger = A hat(S)_z^2 + B (hat(S)_x^2 - hat(S)_y^2) = hat(H)
  $
  即Hamiltonian在时间反演下不变。

  本征态在时间反演下的变换
  $
    hat(P) ket(m) = (-1)^(1 - m) ket(- m)
  $
  对本征态：
  $
    hat(P) ket(psi_+) & = 1/sqrt(2) (hat(P) ket(1) + hat(P) ket(-1)) = 1/sqrt(2) ( - ket(-1) + ket(1)) = ket(psi_+) \
    hat(P) ket(psi_0) & = hat(P) ket(0) = - ket(0) = - ket(psi_0) \
    hat(P) ket(psi_-) & = 1/sqrt(2) (hat(P) ket(1) - hat(P) ket(-1)) = 1/sqrt(2) ( - ket(-1) - ket(1)) = - ket(psi_-)
  $
]
