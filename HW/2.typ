#import "@preview/scripst:1.1.1": *

#show: scripst.with(
  title: [高等量子力学第2次作业],
  author: "AnZrew",
  time: "2025年10月",
)

#exercise(subname: [Steck 5.13])[
  Recall that a coherent state can be defined by $ket(alpha) = D(alpha) ket(0)$, where $ket(0)$ is the harmonic-oscillator ground state, and in the case of purely imaginary $alpha$,
  $
    D(alpha) = e^(i sqrt(2) abs(alpha) x) , i alpha in RR
  $
  Also remember that the coherent state has the expansion
  $
    ket(alpha) = sum_(n=0)^oo (alpha^n)/(sqrt(n!)) e^(-abs(alpha)^2/2) ket(n)
  $
  in the energy basis.

  Use the above coherent-state expressions to derive an expression for the $n$th moment $braket(0, x^n, 0)$ of the   harmonic-oscillator ground state. You should not be setting up any integrals in this problem; work   entirely in Dirac notation.
]

#solution[
  取纯虚的参数$alpha = i y(y in RR)$，有
  $
    D(alpha) = exp(alpha a^dagger - alpha^* a) = exp(i sqrt(2) abs(alpha) x), x = (a^dagger + a)/sqrt(2)
  $
  $
    D(i y) = exp(i y (a^dagger + a)) = exp(i sqrt(2) y x)
  $
  相干态与基态的重叠由展开式给出
  $
    braket(0|alpha) = e^(- abs(alpha)^2/2)
  $
  于是
  $
    braket(0, D(i y), 0) = e^(- y^2/2)
  $
  考虑
  $
    chi(k) = braket(0, e^(i k x), 0) = braket(0, D(i k/sqrt(2)), 0) = e^(- k^2/4)
  $
  由矩生成关系
  $
    braket(0, x^n, 0) = 1/i^n eval(dv(, k, n) chi(k))_(k=0) = 1/i^n eval(dv(, k, n) e^(- k^2/4))_(k=0)
  $
  计算可得
  $
    braket(0, x^(2m+1), 0) = 0\
    braket(0, x^(2m), 0) = (2m - 1)!! / 2^m = (2m)! / (2^m m!)
  $
]

#exercise(subname: [Sakurai 2.7])[
  Consider a particle in three dimensions whose Hamiltonian is given by
  $
    hat(H) = hat(vb(p))^2/(2m) + hat(V)(vb(x))
  $
  By calculating $[x dot p, H]$, obtain
  $
    dv(, t) expval(vb(x) dot hat(vb(p))) = expval(hat(vb(p))^2/m - vb(x) dot grad hat(V)(vb(x)))
  $
  In order for us to identify the preceding relation with the quantum-mechanical analogue of the virial theorem, it is essential that the left-hand side vanish. Under what condition would this happen?
]

#solution[
  用正则对易关系
  $
    [ x_i, p_j ] = i hbar delta_(i j)
  $
  有
  $
    [hat(vb(x)) dot hat(vb(p)), hat(vb(p))^2/(2m)] &= 1/(2m) sum_k [hat(x)_k hat(p)_k, sum_j hat(p)_j^2]\
    &= 1/(2m) sum_(k, j) [hat(x)_k hat(p)_k, hat(p)_j^2]\
    &= 1/(2m) sum_(k, j) (hat(x)_k [hat(p)_k, hat(p)_j^2] + [hat(x)_k, hat(p)_j^2] hat(p)_k)\
    &= 1/(2m) sum_(k, j) (hat(x)_k 0 + 2 i hbar delta_(k j) hat(p)_j hat(p)_k)\
    &= (i hbar)/m sum_k hat(p)_k^2 = (i hbar)/m hat(vb(p))^2
  $
  $
    [hat(vb(x)) dot hat(vb(p)), hat(V)(hat(vb(x)))] &= sum_k [hat(x)_k hat(p)_k, hat(V)(hat(vb(x)))]\
    &= sum_k (hat(x)_k [hat(p)_k, hat(V)(hat(vb(x)))] + [hat(x)_k, hat(V)(hat(vb(x)))] hat(p)_k)\
    &= sum_k (hat(x)_k (- i hbar (pdv(, x_k) hat(V)(hat(vb(x)))) + 0)\
    &= - i hbar sum_k hat(x)_k (pdv(, x_k) hat(V)(hat(vb(x)))) = - i hbar vb(x) dot grad hat(V)(vb(x))
  $
  因此
  $
    [hat(vb(x)) dot hat(vb(p)), hat(H)] = (i hbar)/m hat(vb(p))^2 - i hbar vb(x) dot grad hat(V)(vb(x))
  $
  代入Heisenberg方程
  $
    dv(, t) expval(vb(x) dot hat(vb(p))) = 1/(i hbar)expval([hat(vb(x)) dot hat(vb(p)), hat(H)]) = expval(hat(vb(p))^2/m - vb(x) dot grad hat(V)(vb(x)))
  $
  得证。要使左端为零，必须要求
  $
    dv(, t) expval(vb(x) dot hat(vb(p))) = 0
  $
  这在以下情形成立：
  - 定态：若体系处于能量本征态$ket(E)$则任意不显含时的算符的期望值不随时间变化。
  - 有界运动的长时平均：对束缚态$expval(hat(vb(x)) dot hat(vb(p)))$有界且不随时间发散，取长时间平均则其导数为零。
]

#exercise(subname: [Sakurai 2.12], lab: "Sakurai-2.12")[
  Consider a particle subject to a one-dimensional simple harmonic oscillator potential. Suppose that at $t = 0$ the state vector is given by
  $
    exp((-i hat(p) a)/hbar) ket(0)
  $
  where $p$  is the momentum operator and a is some number with dimension of length. Using the Heisenberg picture, evaluate the expectation value $expval(x)$  for $t >= 0$.
]

#solution[
  简谐振子的Hamilton量为
  $
    hat(H) = hat(p)^2/(2m) + 1/2 m omega^2 hat(x)^2
  $
  在Heisenberg绘景下，位置算符的演化为
  $
    hat(x)^"H" (t) = hat(x) cos(omega t) + hat(p)/(m omega) sin(omega t)
  $
  初态是平移算符$hat(T)(a) = exp((- i hat(p) a)/hbar)$作用在基态$ket(0)$上得到的，即
  $
    hat(T)^dagger (a) hat(x) hat(T)(a) = hat(x) + a\
    hat(T)^dagger (a) hat(p) hat(T)(a) = hat(p)
  $
  因此
  $
    expval(x)(t) & = braket(0, hat(T)^dagger (a) hat(x)^"H" (t) hat(T)(a), 0) \
                 & = braket(0, (hat(x) + a) cos(omega t) + hat(p)/(m omega) sin(omega t), 0) \
                 & = braket(0, hat(x)^"H" (x), 0) + a cos(omega t) \
                 & = a cos(omega t)
  $
]

#exercise(subname: [Sakurai 2.13])[
  + Write down the wave function (in coordinate space) for the state specified in @Sakurai-2.12 at $t=0$. You may use
    $
      braket(x', 0) = pi^(-1/4) x_0^(-1/2) exp(-1/2 (x'/x_0)^2) , x_0 = (hbar/(m omega))^(1/2)
    $
  + Obtain a simple expression for the probability that the state is found in the ground state at $t = 0$. Does this probability change for $t > 0$?
]

#solution[
  初态为
  $
    ket(psi(0)) = e^((- i hat(p) a)/hbar) ket(0)
  $
  + $t=0$的坐标表象波函数为
    $
      psi(x', 0) & = braket(x', psi(0)) = braket(x', e^((- i hat(p) a)/hbar), 0) = braket(x' - a, 0) \
                 & = pi^(-1/4) x_0^(-1/2) exp(-1/2 ((x' - a)/x_0)^2)
    $
  + 在$t=0$处于基态的概率为
    $
      P_0 (0) &= abs(braket(0, psi(0)))^2 = abs(braket(0, e^((- i hat(p) a)/hbar), 0))^2\
      &= abs(braket(0, hat(T)(a), 0))^2 = abs(integral_( - oo)^(oo) braket(0, x') braket(x', hat(T)(a), 0) dd(x'))^2\
      &= abs(integral_( - oo)^(oo) braket(0, x') braket(x' - a, 0) dd(x'))^2\
      &= abs(integral_( - oo)^(oo) (pi^(-1/4) x_0^(-1/2) exp(-1/2 (x'/x_0)^2)) (pi^(-1/4) x_0^(-1/2) exp(-1/2 ((x' - a)/x_0)^2)) dd(x'))^2\
      &= abs(pi^(-1/2) x_0^(-1) integral_( - oo)^(oo) exp(-1/2 (x'/x_0)^2 - 1/2 ((x' - a)/x_0)^2) dd(x'))^2\
      &= abs(pi^(-1/2) x_0^(-1) integral_( - oo)^(oo) exp(- (2x'^2 - 2 a x' + a^2)/(2 x_0^2)) dd(x'))^2\
      &= abs(pi^(-1/2) x_0^(-1) e^(- a^2/(4 x_0^2)) integral_( - oo)^(oo) exp(- (x' - a/2)^2/(x_0^2)) dd(x'))^2\
      &= exp(- a^2/(2 x_0^2))
    $
    而对$t>0$
    $
      P_0 (t) = abs(braket(0, psi(t)))^2 = abs(braket(0, e^((- i hat(p) a)/hbar) 0))^2 = abs(e^((- i E_0 t)/hbar) braket(0, e^((- i hat(p) a)/hbar), 0))^2 = P_0 (0)
    $
    因此
    $
      P_0 (t) = exp(- a^2/(2 x_0^2))
    $
]
