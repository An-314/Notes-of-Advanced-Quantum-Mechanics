#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第10次作业],
  author: "AnZrew",
  time: "2025年12月",
)

#exercise[
  For $hat(H) = vb(hat(P))^2/(2m) + V(vb(hat(X)))$, write down the standard path integral with $integral dd(vb(p))$ and $integral dd(vb(x))$. How to obtain an exactly equivalent path integral with $integral dd(vb(x))$ only?
]

#solution[

  传播子
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) = braket(vb(x)_B, exp(- i/hbar hat(H) (t_B - t_A)), vb(x)_A)
  $
  将时间间隔$[t_A, t_B]$等分为$N$份，每份为$epsilon = (t_A - t_B)/N$，则
  $
    e^(-i/hbar hat(H) (t_B - t_A)) = (e^(- i/hbar hat(H) epsilon))^N = e^(- i/hbar hat(H) epsilon) ... e^(- i/hbar hat(H) epsilon)
  $
  插入完备性关系
  $
    integral dd(x_i) ketbra(x_i) = hat(I), i = 1, 2, ..., N-1
  $
  得到
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) & = integral dd(vb(x)_(N-1)) ... integral dd(vb(x)_1) \
    & braket(vb(x)_B, e^(- i/hbar hat(H) epsilon), vb(x)_(N-1)) ... braket(vb(x)_1, e^(- i/hbar hat(H) epsilon), vb(x)_A) \
  $
  考虑$N$很大即$epsilon$很小的情形，考虑
  $
    e^(A+B) = e^A e^B e^(- 1/2 [A, B])
  $
  有
  $
    exp(- i/hbar hat(H) epsilon) & = exp(- (i epsilon)/hbar (hat(vb(p))^2/(2m) + V(hat(vb(x))))) \
    &= exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)) exp(- (i epsilon)/hbar V(hat(vb(x)))) + O(epsilon^2) \
  $
  因此
  $
    braket(vb(x)_n, e^(- i/hbar epsilon hat(H)), vb(x)_(n-1)) & tilde.eq braket(vb(x)_n, exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)) exp(- (i epsilon)/hbar V(hat(vb(x)))), vb(x)_(n-1)) \
    &= braket(vb(x)_n, exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)), vb(x)_(n-1)) exp(- (i epsilon)/hbar V(vb(x)_(n-1))) \
  $
  插入动量完备性关系
  $
    integral dd(vb(p)) ketbra(vb(p)) = hat(I)
  $
  得到
  $
    braket(vb(x)_n, exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)), vb(x)_(n-1)) &= integral dd(vb(p)_n) braket(vb(x)_n, vb(p)_n) braket(vb(p)_n, exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)), vb(p)_n) braket(vb(p)_n, vb(x)_(n-1)) \
  $
  再用平面波归一化
  $
    braket(vb(x), vb(p)) = 1/sqrt(2 pi hbar) exp(i/hbar vb(p) dot vb(x))
  $
  代入上式，得到
  $
    braket(vb(x)_n, exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)), vb(x)_(n-1)) & = integral dd(vb(p)_n) 1/(2 pi hbar) exp(i/hbar vb(p)_n dot (vb(x)_n - vb(x)_(n-1))- (i epsilon)/hbar vb(p)_n^2/(2m)) \
  $
  从而
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) & = integral (product_(n=1)^(N-1) dd(vb(x)_n)) (product_(n=1)^N dd(vb(p)_n)/(2 pi hbar)) \
    & exp((i)/hbar sum_(n=1)^N (vb(p)_n dot (vb(x)_n - vb(x)_(n-1)) - epsilon (vb(p)_n^2/(2m) + V(vb(x)_(n-1))))) \
  $
  指数的宗量可写成连续形式
  $
    exp((i)/hbar sum_(n=1)^N (vb(p)_n dot (vb(x)_n - vb(x)_(n-1)) - epsilon (vb(p)_n^2/(2m) + V(vb(x)_(n-1))))) \
    -> exp((i)/hbar integral_(t_A)^(t_B) dd(t) (vb(p)(t) dot dot(vb(x)) - H(vb(p), vb(x))))
  $
  从而得到路径积分形式
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) &= integral cal(D)(vb(p)) cal(D)(vb(x)) exp((i)/hbar integral_(t_A)^(t_B) dd(t) (vb(p)(t) dot dot(vb(x)) - H(vb(p), vb(x))))\
    & = integral cal(D)(vb(p)) cal(D)(vb(x)) exp((i)/hbar integral_(t_A)^(t_B) dd(t) (vb(p) dot dot(vb(x)) - (vb(p)^2)/(2m) - V(vb(x))))
  $
  在这里，$cal(D)(vb(p))$和$cal(D)(vb(x))$分别表示动量和坐标的路径积分测度。

  跃迁振幅中，可以将动量路径积分先进行计算
  $
    braket(vb(x)_n, exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)), vb(x)_(n-1)) & = integral dd(vb(p)_n) 1/(2 pi hbar) exp(i/hbar vb(p)_n dot (vb(x)_n - vb(x)_(n-1))- (i epsilon)/hbar vb(p)_n^2/(2m)) \
  $
  利用Gaussian积分公式
  $
    integral dd(p) exp(- a p^2 + b p) = sqrt(pi/a) exp(b^2/(4a)) (Re(a) > 0)
  $
  可以得到
  $
    braket(vb(x)_n, exp(- (i epsilon)/hbar hat(vb(p))^2/(2m)), vb(x)_(n-1)) & = (m/(2 pi hbar i epsilon))^(3/2) exp((i m)/(2 hbar epsilon) (vb(x)_n - vb(x)_(n-1))^2)
  $
  从而有
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) & = lim_(N -> oo) (m/(2 pi hbar i epsilon))^((3 N)/2) (product_(n=1)^(N-1) (m/(2 pi hbar i epsilon))^(3/2) dd(vb(x)_n))\ & exp((i epsilon)/hbar sum_(n=1)^N (m/2 ((vb(x)_n - vb(x)_(n-1))/epsilon)^2 - V(vb(x)_(n-1)))) \
  $
  指数的宗量可写成连续形式
  $
    exp((i epsilon)/hbar sum_(n=1)^N (m/2 ((vb(x)_n - vb(x)_(n-1))/epsilon)^2 - V(vb(x)_(n-1)))) \
    -> exp((i)/hbar integral_(t_A)^(t_B) dd(t) (m/2 dot(vb(x))^2 - V(vb(x))))
  $
  从而得到仅含坐标路径积分的形式
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) & = integral cal(D)(vb(x)) exp((i)/hbar integral_(t_A)^(t_B) dd(t) (m/2 dot(vb(x))^2 - V(vb(x))))
  $
]

#exercise[
  For a more general Hamiltonian, this cannot be done exactly. Usually the saddle point approximation is used, i.e. substitute for $vb(P) = (...)$ using the classical equation of motion (this is only exact in the case considered above). Try this for a relativistic particle $H = sqrt((vb(P)c)^2 + (m c^2)^2)$.
]

#solution[
  对一般$H(p,x)$，相空间离散形式仍是
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) & = integral (product_(n=1)^(N-1) dd(vb(x)_n)) (product_(n=1)^N dd(vb(p)_n)/(2 pi hbar)) \
    & exp((i)/hbar sum_(n=1)^N (vb(p)_n dot (vb(x)_n - vb(x)_(n-1)) - epsilon H(vb(p)_n, vb(x)_(n-1)))) \
  $
  $H$对$vb(p)$不是二次的，无法直接进行Gaussian积分。于是对每个时间片的$vb(p)_n$积分用鞍点近似。对
  $
    Phi(p_n) = vb(p)_n dot (vb(x)_n - vb(x)_(n-1)) - epsilon H(vb(p)_n, vb(x)_(n-1))
  $
  鞍点条件：
  $
    pdv(Phi, vb(p)_n) = 0 => pdv(H, vb(p)_n) = (vb(x)_n - vb(x)_(n-1))/epsilon
  $
  在鞍点$vb(p)_n = vb(p)_n^*$
  $
    Phi(vb(p)_n) & = vb(p)_n^* dot (vb(x)_n - vb(x)_(n-1)) - epsilon H(vb(p)_n^*, vb(x)_(n-1)) = epsilon (vb(p)_n^* dot dot(vb(x)) - H(vb(p)_n^*, vb(x)_(n-1))) \
  $
  这对应经典的Lagrangian$L = vb(p) dot dot(vb(x)) - H(vb(p), vb(x))$，因此鞍点近似后得到纯$vb(x)$的路径积分：
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) & approx integral cal(D)(vb(x)) exp((i)/hbar integral_(t_A)^(t_B) dd(t) L(vb(x), dot(vb(x)))) \
  $
  #newpara()
  对于相对论粒子
  $
    H = sqrt((vb(P)c)^2 + (m c^2)^2)
  $
  鞍点条件为
  $
    dot(vb(x)) = pdv(H, vb(p)) = (vb(p) c^2)/sqrt((vb(p)c)^2 + (m c^2)^2) = (vb(p) c^2)/H
  $
  考虑到
  $
    H = gamma m c^2, vb(p) = gamma m dot(vb(x)), gamma = 1/sqrt(1 - (dot(vb(x))/c)^2)
  $
  可解出
  $
    vb(p)^* = (m dot(vb(x)))/(sqrt(1 - (dot(vb(x))/c)^2))
  $
  因此
  $
    L(vb(x), dot(vb(x))) & = vb(p)^* dot dot(vb(x)) - H(vb(p)^*, vb(x)) \
                         & = - m c^2 sqrt(1 - (dot(vb(x))/c)^2)
  $
  写出鞍点近似后的纯$dd(vb(x))$路径积分为
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) & approx integral cal(D)(vb(x)) exp((i)/hbar integral_(t_A)^(t_B) dd(t) (- m c^2 sqrt(1 - (dot(vb(x))/c)^2))) \
  $
]
