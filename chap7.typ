#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

= 路径积分量子化

== 正则量子化

所谓量子化(quantization)，就是将一个经典理论转化为的量子理论的手续。

目前我们学习的量子力学都是基于*正则量子化*手续。假设经典系统的自由度为一系列的正则坐标 ${q_i}$，对应的正则动量为${p_i}$，它们的时间演化满足Hamilton正则方程
$
  dot(q)_i = pdv(cal(H), p_i), dot(p)_i = -pdv(cal(H), q_i)
$
其中Hamilton量$cal(H)$是$q_i, p_i, t$ 的函数，$H = H(q_i, p_i, t)$。任意力学量$cal(A) = cal(A)(q_i, p_i, t)$的演化方程为
$
  dv(cal(A), t) & = sum_i (pdv(cal(A), q_i) dot(q)_i + pdv(cal(A), p_i) dot(p)_i) + pdv(cal(A), t) \
                & = sum_i (pdv(cal(A), q_i) pdv(cal(H), p_i) - pdv(cal(A), p_i) pdv(cal(H), q_i)) + pdv(cal(A), t)
$
定义任意两个力学量$cal(A)(q_i, p_i, t)$和$cal(B)(q_i, p_i, t)$的Poisson括号为
$
  {cal(A), cal(B)} = sum_i (pdv(cal(A), q_i) pdv(cal(B), p_i) - pdv(cal(A), p_i) pdv(cal(B), q_i))
$
则演化方程可写为
$
  dv(cal(A), t) = {cal(A), cal(H)} + pdv(cal(A), t)
$
若力学量$cal(A)$不显含时间，则
$
  dv(cal(A), t) = {cal(A), cal(H)}
$
尤其是Hamilton正则方程可以写为
$
  dot(q)_i = {q_i, cal(H)}, dot(p)_i = {p_i, cal(H)}
$

#newpara()

*正则量子化*：在经典理论中，基本的Poisson括号为
$
  {q_i, q_j} = 0, {p_i, p_j} = 0, {q_i, p_j} = delta_(i j)
$
所谓正则量子化，就是将经典理论中的*正则坐标和正则动量*替换为相应的*量子力学算符*，
$
  q_i -> hat(q)_i, p_i -> hat(p)_i
$
*基本Poisson括号*替换为算符的*基本对易关系*
$
  [hat(q)_i, hat(q)_j] = 0, [hat(p)_i, hat(p)_j] = 0, [hat(q)_i, hat(p)_j] = i hbar delta_(i j)
$
这组对易关系称为正则对易关系，它是正则量子化的核心。

在Schrödinger绘景中，算符$hat(q)_i$和$hat(p)_i$都不随时间变化，无法与经典理论对应。若采用*Heisenberg绘景*，算符$hat(q)_i$和$hat(p)_i$随时间变化，*基本对易关系*变为*等时对易关系论*
$
  [hat(q)_i (t), hat(q)_j (t)] = 0, [hat(p)_i (t), hat(p)_j (t)] = 0, [hat(q)_i (t), hat(p)_j (t)] = i hbar delta_(i j)
$
任何不含时的力学量$hat(A)$的演化方程为Heisenberg方程
$
  dv(hat(A), t) = 1/(i hbar) [hat(A), hat(H)]
$
此式与经典力学的演化方程对应，而且可以看到，量子力学的对易关系和经典力学的Poisson括号有对应关系
$
  1/(i hbar) [hat(A), hat(B)] <-> {cal(A), cal(B)}
$

#newpara()

用正则量子化方法来构建带电粒子在电磁场中的量子理论。我们需要量子化的是带电粒子的运动，电磁场还是经典的。在经典力学中，质量为$m$、电荷为$q$的带电粒子在电磁场中的运动方程为(采用国际单位制)
$
  m dot.double(x) = q(vb(E) + dot(vb(x)) times vb(B))
$
由此可以构造出Lagrange量
$
  L = 1/2 m dot(vb(x))^2 + q dot(vb(x)) dot vb(A) - q phi
$
其中标量势$phi(vb(x), t)$和矢量势$vb(A)(vb(x),t)$定义为
$
  vb(E) = - grad phi - pdv(vb(A), t), vb(B) = grad times vb(A)
$
对于任意函数$Lambda(vb(x), t)$，若标量势和矢量势做如下变换
$
  phi -> phi - pdv(Lambda, t), vb(A) -> vb(A) + grad Lambda
$
电场$vb(E)$和磁场$vb(B)$是不变的，此即规范变换。在此变换下，Lagrange量变为
$
  L & -> L + q(dot(x) dot grad Lambda + pdv(Lambda, t)) \
    & = L + q (sum_i pdv(x_i, t) pdv(Lambda, x_i) + pdv(Lambda, t)) \
    & = L + q dv(Lambda, t)
$
因而*作用量$S = integral dd(t) L$只变化一个表面项*。

正则动量为
$
  p_i = pdv(L, x_i) = m dot(x)_i + q A_i => vb(p) = m dot(vb(x)) + q vb(A)
$
带电粒子的Hamilton量为
$
  cal(H) = sum_i p_i dot(x)_i - L = (vb(p) - q vb(A))/(2 m) + q phi
$
将$x_i$和$p_i$替换为算符$hat(x)_i$和$hat(p)_i$，满足正则对易关系
$
  [hat(x)_i, hat(x)_j] = 0, [hat(p)_i, hat(p)_j] = 0, [hat(x)_i, hat(p)_j] = i hbar delta_(i j)
$
则量子力学的Hamilton算符为
$
  hat(H) = (hat(vb(p)) - q vb(A)(hat(vb(x)), t))^2/(2 m) + q phi(hat(vb(x)), t)
$
这其中存在的一个问题是：由于标量势和矢量势依赖空间坐标$vb(x)$，因此$hat(vb(p))$与$vb(A)(hat(vb(x)), t)$不对易，所以Hamilton量中的$(hat(vb(p)) - q vb(A))^2$需要进一步明确。常用的方案是
$
  (hat(vb(p)) - q vb(A))^2 & = (hat(vb(p)) - q vb(A)) dot (hat(vb(p)) - q vb(A)) \
                           & = hat(vb(p))^2 -q(hat(vb(p)) dot vb(A) + vb(A) dot hat(vb(p))) + q^2 vb(A)^2
$
此即*对称化排序方案*，它自动保证Hamilton量是厄米算符。

*Heisenberg绘景*：利用Heisenberg方程，可计算坐标的演化方程
$
  dv(hat(x)_i, t) & = 1/(i hbar) [hat(x)_i, hat(H)] \
                  & = 1/(i hbar) [hat(x)_i, (hat(vb(p)) - q vb(A)(hat(vb(x)), t))^2/(2 m)] \
                  & = 1/(2 m i hbar) [hat(x)_i, sum_j (hat(p)_j - q A_j)(hat(p)_j - q A_j)] \
                  & = (hat(p)_i - q A_i)/m
$
进一步计算得到
$
  m dv(hat(x)_i, t, 2) & = dv(hat(p)_i - q A_i, t) \
  & = 1/(i hbar) [hat(p)_i - q A_i, hat(H)] \
  & = q/(2m) sum_(j k) epsilon_(i j k) (hat(p)_j - q A_j) B_k - q/(2m) sum_(j k) epsilon_(i j k) B_j (hat(p)_k - q A_k) - q pdv(A_0, x_i)
$
写成矢量形式，即
$
       dv(hat(vb(x)), t) & = (hat(vb(p)) - q vb(A))/m \
  m dv(hat(vb(x)), t, 2) & = q (vb(E) + 1/2 (dv(hat(vb(x)), t) times vb(B) + vb(B) times dv(hat(vb(x)), t)))
$
#newpara()
*Schrödinger绘景*：Schrödinger方程为
$
  i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t)) = ((hat(vb(p)) - q vb(A)(hat(vb(x)), t))^2/(2 m) + q phi(hat(vb(x)), t)) ket(psi(t))
$
这个方程明*显依赖于标量势和矢量势*。因此产生的一个问题就是：对于同样的电场和磁场，如果选取不同的标量势和矢量势，得到的量子态或者波函数必然是不同的，那么这些*不同的选择会不会改变物理可观测量*？比如，考虑沿$z$方向的均匀磁场，$vb(B) = B vu(z)$，矢量势可以选取*Landau规范*
$
  A_x = - B y, A_y = 0, A_z = 0
$
也可以选择*对称规范*
$
  A_x = - 1/2 B y, A_y = 1/2 B x, A_z = 0
$
考虑不随时间变化的磁场。可以证明，在矢量势的*规范变换*
$
  vb(A) -> tilde(vb(A)) = vb(A) + grad Lambda(vb(x))
$
下，*Schrödinger方程亦是规范不变的*。定义幺正算符
$
  hat(cal(G)) = exp((i q Lambda(hat(vb(x))))/hbar)
$
直接计算可以证明如下变换关系
$
                     hat(cal(G))^dagger hat(vb(x)) hat(cal(G)) & = hat(vb(x)) \
  hat(cal(G))^dagger (hat(vb(p)) - q tilde(vb(A))) hat(cal(G)) & = hat(vb(p)) - q vb(A)
$
定义量子态的变换
$
  ket(tilde(psi)(t)) = hat(cal(G))^dagger ket(psi(t)) = exp((i q Lambda(hat(vb(x))))/hbar) ket(psi(t))
$
则可以推导出变换后的量子态满足的演化方程为
$
  i hbar pdv(, t) ket(tilde(psi)(t)) = ((hat(vb(p)) - q tilde(vb(A))(hat(vb(x)), t))^2/(2 m) + q phi(hat(vb(x)), t)) ket(tilde(psi)(t))
$
在坐标表象下，这意味着两种不同的规范选取对应的波函数只差一个整体相位因子，即
$
  tilde(psi)(vb(x), t) = exp((i q Lambda(vb(x)))/hbar) psi(vb(x), t)
$

== 路径积分

正则量子化方案与经典力学的Hamilton形式对应，因此将Hamilton量放到了一个核心的位置。所以一个很自然的问题就是：可以不可以从经典系统的Lagrange量出发来进行量子化呢？

本小节将学习另一种量子化方案：*路径积分量子化*。

=== 传播子与Green函数

事实上，*路径积分量子化就是利用一种泛函积分的方法计算量子力学中的传播子(propagator)*。

讨论路径积分，先引入传播子的概念。在Schrödinger绘景中，量子力学体系在$t_1$时刻和$t_2$时刻的状态$(t_1 < t_2)$由时间演化算符联系起来
$
  ket(psi(t_2)) = hat(U)(t_2, t_1) ket(psi(t_1))
$
进入坐标表象，得到
$
  braket(vb(x)_2, psi(t_2)) & = braket(vb(x)_2, hat(U)(t_2, t_1), psi(t_1)) \
  & = integral dd(vb(x)_1, 3) braket(vb(x)_2, hat(U)(t_2, t_1), vb(x)_1) braket(vb(x)_1, psi(t_1))
$
定义*传播子*为
$
  K(vb(x)_2, t_2; vb(x)_1, t_1) = braket(vb(x)_2, hat(U)(t_2, t_1), vb(x)_1)
$
则可以写为(*物质波的Huygens原理*)
$
  psi(vb(x)_2, t_2) = integral dd(vb(x)_1, 3) K(vb(x)_2, t_2; vb(x)_1, t_1) psi(vb(x)_1, t_1)
$
假设$t_1$时刻粒子处于位置本征态，$psi(vb(x)_1, t_1) = delta(vb(x)_1 - vb(x)_0)$，根据上式立刻得到$psi(vb(x)_2, t_2) = K(vb(x)_2 , t_2 ; vb(x)_0 , t_1)$。为了方便，不妨仍将$vb(x)_0$记为$vb(x)_1$。那么传播子$K(vb(x)_2, t_2 ; vb(x)_1 , t_1)$的物理意义是：
- 若粒子在$t_1$时刻处于空间$vb(x)_1$处(位置本征态)，则传播子表示在以后任意$t_2$时刻($t_2 > t_1$)粒子处于空间点$vb(x)_2$处的*概率幅*。简而言之，粒子从时空点$(t_1, vb(x)_1)$传播到$(t_2, vb(x)_2)$的概率振幅。
- 若粒子在$t_1$时刻并不处于位置本征态，则利用传播子和$t_1$时刻的波函数可以得到以后任意$t_2$时刻的波函数。这实际上就是*符合因果律的幺正时间演化在坐标表象的体现*。

若Hamilton量$hat(H)$不显含时间，则可利用
$
  hat(H) ket(n) = E_n ket(n), hat(U)(t_2, t_1) = exp(- i/hbar hat(H) (t_2 - t_1))
$
得到
$
  K(vb(x)_2, t_2; vb(x)_1, t_1) &= braket(vb(x)_2, exp(- i/hbar hat(H) (t_2 - t_1)), vb(x)_1) \
  & = sum_(n n') braket(vb(x)_2, ket(n)) braket(n, exp(- i/hbar hat(H) (t_2 - t_1)), n') braket(n', vb(x)_1) \
  & = sum_n psi_n^*(vb(x)_1) psi_n (vb(x)_2) exp(- i/hbar E_n (t_2 - t_1))
$
其中$psi_n(vb(x)) = braket(vb(x), ket(n))$为坐标表象下Hamilton量的本征函数。当$t_2 -> t_1$时
$
  lim_(t_2 -> t_1) K(vb(x)_2, t_2; vb(x)_1, t_1) = sum_n psi_n^*(vb(x)_1) psi_n (vb(x)_2) = delta(vb(x)_2 - vb(x)_1)
$
当$t > t_1$时，易证传播子$K(vb(x), t; vb(x)_1, t_1)$满足Schrödinger方程
$
  i hbar pdv(, t) K(vb(x), t; vb(x)_1, t_1) = (- hbar^2/(2m) laplacian + V(vb(x))) K(vb(x), t; vb(x)_1, t_1)
$
利用传播子的定义以及时间演化算符$hat(U)(t; t_1)$满足的Schrödinger方程即可证明之。对于$t < t_1$情形，传播子无定义。若要拓展其定义，根据因果律可定义
$
  K(vb(x), t; vb(x)_1, t_1) = 0, (t < t_1)
$
此时亦满足Schrödinger方程。但$t = t_1$时情形未知，传播子会出现不连续。根据
$
  lim_(t -> t_1) K(vb(x), t; vb(x)_1, t_1) = delta(vb(x) - vb(x)_1)
$
可以证明传播子满足如下微分方程
$
  (i hbar pdv(, t) + hbar^2/(2m) laplacian - V(vb(x))) K(vb(x), t; vb(x)_1, t_1) = i hbar delta(t - t_1) delta(vb(x) - vb(x)_1)
$
#proof[
  将上式对$t$积分，积分区间为$(- oo, t_1 + epsilon]$，$epsilon = 0^+$，左边
  $
    i hbar eval(K(vb(x), t; vb(x)_1, t_1))_(-oo)^(t_1 + epsilon] = i hbar lim_(t -> t_1) K(vb(x), t; vb(x)_1, t_1) = i hbar delta(vb(x) - vb(x)_1)
  $
  右边的积分也正好等于这个结果。
]
这事实上是说，*传播子是含有源项的Schrödinger方程的Green函数*。

方程右边可以理解为某种“点源”的影响，即：传播子$K(vb(x), t; vb(x)_1; t_1)$是含时Schrödinger方程的一类*Green函数*，即所谓的*推迟Green函数*。所谓推迟，就是要求$t < t_1$ 时 $K(vb(x), t; vb(x)_1; t_1) = 0$。所以传播子的定义可以写为以传播子的定义可以写为
$
  K(vb(x), t; vb(x)_1, t_1) = braket(vb(x), hat(U)(t, t_1), vb(x)_1) theta(t - t_1)
$
其中$theta(x)$为阶跃函数。利用此式可以更直接地证明刚才的微分方程。

在理论物理中，*“传播子” (propagator)这个词通常与“Green函数”等价*，而且存在各种各样的传播子(Green函数)，推迟Green函数只是其中的一种。

#example(subname: [自由粒子])[
  考虑一维自由粒子，Hamilton量为$hat(H) = p^2/2m$，其能量本征态可取为动量本征态。传播子计算如下：
  $
    K(x_2, t_2; x_1, t_1) & = braket(x_2, exp(- i/hbar hat(H) (t_2 - t_1)), x_1) \
                          & = integral dd(p) braket(x_2, exp(- i/hbar hat(H) (t_2 - t_1)), p) braket(p, x_1) \
                          & = integral dd(p) braket(x_2, p) braket(p, x_1) exp(- i/hbar p^2/(2m) (t_2 - t_1)) \
  $
  利用
  $
    braket(x, p) = 1/sqrt(2 pi hbar) exp(i/hbar p x)
  $
  得到
  $
    K(x_2, t_2; x_1, t_1) & = 1/(2 pi hbar) integral_(-oo)^oo dd(p) exp(i/hbar (p (x_2 - x_1) - p^2/(2m) (t_2 - t_1))) \
  $
  利用 Fresnel 积分
  $
    integral_(-oo)^oo dd(x) e^(i a x^2) = sqrt((2 pi i)/a)
  $
  得到
  $
    K(x_2, t_2; x_1, t_1) = sqrt(m/(2 pi i hbar (t_2 - t_1))) exp((i m (x_2 - x_1)^2)/(2 hbar (t_2 - t_1)))
  $
  三维自由粒子也可做同样的计算，结果为
  $
    K(vb(x)_2, t_2; vb(x)_1, t_1) = (m/(2 pi i hbar (t_2 - t_1)))^(3/2) exp((i m (vb(x)_2 - vb(x)_1)^2)/(2 hbar (t_2 - t_1)))
  $
  可以证明，指数因子可以写成$exp(i/hbar S_"cl")$其中$S_"cl"$是经典运动对应的作用量
  $
    S_"cl" = integral_(t_1)^(t_2) dd(t) L = integral_(t_1)^(t_2) dd(t) (1/2 m dot(x)^2) = (m (x_2 - x_1)^2)/(2 (t_2 - t_1))
  $
]

#example(subname: [谐振子])[
  考虑一维谐振子
  $
    hat(H) ket(n) = E_n ket(n), E_n = hbar omega (n + 1/2)
  $
  传播子计算如下：
  $
    K(x_2, t_2; x_1, t_1) & = braket(x_2, exp(- i/hbar hat(H) (t_2 - t_1)), x_1) \
                          & = sum_n braket(x_2, exp(- i/hbar hat(H) (t_2 - t_1)), n) braket(n, x_1) \
                          & = sum_n braket(x_2, n) braket(n, x_1) exp(- i/hbar E_n (t_2 - t_1)) \
                          & = sum_n psi_n^*(x_1) psi_n (x_2) exp(- i/hbar E_n (t_2 - t_1)) \
  $
  利用谐振子本征波函数的表达式
  $
    psi_n(x) = 1/sqrt(2^n n!) ((m omega)/(pi hbar))^(1/4) H_n (xi) exp(- 1/2 xi^2), xi = sqrt((m omega)/hbar) x
  $
  代入得到
  $
    K(x_2, t_2; x_1, t_1) = \
    ((m omega)/(pi hbar))^(1/2) e^(- 1/2 (xi_1^2 + xi_2^2) - i/2 omega (t_2 - t_1)) sum_n 1/(2^n n!) H_n (xi_1) H_n (xi_2) e^(- i omega n (t_2 - t_1))\
  $
  利用Hermite多项式的积分表示
  $
    H_n (xi) = 1/sqrt(pi) (2/i)^n e^(xi^2) integral_(-oo)^(oo) dd(tau) e^(- tau^2 + 2 i xi tau) tau^n
  $
  代入得到
  $
    &K(x_2, t_2; x_1, t_1) \
    =& ((m omega)/(pi^3 hbar))^(1/2) e^(- 1/2 (xi_1^2 + xi_2^2) - i/2 omega (t_2 - t_1)) integral_(-oo)^(oo) dd(tau) integral_(-oo)^(oo) dd(tau') e^(- (tau^2 + tau'^2) + 2 i (xi_2 tau + xi_1 tau')) sum_(n=0)^(oo) (-1)^n/n! (2 tau tau')^n e^(- i n omega (t_2 - t_1)) \
    = & ((m omega)/(pi^3 hbar))^(1/2) e^(- 1/2 (xi_1^2 + xi_2^2) - i/2 omega (t_2 - t_1)) integral_(-oo)^(oo) dd(tau) e^(- tau^2 + 2 i xi_2 tau) integral_(-oo)^(oo) dd(tau') exp(- tau'^2 + 2 i xi_1 tau' - 2 tau tau' e^(- i omega (t_2 - t_1))) \
  $
  对$tau, tau'$的积分可利用Gauss积分
  $
    integral_(-oo)^(oo) dd(x) e^(- a x^2 + b x) = sqrt(pi/a) e^(b^2/(4 a)) (Re(a) > 0)
  $
  最终计算得到
  $
    &k(x_2, t_2; x_1, t_1) \
    = &sqrt((m omega)/(2 pi i hbar sin(omega (t_2 - t_1)))) exp((i m omega)/(2 hbar sin(omega T)) ((x_2^2 + x_1^2) cos(omega T) - 2 x_2 x_1)) , T = t_2 - t_1
  $
]

=== 传播子与统计力学

考虑传播子中$vb(x)_1 = vb(x)_2 = vb(x)$的(*封闭*)情形，进一步对坐标$vb(x)$积分，定义
$
  G(t) = integral dd(vb(x), 3) K(vb(x), t; vb(x), 0)
$
这里考虑Hamilton量不含时，令$t_2 = t, t_1 = 0$。根据传播子的定义，得
$
  G(t) = integral dd(vb(x), 3) braket(vb(x), exp(- i/hbar hat(H) t), vb(x))
$
这实际上是*求迹*
$
  G(t) = tr(exp(- i/hbar hat(H) t))
$
而求迹与表象无关，因此可以用任意表象求迹。实际上，利用能量本征态$ket(n)$得到
$
  G(t) & = integral dd(vb(x), 3) braket(vb(x), exp(- i/hbar hat(H) t), vb(x)) \
       & = sum_n integral dd(vb(x), 3) braket(vb(x), e^(- i/hbar hat(H) t), n) braket(n, vb(x)) \
       & = sum_n e^(- i/hbar E_n t) integral dd(vb(x), 3) braket(vb(x), n) braket(n, vb(x)) \
       & = sum_n e^(- i/hbar E_n t) \
       & = sum_n braket(n, e^(- i/hbar hat(H) t), n) \
$
如果做如下替换
$
  (i t)/hbar -> beta = 1/(k_B T)
$
则函数$G(t)$变为
$
  G(t) -> Z(beta) = sum_n e^(- beta E_n)
$
这实际上就是*平衡态量子统计力学的配分函数*，前面的求迹则给出更一般的表达式
$
  Z(beta) = tr(e^(- beta hat(H)))
$
#newpara()
如果对$G(t)$做Fourier变换
$
  tilde(G)(E) = 1/(i hbar) integral_0^(oo) dd(t) G(t) e^(i/hbar E t)
$
利用前面的结果，得到
$
  tilde(G)(E) & = 1/(i hbar) integral_0^(oo) dd(t) sum_n e^(- i/hbar E_n t) e^(i/hbar E t) \
              & = 1/(i hbar) sum_n integral_0^(oo) dd(t) e^(i/hbar (E - E_n) t) \
$
若将$E$加上无穷小的虚部 (加上收敛因子)， $E -> E + i epsilon$，则积分可计算为
$
  tilde(G)(E) = sum_n 1/(E - E_n + i epsilon)
$
如果能将$G(t)$计算出来，则可以通过Fourier变换得到能级结构。

=== 从传播子到路径积分

考虑*传播子*
$
  K(vb(x)_3, t_3; vb(x)_1, t_1) = braket(vb(x)_3, hat(U)(t_3, t_1), vb(x)_1)
$
利用时间演化算符的合成性质
$
  hat(U)(t_3, t_1) = hat(U)(t_3, t_2) hat(U)(t_2, t_1)
$
得到$(t_1 < t_2 < t_3)$
$
  K(vb(x)_3, t_3; vb(x)_1, t_1) &= braket(vb(x)_3, hat(U)(t_3, t_2) hat(U)(t_2, t_1), vb(x)_1) \
  & = integral dd(vb(x)_2, 3) braket(vb(x)_3, hat(U)(t_3, t_2), vb(x)_2) braket(vb(x)_2, hat(U)(t_2, t_1), vb(x)_1) \
  & = integral dd(vb(x)_2, 3) K(vb(x)_3, t_3; vb(x)_2, t_2) K(vb(x)_2, t_2; vb(x)_1, t_1)
$
受此启发，我们考虑从时空点$(vb(x)_1, t_1)$到$(vb(x)_N, t_N)$的传播子
$
  K(vb(x)_N, t_N; vb(x)_1, t_1) = braket(vb(x)_N, hat(U)(t_N, t_1), vb(x)_1)
$
将时间间隔$[t_1, t_N]$等分为$N − 1$个小间隔，每个小间隔为
$
  t_i - t_(i-1) = Delta t = (t_N - t_1)/(N - 1), i = 2, 3, ..., N
$
利用时间演化算符的合成性质
$
  hat(U)(t_N, t_1) = hat(U)(t_N, t_(N-1)) hat(U)(t_(N-1), t_(N-2)) ... hat(U)(t_2, t_1), t_1 < t_2 < ... < t_N
$
在相邻两个时间演化算符之间插入完备性关系
$
  integral dd(vb(x)_i, 3) ket(vb(x)_i) bra(vb(x)_i) = hat(I), i = 2, 3, ..., N-1
$
得到
$
  braket(vb(x)_N, hat(U)(t_N, t_1), vb(x)_1) & = integral dd(vb(x)_(N-1), 3) ... integral dd(vb(x)_2, 3) \
                                             & braket(vb(x)_N, hat(U)(t_N, t_(N-1)), vb(x)_(N-1)) \
                                             & braket(vb(x)_(N-1), hat(U)(t_(N-1), t_(N-2)), vb(x)_(N-2)) ... \
                                             & braket(vb(x)_2, hat(U)(t_2, t_1), vb(x)_1) \
$
即
$
  K(vb(x)_N, t_N; vb(x)_1, t_1) & = integral dd(vb(x)_(N-1), 3) ... integral dd(vb(x)_2, 3) \
  & K(vb(x)_N, t_N; vb(x)_(N-1), t_(N-1)) K(vb(x)_(N-1), t_(N-1); vb(x)_(N-2), t_(N-2)) ... K(vb(x)_2, t_2; vb(x)_1, t_1) \
$
若$hat(H)$不含时，利用Heisenberg绘景中位置算符的本征态
$
  ket(vb(x)\, t) = e^(i/hbar hat(H) t) ket(vb(x))
$
则传播子即为不同时空点之间的*跃迁振幅*
$
  braket(vb(x)_N, hat(U)(t_N, t_1), vb(x)_1) = braket(vb(x)_N\, t_N, vb(x)_1\, t_1)
$
利用完备性关系
$
  integral dd(vb(x), 3) ketbra(vb(x)\, t) = hat(I)
$
得到
$
  braket(vb(x)_N\, t_N, vb(x)_1\, t_1) & = integral dd(vb(x)_(N-1), 3) ... integral dd(vb(x)_2, 3) \
  & braket(vb(x)_N\, t_N, vb(x)_(N-1)\, t_(N-1)) braket(vb(x)_(N-1)\, t_(N-1), vb(x)_(N-2)\, t_(N-2)) ... braket(vb(x)_2\, t_2, vb(x)_1\, t_1) \
$
以上结果表明，*粒子从一个时空点传播到另一个时空点，所有的路径$vb(x)(t)$都是有可能的！*
#figure(
  image("pic/2025-12-21-19-21-32.png", width: 60%),
  numbering: none,
)

=== 路径积分量子化

Fermman通过对双缝干涉实验的思考，结合力学中的最小作用量原理和光学中的Huygens原理，提出了*提出了路径积分量子化*。Fermman认为，*粒子从时空点$(vb(x)_A, t_A)$传播到$(vb(x)_B, t_B)$的几率振幅是所有可能的路径的几率振幅相加*，这样很自然地满足了量子力学的态叠加原理和传播子的合成性质。径积分量子化方案(Fermman公设)的具体内容是：
- 找出连接$(x_A, t_A)$和$(x_B, t_B)$的全部路径$(t_B > t_A)$
- 对每条路径$vb(x)(t)$，计算对应的作用量
  $
    S[vb(x)(t)] = integral_(t_A)^(t_B) dd(t) L(vb(x), dot(vb(x)), t)
  $
  该路径对应的*几率振幅*正比于$exp(i/hbar S[vb(x)(t)])$
- $K(x_B, t_B; x_A, t_A)$即为*所有路径的几率振幅相加*
  $
    K(vb(x)_B, t_B; vb(x)_A, t_A) = c sum_"all path" exp(i/hbar S[vb(x)(t)])
  $
  其中$C$为适当的归一化常数
由于粒子的路径是连续变化的，作用量$S$是$vb(x)(t)$的泛函，上述求和实际上是对$vb(x)(t)$的*泛函积分*
$
  K(vb(x)_B, t_B; vb(x)_A, t_A) = integral_((vb(x)_A, t_A))^((vb(x)_B, t_B)) cal(D)[vb(x)(t)] exp(i/hbar S[vb(x)(t)])
$
但是，常数$C$的取值或者说*积分测度*$D[vb(x)(t)]$的定义尚未明确。

*推导路径积分*：考虑在势场$V(x)$中运动的粒子，Hamilton量为
$
  hat(H) = hat(vb(p))^2/(2m) + V(hat(vb(x)))
$
则传播子为
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
  &= (m/(2 pi hbar i epsilon))^(3/2) exp((i m (vb(x)_n - vb(x)_(n-1))^2)/(2 hbar epsilon)) exp(- (i epsilon)/hbar V(vb(x)_(n-1))) \
  &= (m/(2 pi hbar i epsilon))^(3/2) exp((i epsilon)/hbar (m/2 ((vb(x)_n - vb(x)_(n-1))/epsilon)^2 - V(vb(x)_(n-1))))
$
当$N -> oo, epsilon -> 0$时，高阶小量可忽略，得到
$
  K(vb(x)_B, t_B; vb(x)_A, t_A) & = lim_(N -> oo) (m/(2 pi hbar i epsilon))^(3 N/2) (product_(n=1)^(N-1) (m/(2 pi hbar i epsilon))^(3/2) dd(vb(x)_n))\ & exp((i epsilon)/hbar sum_(n=1)^N (m/2 ((vb(x)_n - vb(x)_(n-1))/epsilon)^2 - V(vb(x)_(n-1)))) \
$
指数的宗量可写成连续形式
$
  (i epsilon)/hbar sum_(n=1)^N (m/2 ((vb(x)_n - vb(x)_(n-1))/epsilon)^2 - V(vb(x)_(n-1))) & = (i)/hbar integral_(t_A)^(t_B) dd(t) (m/2 dot(vb(x))^2 - V(vb(x))) \
  & = (i)/hbar integral_(t_A)^(t_B) dd(t) L(vb(x), dot(vb(x)), t) \
  & = (i)/hbar S[vb(x)(t)]
$
定义*泛函积分测度*
$
  cal(D)[vb(x)(t)] = lim_(N -> oo) (m/(2 pi hbar i epsilon))^(3/2) (product_(n=1)^(N-1) (m/(2 pi hbar i epsilon))^(3/2) dd(vb(x)_n))
$
即可得到Fermman的路径积分表达式
$
  K(vb(x)_B, t_B; vb(x)_A, t_A) = integral_((vb(x)_A, t_A))^((vb(x)_B, t_B)) cal(D)[vb(x)(t)] exp((i)/hbar S[vb(x)(t)])
$
路径积分的积分变量是连续函数，是无穷维积分(甚至不可数)。

注：
- 归一化常数$C$和积分测度中的常数与相互作用无关；
  $
    (m/(2 pi hbar i epsilon))^((3N)/2) = ((m N)/(2 pi hbar i(t_B - t_A)))^(3/2) (N->oo)
  $
  在计算物理可观测量时，这些常数往往会相互抵消，因此不必过分关注它们的具体数值；
- 路径积分中不需要态和算符，只需要经典的Lagrange量和作用量。路径积分量子化与Schrödinger方程(正则量子化) 等价，实际上，上面已经从Schrödinger方程推导出了路径积分，反之亦然；
- 路径积分形式优美简洁，但是却难于计算(无穷维积分)。只有当作用量为$vb(x)(t)$的二次型时(Gauss型)，路径积分才可以被精确计算。

=== 经典极限

当$hbar -> 0$时，路径积分应回到经典极限。当$hbar -> 0$时，$exp(i/hbar S)$是剧烈震荡的函数，大部分相邻路径的贡献之间将导致干涉相消。考虑某个路径$vb(x)(t)$，其相邻路径为$vb(x)(t) + delta vb(x)(t)$(保持两个端点不变)，作用量的变化为
$
  S[vb(x)(t) + delta vb(x)(t)] = S[vb(x)(t)] + delta S
$
对于$delta S >> hbar$的那些路径(即宏观可区分的路径)，相位差$(delta S)/hbar >> 1$，彼此产生*相消干涉*。但是，对于$delta S = 0$的路径，即满足经典力学规律的路径$vb(x)_"c" (t)$的路径，它和附近路径之间的相位差很小，彼此产生*相长干涉*。

在数学上，可以用*驻点相位近似(stationary phase approximation)*或者*鞍点近似(saddle point approximation)*更为严格地说明。考虑积分
$
  I = 1/sqrt(2 pi hbar) integral_(-oo)^oo dd(x) exp((i)/hbar f(x))
$
可以证明，当$hbar -> 0$时，有
$
  I = 1/sqrt(f''(x_"c")) exp((i)/hbar f(x_"c")) (1 + hbar(...) + hbar^2(...) + ...)
$
其中，$x_"c"$是函数$f(x)$的最小值点。

=== 路径积分的计算

- *计算方法 1*：积分变量是连续函数$vb(x)(t)$，将时间*离散化*，变成离散的积分(量子场论：时空离散化 $->$ 时空格点、格点规范场论)
- *计算方法 2*：利用边界条件，将连续函数$vb(x)(t)$展开为*离散级数*，可将积分变量从连续函数变为离散的展开系数(模式)
  $
    vb(x)(t) = sum_n c_n vb(u)_n (t), cal(D)[vb(x)(t)] = product_n dd(c_n)
  $
- *计算方法 3*：这里介绍一种作用量为二次型的路径积分的计算方法：*经典贡献 + 量子涨落*。

考虑一维粒子，最一般的二次型(Gauss型)Lagrange量为
$
  L = a(t) dot(x)^2 + b(t) x dot(x) + c(t) x^2 + d(t) dot(x) + e(t) x + f(t)
$
假设粒子从时空点$(x_A, t_A)$运动到$(x_B, t_B)$的*经典路径*为$x_"c" (t)$，可将任意路径$x(t)$写为
$
  x(t) = x_"c" (t) + y(t)
$
其中$y(t)$为“围绕经典路径的*量子涨落*”，满足边界条件
$
  y(t_A) = 0, y(t_B) = 0
$
Lagrange可以用$x_"c" (t)$和$y(t)$重写为
$
  L[x(t)] = L^((0)) + L^((1)) + L^((2))
$
其中$L^((n))$包含了涨落$y(t)$的$n$次方项
$
  L^((0)) & = a dot(x)_"c"^2 + b x_"c" dot(x)_"c" + c x_"c"^2 + d dot(x)_"c" + e x_"c" + f \
  L^((1)) & = 2a dot(x)_"c" dot(y) + b x_"c" dot(y) + b dot(x)_"c" y + 2 c x_"c" y + d dot(y) + e y \
  L^((2)) & = a dot(y)^2 + b y dot(y) + c y^2 \
$
相应地，作用量可以写为
$
  S[x(t)] = integral_(t_A)^(t_B) dd(t) L[x(t)] = S^((0)) + S^((1)) + S^((2))
$
其中，经典路径对应的作用量为
$
  S^((0)) = S_"cl" = S[x_"c" (t)] = integral_(t_A)^(t_B) dd(t) L^((0))
$
可以证明，涨落$y(t)$的一次方项贡献为零。具体计算为
$
  S^((1)) & = integral_(t_A)^(t_B) dd(t) L^((1)) \
          & = integral_(t_A)^(t_B) dd(t) (dot(y) (2 a dot(x)_"c" + b x_"c" + d) + y (b dot(x)_"c" + 2 c x_"c" + e)) \
          & = - integral_(t_A)^(t_B) dd(t) y (dv(, t) (2 a dot(x)_"c" + b x_"c" + d) - (b dot(x)_"c" + 2 c x_"c" + e)) \
          & = - integral_(t_A)^(t_B) dd(t) y (dv(, t) (pdv(L^((0)), dot(x)_"c")) - pdv(L^((0)), x_"c")) = 0
$
最后，计算路径积分
$
  K(x_B, t_B; x_A, t_A) & = integral_((x_A, t_A))^((x_B, t_B)) cal(D)[x(t)] exp((i)/hbar S[x(t)]) \
$
对积分变量做平移$x(t) -> y(t)$，不改变积分测度，$cal(D)[x(t)] = cal(D)[y(t)]$，得到

$
  K(x_B, t_B; x_A, t_A) & = integral_((0, t_A))^((0, t_B)) cal(D)[y(t)] exp((i)/hbar S[x_"c" (t) + y(t)]) \
                        & = exp((i)/hbar S_"cl") integral_((0, t_A))^((0, t_B)) cal(D)[y(t)] exp((i)/hbar S^((2))[y(t)]) \
$
其中
$
  S^((2)) = integral_(t_A)^(t_B) dd(t) L^((2)) = integral_(t_A)^(t_B) dd(t) (a dot(y)^2 + b y dot(y) + c y^2)
$
所以，*Gauss型Lagrange量对应的几率振幅为*
$
  K(x_B, t_B; x_A, t_A) & = F(t_A,t_B) exp((i)/hbar S_"cl") \
$
其中
$
  F(t_A,t_B) = integral_((0, t_A))^((0, t_B)) cal(D)[y(t)] exp((i)/hbar integral_(t_A)^(t_B) dd(t) (a dot(y)^2 + b y dot(y) + c y^2))
$
是来自量子涨落的贡献。它可以通过方法1和方法2直接计算。下面以自由粒子和谐振子为例介绍一种简便的间接方法(Fermman)。
