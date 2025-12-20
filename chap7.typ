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
考虑不随时间变化的磁场。可以证明，在矢量势的规范变换
$
  vb(A) -> tilde(vb(A)) = vb(A) + grad Lambda(vb(x))
$
下，Schrödinger方程亦是规范不变的。定义幺正算符
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

正则量子化方案与经典力学的Hamilton形式对应，因此将Hamilton量放到了一个核心的位置。所以一个很自然的问题就是：可以不可以从经典系统的Lagrange量出发来进行量子化呢？确实是可以的。

本小节将学习另一种量子化方案：*路径积分量子化*。

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

若粒子在$t_1$时刻处于空间$vb(x)_1$ 处 (位置本征态)，则传播子表示在以后任意 $t_2$ 时刻 ($t_2 > t_1$) 粒子处于空间点 $vb(x)_2$ 处的概率幅。简而言之，粒子从时空点 ($t_1$ , $vb(x)_1$ ) 传播到 ($t_2$ , $vb(x)_2$ ) 的概率振幅。若粒子在 $t_1$ 时刻并不处于位置本征态，则利用传播子和 $t_1$ 时刻的波函数可以得到以后任意 $t_2$ 时刻的波函数。这实际上就是符合因果律的幺正时间演化在坐标表象的体现。
