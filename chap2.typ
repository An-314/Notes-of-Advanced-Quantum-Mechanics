#import "@preview/scripst:1.1.1": *

= 量子动力学

在经典力学中，动力学就是研究系统的*时间演化*，即在系统内部的*相互作用*以及系统与外界之间的相互作用的支配下，系统的状态如何随时间演化。在正统的量子力学中，时间演化分为两种，一是在相互作用支配下的动力学演化，二是量子测量过程。关于量子测量是否可以纳入动力学演化，仍然是一个 open 的问题。本章讨论量子系统的动力学演化，即“计算量子力学” 。

主要探讨如下内容：
- Schrödinger绘景：时间演化算符和传播子
- 绘景理论：Heisenberg绘景和相互作用绘景
- 绝热演化：量子绝热近似和几何相位
- 密度矩阵：纯态和混合态，量子开放系统

== 时间演化算符

量子系统的状态$ket(ψ(t))$的动力学演化满足Schrödinger方程
$
  i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
$
假设$t_0$时刻系统的初始状态为$ket(psi(t_0))$，则$t$时刻$(t > t_0)$的状态为系统的状态$ket(psi(t))$在形式上可以写为
$
  ket(psi(t)) = hat(U)(t, t_0) ket(psi(t_0))
$
其中，算符$hat(U)(t, t_0)$称为*时间演化算符*。代入Schrödinger方程，可以得到时间演化算符满足的方程
$
  i hbar pdv(, t) hat(U)(t, t_0) = hat(H) hat(U)(t, t_0)
$
此即时间演化算符的演化方程。其Hermite共轭方程为
$
  - i hbar pdv(, t) hat(U)^(dagger)(t, t_0) = hat(U)^(dagger)(t, t_0) hat(H)^dagger
$
因此可以计算得到
$
  i hbar pdv(, t) (hat(U)^(dagger)(t, t_0) hat(U)(t, t_0)) = hat(U)^(dagger)(t, t_0) (hat(H) - hat(H)^dagger) hat(U)(t, t_0)
$
考虑到初始条件$hat(U)(t_0,t_0)=hat(I)$，只要Hamilton量是Hermite算符，就有
$
  hat(U)^(dagger)(t, t_0) hat(U)(t, t_0) = hat(I)
$
因此，时间演化算符是*幺正算符*。

进一步地，我们可以计算
$
  braket(psi(t)) = braket(psi(t_0), hat(U)^(dagger)(t, t_0) hat(U)(t, t_0), psi(t_0)) = braket(psi(t_0))
$
因此，时间演化的幺正性保证了量子态的模方保持不变。对于非相对论粒子系统来说，就是*概率守恒*。

上述演绎也可以反过来，即事先假定演化是幺正的。假设$t_0$时刻系统的初始状态为$ket(ψ(t_0))$，将任意 $t > t_0$ 时刻的状态 $ket(ψ(t))$ 写为
$
  ket(psi(t)) = hat(U)(t, t_0) ket(psi(t_0)) <-> bra(psi(t_0)) = bra(psi(t)) hat(U)^(dagger)(t, t_0)\
  braket(psi(t)) = braket(psi(t_0), hat(U)^(dagger)(t, t_0) hat(U)(t, t_0), psi(t_0)) = braket(psi(t_0))
$
因此，如果时间演化算符是幺正算符，
$
  hat(U)^(dagger)(t, t_0) hat(U)(t, t_0) = hat(I)
$
那么就可以保证量子态的模方$braket(ψ(t))$不随时间变化。

根据时间演化算符的定义，容易证明其满足初始条件
$
  hat(U)(t_0, t_0) = hat(I)
$
以及合成性质
$
  hat(U)(t_2, t_0) = hat(U)(t_2, t_1) hat(U)(t_1, t_0), (t_2 > t_1 > t_0)
$
考虑无穷小时间演化：$t -> t + dd(t)$，即
$
  ket(psi(t + dd(t))) = hat(U)(t + dd(t), t) ket(psi(t))
$
将时间演化算符$hat(U)(t + dd(t), t)$在$dd(t)$处展开，得到
$
  hat(U)(t + dd(t), t) = hat(I) - i hat(Omega)(t) dd(t) + O(dd(t)^2)
$
由于假设$hat(U)$是幺正的，$hat(U)^(dagger) hat(U) = hat(I)$，因此
$
  (1 + i hat(Omega)^(dagger) dd(t))(1 - i hat(Omega) dd(t)) = hat(I) + i (hat(Omega)^(dagger) - hat(Omega)) dd(t) + O(dd(t)^2) = hat(I)
$
所以*生成元*$hat(Omega)$是Hermite算符，即$hat(Omega) = hat(Omega)^(dagger)$。由于$hat(Omega)$具有频率的量纲，引入具有能量量纲的Hamilton算符$hat(H)$，使得
$
  hat(Omega) = hat(H)/hbar
$
这与经典力学类似，Hamilton量亦是时间演化的生成元。无穷小时间演化算符可写为
$
  hat(U)(t + dd(t), t) = hat(I) - i/(hbar) hat(H) dd(t)
$
其中$hbar$实际上就是约化Planck常数，这保证最终得到的时间演化
方程就是Schrödinger方程。

根据时间演化算符的合成性质
$
  hat(U)(t_2, t_0) = hat(U)(t_2, t_1) hat(U)(t_1, t_0)
$
得到
$
  hat(U)(t + dd(t), t_0) = hat(U)(t + dd(t), t) hat(U)(t, t_0)
$
代入上式得到
$
  hat(U)(t + dd(t), t_0) = (hat(I) - i/(hbar) hat(H) dd(t)) hat(U)(t, t_0)\
  i hbar (hat(U)(t + dd(t), t_0) - hat(U)(t, t_0))/dd(t) = hat(H) hat(U)(t, t_0)
$
由于$dd(t)$是无穷小量，上式左边改写为求导，就得到了时间演化算符满足的运动方程
$
  i hbar pdv(, t) hat(U)(t, t_0) = hat(H) hat(U)(t, t_0)
$
此方程两边都作用初态$ket(ψ(t_0))$，就得到了Schrödinger方程
$
  i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
$
在上述“推导”中，利用了两个假设：
1. 时间演化算符是幺正算符；
2. 时间演化算符的生成元就是系统的Hamilton量。

只要求出了时间演化算符，就在形式上求出了量子态的时间演化。下面讨论演化方程
$
  i hbar pdv(, t) hat(U)(t, t_0) = hat(H) hat(U)(t, t_0)
$
的解。这是关于时间的一阶微分方程，初始条件为$hat(U)(t_0, t_0) = hat(I)$。如果Hamilton量$hat(H)(t)$已知，则可求解出时间演化算符$hat(U)(t, t_0)$。

最简单的情形，哈密顿量不显含时间，易证
$
  hat(U)(t, t_0) = exp(- i/(hbar) hat(H) (t - t_0))
$
对于哈密顿量不含时的情形，体系的时间演化归结为求解哈密顿量的本征态，即定态薛定谔方程
$
  hat(H) ket(psi_n) = E_n ket(psi_n)
$
一般情况下，哈密顿量的本征态可能存在简并，此时可引入*力学量完备集*${hat(H), hat(A), hat(B), ...}$，使得
$
  [hat(H), hat(A)] = [hat(H), hat(B)] = [hat(A), hat(B)] = ... = 0\
  [hat(A), hat(B)] = [hat(A), hat(B)] = ... = 0
$
这组力学量完备集存在完备的共同本征态，记为$ket(n)$（$n$代表一组完备的量子数），则时间演化算符可计算为（*能量表象*）
$
  e^(- (i hat(H) (t - t_0))/hbar) = sum_n sum_n' ket(n) braket(n, e^(- (i hat(H) (t - t_0))/hbar), n') bra(n') = sum_n e^(- (i E_n (t - t_0))/hbar) ketbra(n)
$
将初态用这组完备本征态${ket(n)}$展开，$ket(ψ(t_0))) = sum_n c_n ket(n)$，得到
$
  ket(psi(t)) & = hat(U)(t, t_0) ket(psi(t_0)) \
              & = sum_n e^(- (i E_n (t - t_0))/hbar) sum_n c_n ket(n) \
              & = sum_n c_n e^(- (i E_n (t - t_0))/hbar) ket(n)
$
这是我们熟悉的结果。解析计算的难度在于：
+ 哈密顿量的本征值和本征态是否存在解析解；
+ 最后的级数求和是否能存在有限结果。
满足这两个条件的量子系统是很少的（基本都在教科书上）。

如果哈密顿量显含时间，但不同时刻的哈密顿量对易，即
$
  [hat(H)(t_1), hat(H)(t_2)] = 0, (forall t_1 , t_2)
$
则演化方程的解也比较简单，易证明为
$
  hat(U)(t, t_0) = exp(- i/hbar integral_(t_0)^(t) hat(H)(t') dd(t'))
$

对于最一般的情形，即哈密顿量显含时间，且不同时刻的哈密顿量不对易，此时演化方程的形式解为
$
  hat(U)(t, t_0) = 1 + sum_(n=1)^(oo) (- i/hbar)^n integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) dd(t_2) ... integral_(t_0)^(t_(n-1)) dd(t_n) hat(H)(t_1) hat(H)(t_2) ... hat(H)(t_n)
$
此解称为Dyson 级数。式中各个积分中的时间变量满足
$
  t >= t_1 >= t_2 >= ... >= t_n >= t_0
$
这个解可以用迭代方法证明。

将演化方程
$
  i hbar pdv(, t) hat(U)(t, t_0) = hat(H)(t) hat(U)(t, t_0)
$
中的变量$t$改写为$t_1$，两边对$t_1$在区间$[t_0, t]$积分，得到
$
  i hbar eval((hat(U)(t_1, t_0)))_(t_0)^(t) = integral_(t_0)^(t) hat(H)(t_1) hat(U)(t_1, t_0) dd(t_1)
$
利用初始条件$hat(U)(t_0, t_0) = hat(I)$，得到
$
  hat(U)(t, t_0) = hat(I) - i/hbar integral_(t_0)^(t) hat(H)(t_1) hat(U)(t_1, t_0) dd(t_1)
$
这是一个*积分方程*，与微分方程 + 初始条件完全等价。 可以利用迭代法写出形式解。
$
  hat(U)(t_1, t_0) = hat(I) - i/hbar integral_(t_0)^(t_1) hat(H)(t_2) hat(U)(t_2, t_0) dd(t_2)
$
再次代入上式，得到
$
  hat(U)(t, t_0) &= hat(I) - i/hbar integral_(t_0)^(t) hat(H)(t_1) (hat(I) - i/hbar integral_(t_0)^(t_1) hat(H)(t_2) hat(U)(t_2, t_0) dd(t_2)) dd(t_1) \
  &=hat(I) - i/hbar integral_(t_0)^(t) hat(H)(t_1) dd(t_1) + (- i/hbar)^2 integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) hat(H)(t_1) hat(H)(t_2) hat(U)(t_2, t_0) dd(t_2)
$
这样不断迭代下去，我们就得到 Dyson 级数解
$
  hat(U)(t, t_0) = hat(I) + sum_(n=1)^(oo) (- i/hbar)^n integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) dd(t_2) ... integral_(t_0)^(t_(n-1)) dd(t_n) hat(H)(t_1) hat(H)(t_2) ... hat(H)(t_n)
$
需要注意的是，每一项中的哈密顿量算符按照*时序*排列，一般不能交换顺序。但是，此表达式中每个积分上限与前一个积分变量互相依赖，处理起来很困难。最简单的想法是将积分上限都写为$t$，则得到
$
  hat(I) + (- i/hbar) integral_(t_0)^(t) hat(H)(t_1) dd(t_1) + (- i/hbar)^2 (integral_(t_0)^(t) dd(t_1) hat(H)(t_1))^2 + ...
$
这个结果显然是不对的，只有在不同时刻的哈密顿量对易的情形才成立。例如二阶项的结果为
$
  (- i/hbar)^2 integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t) dd(t_2) hat(H)(t_1) hat(H)(t_2)
$
其中要求 t1 > t2 。要将积分上限都扩展为$t$，我们需要引入时序乘积(又称为编时乘积) 的概念。时序乘积$T$使若干个含时算符的乘积从左到右按照时间的大小降序排列。最简单的例子是
$
  T[hat(H)(t_1) hat(H)(t_2)] = theta(t_1 - t_2) hat(H)(t_1) hat(H)(t_2) + theta(t_2 - t_1) hat(H)(t_2) hat(H)(t_1)
$
其中，$theta(x)$为阶越函数。利用时序乘积的定义，我们可以计算
$
  T[integral_(t_0)^t hat(H)(t') dd(t')]^2 &= T integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t) dd(t_2) hat(H)(t_1) hat(H)(t_2) \
  &= integral_(t_0)^t dd(t_1) integral_(t_0)^(t_1) dd(t_2) T[hat(H)(t_1) hat(H)(t_2)]\
  &= integral_(t_0)^t dd(t_1) integral_(t_0)^(t_1) dd(t_2) (theta(t_1 - t_2) hat(H)(t_1) hat(H)(t_2) + theta(t_2 - t_1) hat(H)(t_2) hat(H)(t_1))\
  &= integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) dd(t_2) hat(H)(t_1) hat(H)(t_2) + integral_(t_0)^(t) dd(t_2) integral_(t_0)^(t_2) dd(t_1) hat(H)(t_2) hat(H)(t_1)\
  &= 2 integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) dd(t_2) hat(H)(t_1) hat(H)(t_2)
$
对于更一般的情形，时序乘积定义为
$
  T[hat(H)(t_1) hat(H)(t_2) ... hat(H)(t_n)] = sum theta(tau_1 - tau_2) theta(tau_2 - tau_3) ... theta(tau_(n-1) - tau_n) hat(H)(tau_1) hat(H)(tau_2) ... hat(H)(tau_n)
$
其中$τ_1 , τ_2 , ... , τ_n$是$t_1 , t_2 , ... , t_n$的任意排列，上述求和共包含$n!$项。类似地，我们可以计算
$
  T [integral_(t_0)^(t) hat(H)(t') dd(t')]^n = n! integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) dd(t_2) ... integral_(t_0)^(t_(n-1)) dd(t_n) hat(H)(t_1) hat(H)(t_2) ... hat(H)(t_n)
$
所以，Dyson 级数解中的第$n$阶项可以写为
$
  (- i/hbar)^n integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) dd(t_2) ... integral_(t_0)^(t_(n-1)) dd(t_n) hat(H)(t_1) hat(H)(t_2) ... hat(H)(t_n) = 1/n! (- i/hbar)^n T [integral_(t_0)^(t) hat(H)(t') dd(t')]^n
$
利用这个结果，我们可以将 Dyson 级数解写成简洁的形式
$
  hat(U)(t, t_0) & = T{1 + sum_(n=1)^(oo) 1/n! (- i/hbar)^n [integral_(t_0)^(t) hat(H)(t') dd(t')]^n} \
                 & =T exp(- i/hbar integral_(t_0)^(t) hat(H)(t') dd(t'))
$

== 绘景理论

我们一直都听说，海森堡、波恩和约当在几篇开创性论文中提出了*矩阵力学*。如果认为矩阵力学就是薛定谔方程在离散表象的具体形式，那么我们去阅读那几篇原始论文，必定会一头雾水。因为*表象*(representation) 只是量子力学的一个方面，它决定了量子态和力学量的具体表现形式。而量子力学还存在另一个重要方面，就是*绘景*(picture)，它决定了运动方程的具体表现形式。

假设量子系统的哈密顿量$hat(H)$不显含时间，其本征值和本征态已经求解出来，
$
  hat(H) ket(n) = E_n ket(n)
$
这些能量本征态${ket(n)}$是完备的，构成能量表象。我们来写出薛定谔方程在能量表象下的具体形式
$
  i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
$
的左边计算为
$
  i hbar pdv(, t) ket(psi(t)) & = i hbar pdv(, t) sum_n ketbra(n) ket(psi(t)) \
                              & = sum_n ket(n) (i hbar pdv(, t) braket(n, psi(t))) \
                              & = sum_n ket(n) i hbar pdv(, t) psi_n (t)
$
右边计算为
$
  hat(H) ket(psi(t)) & = sum_n hat(H) ket(n) hat(H) sum_m ketbra(m) ket(psi(t)) \
                     & = sum_n ket(n) sum_m braket(n, hat(H), m) psi_m(t) \
                     & = sum_n ket(n) E_n psi_n(t)
$
利用正交性得到
$
  i hbar pdv(, t) psi_n(t) = E_n psi_n(t)
$
是薛定谔方程在能量表象中的形式，假设初始条件为$ψ_n (t_0) = c_n$，其解很容易求出为
$
  psi_n(t) = c_n e^(-(i E_n (t - t_0)) / hbar)
$
这是我们熟悉的结果。如果仔细检查一下计算过程，可以看到，我们假定基矢$ket(n)$是不随时间变化的。现在，我们考虑“随动”的基矢
$
  ket(tilde(n)(t)) = e^(-(i E_n (t - t_0)) / hbar) ket(n)
$
这些基矢仍然满足
$
  hat(H) ket(tilde(n)(t)) = E_n ket(tilde(n)(t))\
  braket(tilde(n)(t), tilde(m)(t)) = delta_(n,m)\
  sum_n ketbra(tilde(n)(t)) = hat(I)
$
利用这组“随动”的基矢，薛定谔方程的左边计算为
$
  i hbar pdv(, t) ket(psi(t)) & = i hbar pdv(, t) sum_n ketbra(tilde(n)(t)) ket(psi(t)) \
  & = sum_n ((i hbar pdv(, t)) braket(tilde(n)(t), psi(t))) ket(tilde(n)(t)) + ketbra(tilde(n)(t)) (i hbar pdv(, t) ket(psi(t)))\
$
