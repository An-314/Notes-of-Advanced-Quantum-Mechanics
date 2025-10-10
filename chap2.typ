#import "@preview/scripst:1.1.1": *

= 量子动力学

在经典力学中，动力学就是研究系统的*时间演化*，即在系统内部的*相互作用*以及系统与外界之间的相互作用的支配下，系统的状态如何随时间演化。

在正统的量子力学中，时间演化分为两种，一是在相互作用支配下的动力学演化，二是量子测量过程。关于量子测量是否可以纳入动力学演化，仍然是一个 open 的问题。本章讨论量子系统的动力学演化，即“计算量子力学” 。

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
考虑到初始条件$hat(U)(t_0,t_0)=hat(I)$，*只要Hamilton量是Hermite算符*，就有
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

根据时间演化算符的定义，容易证明其满足*初始条件*
$
  hat(U)(t_0, t_0) = hat(I)
$
以及*合成性质*
$
  hat(U)(t_2, t_0) = hat(U)(t_2, t_1) hat(U)(t_1, t_0), (t_2 > t_1 > t_0)
$
#newpara()
考虑*无穷小时间演化*：$t -> t + dd(t)$，即
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
所以*生成元*$hat(Omega)$是Hermite算符，即$hat(Omega) = hat(Omega)^(dagger)$。由于$hat(Omega)$具有频率的量纲，引入具有能量量纲的*Hamilton算符*$hat(H)$，使得
$
  hat(Omega) = hat(H)/hbar
$
这与经典力学类似，*Hamilton量亦是时间演化的生成元*。无穷小时间演化算符可写为
$
  hat(U)(t + dd(t), t) = hat(I) - (i hat(H) dd(t))/(hbar)
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
由于$dd(t)$是无穷小量，上式左边改写为求导，就得到了*时间演化算符满足的运动方程*
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
的解。这是关于时间的*一阶微分方程*，初始条件为$hat(U)(t_0, t_0) = hat(I)$。如果Hamilton量$hat(H)(t)$已知，则可求解出时间演化算符$hat(U)(t, t_0)$。

最简单的情形，*Hamilton量不显含时间*，易证
$
  hat(U)(t, t_0) = exp(- i/(hbar) hat(H) (t - t_0))
$
#newpara()
对于*Hamilton量不含时*的情形，体系的时间演化归结为求解Hamilton量的本征态，即*定态Schrödinger方程*
$
  hat(H) ket(psi_n) = E_n ket(psi_n)
$
一般情况下，Hamilton量的本征态可能存在简并，此时可引入*力学量完备集*${hat(H), hat(A), hat(B), ...}$，使得
$
  [hat(H), hat(A)] = [hat(H), hat(B)] = [hat(H), hat(C)] = ... = 0\
  [hat(A), hat(B)] = [hat(B), hat(C)] = [hat(C), hat(A)] = ... = 0\
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
+ Hamilton量的本征值和本征态是否存在解析解；
+ 最后的级数求和是否能存在有限结果。
满足这两个条件的量子系统是很少的（基本都在教科书上）。

如果*Hamilton量显含时间，但不同时刻的Hamilton量对易*，即
$
  [hat(H)(t_1), hat(H)(t_2)] = 0, (forall t_1 , t_2)
$
则演化方程的解也比较简单，易证明为
$
  hat(U)(t, t_0) = exp(- i/hbar integral_(t_0)^(t) hat(H)(t') dd(t'))
$
#example(subname: [粒子自旋磁矩])[
  粒子自旋磁矩与沿$z$方向的变化磁场相互作用，Hamilton量为
  $
    hat(H)(t) = - mu_0 B(t) sigma_z
  $
  时间演化算符为
  $
    hat(U)(t, t_0) = exp((i mu_0 sigma_z)/hbar integral_(t_0)^(t) B(t') dd(t'))
  $
]
#newpara()
对于最一般的情形，即*Hamilton量显含时间，且不同时刻的Hamilton量不对易*，此时演化方程的形式解为
$
  hat(U)(t, t_0) = 1 + sum_(n=1)^(oo) (- i/hbar)^n integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t_1) dd(t_2) ... integral_(t_0)^(t_(n-1)) dd(t_n) hat(H)(t_1) hat(H)(t_2) ... hat(H)(t_n)
$
此解称为*Dyson 级数*。式中各个积分中的时间变量满足
$
  t >= t_1 >= t_2 >= ... >= t_n >= t_0
$
这个解可以用迭代方法证明。

#proof[
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
]
需要注意的是，每一项中的Hamilton量算符按照*时序*排列，一般不能交换顺序。但是，此表达式中每个积分上限与前一个积分变量互相依赖，处理起来很困难。最简单的想法是将积分上限都写为$t$，则得到
$
  hat(I) + (- i/hbar) integral_(t_0)^(t) hat(H)(t_1) dd(t_1) + (- i/hbar)^2 (integral_(t_0)^(t) dd(t_1) hat(H)(t_1))^2 + ...
$
这个结果显然是不对的，只有在不同时刻的Hamilton量对易的情形才成立。例如二阶项的结果为
$
  (- i/hbar)^2 integral_(t_0)^(t) dd(t_1) integral_(t_0)^(t) dd(t_2) hat(H)(t_1) hat(H)(t_2)
$
其中要求$t_1 > t_2$。要将积分上限都扩展为$t$，我们需要引入*时序乘积*(又称为*编时乘积*) 的概念。时序乘积$T$使若干个含时算符的乘积从左到右按照时间的大小降序排列。最简单的例子是
$
  T[hat(H)(t_1) hat(H)(t_2)] = theta(t_1 - t_2) hat(H)(t_1) hat(H)(t_2) + theta(t_2 - t_1) hat(H)(t_2) hat(H)(t_1)
$
其中，$theta(x)$为*阶越函数*。利用时序乘积的定义，我们可以计算
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

我们一直都听说，Heisenberg、Born和Jordan在几篇开创性论文中提出了*矩阵力学*。如果认为矩阵力学就是Schrödinger方程在离散表象的具体形式，那么我们去阅读那几篇原始论文，必定会一头雾水。因为*表象*(representation) 只是量子力学的一个方面，它决定了量子态和力学量的具体表现形式。而量子力学还存在另一个重要方面，就是*绘景*(picture)，它决定了运动方程的具体表现形式。

假设量子系统的Hamilton量$hat(H)$不显含时间，其本征值和本征态已经求解出来，
$
  hat(H) ket(n) = E_n ket(n)
$
这些能量本征态${ket(n)}$是完备的，构成*能量表象*。我们来写出Schrödinger方程在能量表象下的具体形式
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
                     & = sum_n ket(n) sum_m braket(n, hat(H), m) psi_m (t) \
                     & = sum_n ket(n) E_n psi_n (t)
$
利用正交性得到
$
  i hbar pdv(, t) psi_n (t) = E_n psi_n (t)
$
*是Schrödinger方程在能量表象中的形式*，假设初始条件为$ψ_n (t_0) = c_n$，其解很容易求出为
$
  psi_n (t) = c_n e^(-(i E_n (t - t_0)) / hbar)
$
这是我们熟悉的结果。如果仔细检查一下计算过程，可以看到，*我们假定基矢$ket(n)$是不随时间变化的*。

现在，我们考虑“随动”的基矢
$
  ket(tilde(n)(t)) = e^(-(i E_n (t - t_0)) / hbar) ket(n)
$
这些基矢仍然满足
$
  hat(H) ket(tilde(n)(t)) = E_n ket(tilde(n)(t))\
  braket(tilde(n)(t), tilde(m)(t)) = delta_(n,m)\
  sum_n ketbra(tilde(n)(t)) = hat(I)
$
利用这组“随动”的基矢，Schrödinger方程的左边计算为
$
  i hbar pdv(, t) ket(psi(t)) & = i hbar pdv(, t) sum_n ketbra(tilde(n)(t)) ket(psi(t)) \
  & = sum_n ((i hbar pdv(, t) ket(tilde(n)(t))) braket(tilde(n)(t), psi(t)) + ket(psi(t)) (i hbar pdv(, t) ketbra(tilde(n)(t))))\
  &= sum_n ket(tilde(n)(t)) (E_n tilde(psi)_n (t) + i hbar pdv(, t) tilde(psi)_n (t))
$
其中$tilde(psi)_n (t)$定义为
$
  tilde(psi)_n (t) = braket(tilde(n)(t), psi(t)) = e^((i E_n (t - t_0)) / hbar) braket(n, psi(t)) = e^((i E_n (t - t_0)) / hbar) psi_n (t)
$
Schrödinger方程的右边计算为
$
  hat(H) ket(psi(t)) & = sum_n ketbra(tilde(n)(t)) hat(H) sum_m ketbra(tilde(m)(t)) ket(psi(t)) \
                     & = sum_n ket(tilde(n)(t)) sum_m braket(tilde(n)(t), hat(H), tilde(m)(t)) tilde(psi)_m (t) \
                     & = sum_n ket(tilde(n)(t)) E_n tilde(psi)_n (t)
$
将左边和右边计算的结果带入Schrödinger方程，得到
$
  i hbar pdv(, t) tilde(psi)_n (t) = 0
$
即：*波函数不随时间演化*！

波函数不随时间演化，那么什么在随时间变化呢？我们来考虑力学量$hat(A)$的矩阵元，采用“随动”基矢表示为
$
  tilde(A)_(m n)(t) = braket(tilde(m)(t), hat(A), tilde(n)(t))
$
假设$hat(A)$不显含时间，计算其对时间的导数
$
  dv(, t) tilde(A)_(m n)& = dv(, t) braket(tilde(m)(t), hat(A), tilde(n)(t)) \
  & = (dv(, t) bra(tilde(m)(t))) hat(A) ket(tilde(n)(t)) + bra(tilde(m)(t)) hat(A) (dv(, t) ket(tilde(n)(t)))\
  & = (i E_m) / hbar braket(tilde(m)(t), hat(A), tilde(n)(t)) - (i E_n) / hbar braket(tilde(m)(t), hat(A), tilde(n)(t))\
  & = 1/(i hbar) braket(tilde(m)(t), [hat(A), hat(H)], tilde(n)(t))
$
写成矩阵形式，就是Born和Jordan在1925年提出的*Heisenberg方程*
$
  dv(, t) tilde(A) = 1/(i hbar) [tilde(A), hat(H)]
$
这就是矩阵力学的基本方程。

仔细考察上面的计算过程，我们不难发现，选取“随动”基矢相当于对量子态和力学量做了如下变换：
$
  ket(psi(t)) & -> e^((i hat(H) (t - t_0)) / hbar) ket(psi(t)) = ket(tilde(psi)(t)) \
       hat(A) & -> e^((i hat(H) (t - t_0)) / hbar) hat(A) e^(-(i hat(H) (t - t_0)) / hbar) = tilde(A)(t)
$
由于$ket(tilde(psi)(t))$满足时间演化
$
  ket(psi(t)) = e^(- (i hat(H) (t - t_0)) / hbar) ket(psi(t_0))
$
因此，在此变换下，$ket(psi(t)) -> ket(psi(t_0))$，即：*量子态不随时间演化*。也可以从用力学量的平均值来理解。在$t$时刻，力学量$hat(A)$的平均值为
$
  braket(psi(t), hat(A), psi(t)) = braket(psi(t_0), e^((i hat(H) (t - t_0)) / hbar) hat(A) e^(-(i hat(H) (t - t_0)) / hbar), psi(t_0)) = braket(tilde(psi)(t_0), tilde(A)(t), tilde(psi)(t_0))
$
平均值随着时间演化，按照左边，我们认为是由于量子态在演化，也可以按照右边，等效地认为是由于*力学量在演化*。

在上述例子的基础上，我们考虑更一般的*含时幺正变换*$hat(Q)(t)$
$
  hat(Q)^dagger (t) hat(Q)(t) = hat(Q)(t) hat(Q)^dagger (t) = hat(I)
$
它将量子态$ket(psi(t))$和力学量$hat(A)$变换为
$
  ket(psi(t)) -> ket(psi^Q (t)) = hat(Q) ket(psi(t))\
  hat(A) -> hat(A)^Q = hat(Q) hat(A) hat(Q)^dagger = hat(Q)hat(A) hat(Q)^(-1)
$
可以证明，在含时幺正变换下，量子系统的所有观测结果都是不变的，包括力学量的本征值、取值概率和平均值。
- *力学量的平均值不变*
  $
    braket(psi(t), hat(A), psi(t)) = braket(psi(t), hat(Q)^dagger (t) hat(Q)(t) hat(A) hat(Q)^dagger (t) hat(Q)(t), psi(t)) = braket(psi^Q (t), hat(A)^Q, psi^Q (t))
  $
- *力学量之间的对易关系不变*

  假设$[hat(A), hat(B)] = hat(C)$，则
  $
    [hat(A)^Q, hat(B)^Q] & = hat(Q)(t) [hat(A), hat(B)] hat(Q)^dagger (t) \
                         & = hat(Q)(t) hat(C) hat(Q)^dagger (t) = hat(C)^Q
  $
  注意：两个时间必须相同才成立，即等时对易关系。
- *力学量的本征值和取值概率不变*

  假设$hat(A) ket(a) = a ket(a)$，则
  $
       & hat(Q)^dagger hat(A)^Q hat(Q) ket(a) = a ket(a) \
    => & hat(A)^Q (hat(Q) ket(a)) = a (hat(Q) ket(a)) \
    => & hat(A)^Q ket(a^Q) = a ket(a^Q)
  $
- *力学量的取值几率不变*

  系统处于量子态$ket(psi(t))$时，对力学量$hat(A)$进行测量。将$ket(psi(t))$在$hat(A)$的本征态${ket(a)}$上展开
  $
    ket(psi(t)) = sum_a c_a (t) ket(a)
  $
  则得到结果$a$的几率为$abs(c_a (t))^2$。进行$Q$变换以后，将$ket(psi^Q (t))$用$hat(A)^Q$的本征值展开
  $
    ket(psi^Q (t)) = sum_a c_a^Q (t) ket(a^Q)
  $
  两边同时作用$hat(Q)^dagger (t)$，得到
  $
    c^Q_a (t) = braket(a^Q, psi^Q (t)) = braket(a, hat(Q)^dagger (t) hat(Q)(t), psi(t)) = braket(a, psi(t)) = c_a (t)
  $

对系统进行$Q$变换以后，*量子态$ket(psi^Q (t))$的演化方程*也将发生变化。直接计算得到
$
  i hbar pdv(, t) ket(psi^Q (t)) & = i hbar pdv(, t) (hat(Q) ket(psi(t))) \
                                 & = hat(Q) i hbar pdv(, t) ket(psi(t)) + i hbar pdv(, t) hat(Q) ket(psi(t)) \
                                 & = hat(Q) hat(H) ket(psi(t)) + i hbar pdv(hat(Q), t) ket(psi(t)) \
                                 & = (hat(Q) hat(H) hat(Q)^dagger + i hbar pdv(hat(Q), t) hat(Q)^dagger) ket(psi^Q (t))
$
所以，量子态$ket(psi^Q (t))$满足的演化方程为
$
  i hbar pdv(, t) ket(psi^Q (t)) = (hat(H)^Q (t) + i hbar pdv(hat(Q), t) hat(Q)^dagger) ket(psi^Q (t))
$
进行$Q$变换以后，*力学量$A^Q$也将随时间演化*。直接计算得到
$
  dv(, t) hat(A)^Q & = dv(, t) (hat(Q) hat(A) hat(Q)^dagger) \
  & = (dv(, t) hat(Q)) hat(A) hat(Q)^dagger + hat(Q) (pdv(, t) hat(A)) hat(Q)^dagger + hat(Q) hat(A) (dv(, t) hat(Q)^dagger) \
  & = dv(hat(Q), t) hat(Q)^dagger hat(Q) hat(A) hat(Q)^dagger + hat(Q) hat(A) hat(Q)^dagger hat(Q) dv(hat(Q)^dagger, t) + hat(Q) pdv(hat(A), t) hat(Q)^dagger
$
以及
$
  dv(hat(Q), t) hat(Q)^dagger + hat(Q) pdv(hat(Q)^dagger, t) = dv(hat(Q) hat(Q)^dagger, t) = 0
$
所以，力学量$hat(A)^Q$满足的演化方程为
$
  dv(, t) hat(A)^Q = [dv(hat(Q), t) hat(Q)^dagger, hat(A)^Q] + hat(Q) pdv(hat(A), t) hat(Q)^dagger
$
根据以上一般性的讨论，我们发现，从Schrödinger方程
$
  i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
$
出发，通过含时幺正变换$hat(Q)(t)$，可以得到量子体系时间演化的其他等价形式。这个含时幺正变换 $hat(Q)(t)$ 被称为绘景变换，变换后的时间演化形式称为 $Q$ 绘景。

我们将没有经过变换的时间演化形式，即原始的Schrödinger方程演化形式，称为Schrödinger绘景。在理论上，如果找到一个含时幺正变换$hat(Q)$，就可以创建一个新的绘景。除了Schrödinger绘景，常用的其他绘景实际上只有两种：Heisenberg绘景和相互作用绘景（又称为Dirac绘景）。

#theorem(subname: [绘景变换])[
  量子力学中，所有含时幺正变换$hat(Q)(t)$都定义了一种新的绘景。不同绘景下的量子态和力学量的取值概率、平均值和对易关系都是相同的。

  力学量$hat(A)$和量子态$ket(psi(t))$在$Q$绘景下的变换为
  $
    ket(psi(t)) -> ket(psi^Q (t)) = hat(Q) ket(psi(t))\
    hat(A) -> hat(A)^Q = hat(Q) hat(A) hat(Q)^dagger = hat(Q)hat(A) hat(Q)^(-1)
  $
  量子态$ket(psi^Q (t))$和力学量$hat(A)^Q$满足的演化方程为
  $
    i hbar pdv(, t) ket(psi^Q (t)) = (hat(H)^Q (t) + i hbar pdv(hat(Q), t) hat(Q)^dagger) ket(psi^Q (t))\
    dv(, t) hat(A)^Q = [dv(hat(Q), t) hat(Q)^dagger, hat(A)^Q] + hat(Q) pdv(hat(A), t) hat(Q)^dagger
  $
]

=== Heisenberg绘景

在Schrödinger绘景中，力学量不随时间演化，量子态的演化为
$
  ket(psi(t)) = hat(U)(t, t_0) ket(psi(t_0)) <=> i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
$
其中$hat(U)(t, t_0)$ 就是前面研究过的时间演化算符。取含时幺正变换为
$
  hat(Q)(t) = hat(U)^dagger (t, t_0) = hat(U)^(-1) (t, t_0)
$
得到的绘景称为*Heisenberg绘景*。量子态和力学量变换为
$
  ket(psi^"H" (t)) = hat(U)^dagger (t, t_0) ket(psi(t)) = ket(psi(t_0))\
  hat(A)^"H" (t) = hat(U)^dagger (t, t_0) hat(A) hat(U)(t, t_0)
$
通常用上标$"H"$表示Heisenberg绘景 (必要时，用上标$"S"$表示Schrödinger绘景)。

*在Heisenberg绘景中，量子态不随时间演化，力学量随时间演化*。因此我们需要推导出力学量的演化方程。

前面已经得出，对于任意含时幺正变换$hat(Q)$，力学量$hat(A)^Q$的演化方程为
$
  dv(, t) hat(A)^Q = [dv(hat(Q), t) hat(Q)^dagger, hat(A)^Q] + hat(Q) pdv(hat(A), t) hat(Q)^dagger
$
取$hat(Q)(t) = hat(U)^dagger (t, t_0)$，并利用演化方程
$
  dv(hat(Q), t) hat(Q)^dagger = - hat(Q) dv(hat(Q)^dagger, t) = - hat(U)^dagger dv(hat(U), t) = - 1/(i hbar) hat(U)^dagger hat(H) hat(U) = - 1/(i hbar) hat(H)^"H" (t)
$
带入前式得到力学量的演化方程，即*Heisenberg方程*
$
  dv(, t) hat(A)^"H" (t) = 1/(i hbar) [hat(A)^"H" (t), hat(H)^"H" (t)] + hat(U)^dagger pdv(hat(A), t) hat(U)
$
需要注意的是，在定义Heisenberg绘景时，我们并未限定Schrödinger绘景中的Hamilton量$hat(H)$不显含时间。如果$hat(H)$不显含时间，则
$
  hat(U)(t, t_0) = exp(- (i hat(H) (t - t_0)) / hbar) => hat(Q)(t) = exp((i hat(H) (t - t_0)) / hbar)
$
在Heisenberg绘景中，Hamilton量为
$
  hat(H)^"H" (t) = hat(U)^dagger (t, t_0) hat(H) hat(U)(t, t_0) = hat(H)
$
此时两种绘景中的Hamilton量一致
$
  hat(H)^"H" = hat(H)^"S" = hat(H)
$

==== 算力学量的时间演化

考虑不显含时间的力学量$hat(A)$，在Heisenberg绘景中求解其时间演化，
最一般的方法是求解Heisenberg方程
$
  dv(, t) hat(A)^"H" (t) = 1/(i hbar) [hat(A)^"H" (t), hat(H)^"H" (t)]
$
这是一阶微分方程，初始条件为$hat(A)^"H" (t_0) = hat(A)$。计算Heisenberg方程右边的对易关系，可以利用技巧
$
  [hat(A)^"H" (t), hat(H)^"H" (t)] & = hat(U)^dagger [hat(A), hat(H)] hat(U)
$
对于$hat(H)^"S" (t)$含时的系统，即使$hat(H)^"H" (t)$很难求出，Heisenberg方程仍然有可能写出并求解。

对于Hamilton量$hat(H)^"S"$不含时的情形，也可以利用Baker-Hausdorff公式
$
  e^(hat(X)) hat(Y) e^(-hat(X)) = hat(Y) + [hat(X), hat(Y)] + 1/2! [hat(X), [hat(X), hat(Y)]] + 1/3! [hat(X), [hat(X), [hat(X), hat(Y)]]] + ...
$
#proof[
  引入函数$f(lambda) = e^(lambda hat(X)) hat(Y) e^(- lambda hat(X))$，则
  $
    f(0) = hat(Y)\
    f(1) = e^(hat(X)) hat(Y) e^(-hat(X))
  $
  对$lambda$求导得
  $
    dv(f, t) &= e^(lambda hat(X)) hat(X) hat(Y) e^(- lambda hat(X)) - e^(lambda hat(X)) hat(Y) hat(X) e^(- lambda hat(X)) \
    &= e^(lambda hat(X)) [hat(X), hat(Y)] e^(- lambda hat(X))\
    dv(f, lambda, 2) &= e^(lambda hat(X)) [hat(X), [hat(X), hat(Y)]] e^(- lambda hat(X))\
  $
  这样，我们可以写出Taylor展开
  $
    f(lambda) = f(0) + sum_(n=1)^(oo) lambda^n/n! eval(dv(f, lambda, n))_(lambda=0)
  $
  最后令$lambda=1$得到到Baker-Hausdorff公式。
]
应用于力学量$hat(A)$，根据$hat(A)^"H" (t) = hat(U)^dagger (t, t_0) hat(A) hat(U)(t, t_0)$，得到
$
  hat(A)^"H" (t) = hat(A) + (i t)/hbar [hat(H), hat(A)] + 1/2! ((i t)/hbar)^2 [hat(H), [hat(H), hat(A)]] + 1/3! ((i t)/hbar)^3 [hat(H), [hat(H), [hat(H), hat(A)]]] + ...
$
不失一般性，从这里开始，令$t_0 = 0$。一般来说，对于如下的情况上述方法可以奏效：
+ 上述级数中的对易子在某一阶等于零，从而只需要计算有限阶
+ 上述级数中的对易子出现周期性结果，从而可以写成级数形式并求和

==== 谐振子

谐振子是量子理论中非常基本且重要的一个模型。考虑一维谐振子，其Hamilton量为
$
  hat(H) = hat(p)^2/(2 m) + 1/2 m omega^2 hat(x)^2
$
其中$hat(q)$和$hat(p)$应该更一般地理解为系统的广义坐标和广义动量，满足对易关系
$
  [hat(q), hat(p)] = i hbar
$
在本科量子力学中已经用波动力学解法和代数解法求解了其能量本征值问题
$
  E_n = (n + 1/2) hbar omega, (n = 0, 1, 2, ... )\
$
由于Hamilton量不含时，Heisenberg绘景中的时间演化可以直接利用Baker-Hausdorff公式计算
$
  hat(q)^"H" (t) & = hat(q) + (i t)/hbar [hat(H), hat(q)] + 1/2! ((i t)/hbar)^2 [hat(H), [hat(H), hat(q)]] + 1/3! ((i t)/hbar)^3 [hat(H), [hat(H), [hat(H), hat(q)]]] + ...\
$
利用
$
  [hat(H), hat(q)] = - (i hbar)/m hat(p)\
  [hat(H), hat(p)] = i hbar m omega^2 hat(q)
$
得到
$
  hat(q)^"H" (t) &= hat(q) (1 - 1/2! omega^2 t^2 + 1/4! omega^4 t^4 - ...) + hat(p)/(m omega) (omega t - 1/3! omega^3 t^3 + 1/5! omega^5 t^5 - ...)\
$
$
  hat(q)^"H" (t) & = hat(q) cos(omega t) + hat(p)/(m omega) sin(omega t)
$
同样可以得到
$
  hat(p)^"H" (t) = hat(p) cos(omega t) - m omega hat(q) sin(omega t)
$
另一种方法是利用Heisenberg方程求解。得到
$
  dv(, t) hat(q)^"H" (t) &= 1/(i hbar) [hat(q)^"H" (t), hat(H)] = i/hbar hat(U)^dagger [hat(q), hat(H)] hat(U) = (hat(p)^"H" (t))/m\
  dv(, t) hat(p)^"H" (t) &= 1/(i hbar) [hat(p)^"H" (t), hat(H)] = i/hbar hat(U)^dagger [hat(p), hat(H)] hat(U) = - m omega^2 hat(q)^"H" (t)
$
这是关于$hat(q)^"H" (t)$和$hat(p)^"H" (t)$的耦合一阶微分方程组。为了解耦合，将其改写成如下的形式
$
  dv(hat(q)^"H" (t), t) = - i omega ((i hat(p)^"H" (t))/(m omega))\
  dv(, t)(i hat(p)^"H" (t))/(m omega) = i omega hat(q)^"H" (t)
$
将两个方程分别相加和相减，得到
$
  dv(, t)(hat(q)^"H" (t) + (i hat(p)^"H" (t))/(m omega)) = - i omega (hat(q)^"H" (t) + (i hat(p)^"H" (t))/(m omega))\
  dv(, t)(hat(q)^"H" (t) - (i hat(p)^"H" (t))/(m omega)) = i omega (hat(q)^"H" (t) - (i hat(p)^"H" (t))/(m omega))
$
因此，可引入新算符
$
  hat(a)_"H" (t) = C (hat(q)^"H" (t) + (i hat(p)^"H" (t))/(m omega))\
  hat(a)^dagger_"H" (t) = C (hat(q)^"H" (t) - (i hat(p)^"H" (t))/(m omega))\
$
其中$C$为任意实常数。算符$hat(a)_"H" (t)$和$hat(a)^dagger_"H" (t)$满足方程
$
  dv(, t) hat(a)_"H" (t) = - i omega hat(a)_"H" (t)\
  dv(, t) hat(a)^dagger_"H" (t) = i omega hat(a)^dagger_"H" (t)
$
其解很容易求得为
$
  hat(a)_"H" (t) = hat(a) e^(- i omega t)\
  hat(a)^dagger_"H" (t) = hat(a)^dagger e^(i omega t)
$
利用
$
  hat(q) = 1/(2 C) (hat(a)_"H" (t) + hat(a)^dagger_"H" (t))\
  hat(p) = (m omega)/(2 i C) (hat(a)_"H" (t) - hat(a)^dagger_"H" (t))
$
即可求得$hat(q)^"H" (t)$和$hat(p)^"H" (t)$的解，与之前的结果一致。

== 谐振子的相干态

在前面的计算中，令
$
  C = sqrt((m omega)/(2 hbar))
$
得到
$
  hat(a) = sqrt((m omega)/(2 hbar)) (hat(q) + (i hat(p))/(m omega))\
  hat(a)^dagger = sqrt((m omega)/(2 hbar)) (hat(q) - (i hat(p))/(m omega))
$
算符$hat(a)$和$hat(a)^dagger$*无量纲*，称为*下降算符和上升算符*，或者说*湮灭算符和产生算符*，满足对易关系
$
  [hat(a), hat(a)^dagger] = 1
$
谐振子的Hamilton量可表达为
$
  hat(H) = hbar omega (hat(a)^dagger hat(a) + 1/2)
$
#newpara()

引入粒子数算符
$
  hat(N) = hat(a)^dagger hat(a)
$
可以证明 (本科量子力学)，$hat(N)$的本征值为非负整数，即
$
  hat(N) ket(n) = n ket(n), (n = 0, 1, 2, ... )
$
算符$hat(a)$和$hat(a)^dagger$的作用为
$
  hat(a) ket(n) = sqrt(n) ket(n - 1)\
  hat(a)^dagger ket(n) = sqrt(n + 1) ket(n + 1)
$
尤其是（“湮灭真空”）
$
  hat(a) ket(0) = 0
$
显然，$hat(N)$的本征态就是能量本征态
$
  hat(H) ket(n) = E_n ket(n)\
  E_n = (n + 1/2) hbar omega, (n = 0, 1, 2, ... )\
$
所有的激发态可以通过将产生算符$hat(a)^dagger$作用在基态$ket(0)$（“真空态”）不断作用得到。利用
$
  ket(n+1) = 1/sqrt(n + 1) hat(a)^dagger ket(n)
$
得到
$
  ket(n) = 1/sqrt(n!) (hat(a)^dagger)^n ket(0)
$
进一步地，进入坐标表象，可以求得各个能级的波函数$psi_n (q)= braket(q, n)$。对于基态，利用$hat(a) ket(0) = 0$的得到
$
  bra(q) (hat(q) + (i hat(p))/(m omega)) ket(0) = 0\
  integral dd(q') bra(q) (hat(q) + (hbar)/(m omega) hat(p)) ketbra(q') ket(0) = 0\
$
利用
$
  braket(q, hat(q), q') = q delta(q - q')\
  braket(q, hat(p), q') = - i hbar pdv(, q') delta(q - q')
$
得到
$
  (q + (hbar)/(m omega) dv(, q)) psi_0 (q) = 0
$
其归一化的解即为基态波函数
$
  psi_0 (q) = ((m omega)/(pi hbar))^(1/4) exp(- (m omega q^2)/(2 hbar))
$
对于任意激发态$ket(n)$，得到
$
  psi_n (q) & = 1/sqrt(n!) braket(q, (hat(a)^dagger)^n ket(0)) \
            & = 1/sqrt(n!) (sqrt((m omega)/(2 hbar)))^n braket(q, (hat(q) - (i hat(p))/(m omega))^n, 0) \
$
插入$n$个$q$表象的完备性关系
$
  psi_n (x) = 1/sqrt(n!) (sqrt((m omega)/(2 hbar)))^n integral dd(q_1) dd(q_2) ... dd(q_n) \ braket(q, (hat(q) - (i hat(p))/(m omega)), q_1) braket(q_1, (hat(q) - (i hat(p))/(m omega)), q_2) ... braket(q_(n-1), (hat(q) - (i hat(p))/(m omega)), q_n) braket(q_n, 0) \
$
将基本矩阵元代入，积分后得到
$
  psi_n (q) = 1/sqrt(2^n n!) ((m omega)/(pi hbar))^n (q - (hbar)/(m omega) dv(, q))^n psi_0 (q)
$

#newpara()

*相干态*的引出：在前面海森堡绘景的计算中，谐振子的坐标和动量算符的解为
$
  hat(q)^"H" (t) = hat(q) cos(omega t) + hat(p)/(m omega) sin(omega t)\
  hat(p)^"H" (t) = hat(p) cos(omega t) - m omega hat(q) sin(omega t)
$
虽然在形式上与经典谐振子的解一致，但是若初态为能量本征态(定态)，则坐标和动量的平均值为
$
  braket(n, hat(q)^"H" (t), n) = 0\
  braket(n, hat(p)^"H" (t), n) = 0
$
所以，在力学量平均值的意义上，量子谐振子与经典谐振子完全不一样 (基态$n = 0$是例外)。

问：能不能找到一种初态$ket(z)$，力学量平均值的演化与经典谐振子一致？(考虑能量时须将零点能去除)

这样一种初态真的存在，它就是湮灭算符$hat(a)$的本征态
$
  hat(a) ket(z) = z ket(z) <=> bra(z) hat(a)^dagger = z^* bra(z)
$
由于$hat(a)$不是厄米算符，所以$z$一般是复数。显然，$ket(z)$可以用能量本征态展开
$
  ket(z) = sum_(n=0)^(oo) ketbra(n, z) ket(n)
$
下面将会看到，除了特例$z=0$外，$ket(z)$是所有的能量本征态的相干叠加。利用能量本征态的表达式
$
  ket(n) = 1/sqrt(n!) (hat(a)^dagger)^n ket(0)\
  bra(n) = 1/sqrt(n!) bra(0) (hat(a))^n
$
得到展开系数
$
  braket(n, z) = 1/sqrt(n!) bra(0) (hat(a))^n ket(z) = z^n/sqrt(n!) braket(0, z)
$
考虑态$ket(z)$的归一化条件 $braket(z)=1$，得到
$
  1 = braket(z) = sum_(n=0)^(oo) abs(braket(n, z))^2 = abs(braket(0, z))^2 sum_(n=0)^(oo) abs(z)^(2 n)/n! = abs(braket(0, z))^2 e^(abs(z)^2)
$
不失一般性，取
$
  braket(0, z) = e^(- abs(z)^2/2)
$
就有
$
  braket(n, z) = z^n/sqrt(n!) e^(- abs(z)^2/2)
$
所以$ket(z)$可以用能量本征态表达为
$
  ket(z) = e^(- abs(z)^2/2) sum_(n=0)^(oo) z^n/sqrt(n!) ket(n)
$
这样的态被称为*相干态*，本征值$z$可以取所有的复数。注意：只有$z = 0$的相干态才与基态$ket(0)$重合，$z = 1, 2, 3, dots$的那些态与能量本征态$ket(n)$是不同的。

