#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

= 量子散射理论

散射或者碰撞是自然界中常见的物理现象，也是物理学家探索微观世界的物质结构、运动规律和基本相互作用的最重要的实验手段。本章我们探讨的是非相对论量子力学框架内的势散射理论，但是其中的一些基本概念是普遍的。量子势散射理论有两种物理图像：
- *定态散射图像*：求解粒子在势场$V(vb(r))$中的定态问题。粒子的入射能量$E$就是定态Schrödinger方程中的能量本征值，对应的是总Hamilton量的*连续谱*。利用散射态波函数在$abs(vb(r))->oo$时的渐进行为得到*散射振幅*。
- *含时散射图像*：$t -> -oo$时粒子波包从无穷远处入射，进入相互作用区域后在势场$V(vb(r))$的作用下进行时间演化，最后从各个方向出射，从而被探测器探测到。利用含时演化理论计算散射振幅。
我们将证明，两种物理图像是等价的。

== 定态散射理论

最广义的散射现象是两个粒子发生碰撞的现象，两个粒子都可以具有内部结构和内部状态 (例如自旋)。由于两者之间的相互作用，它们的运动状态都要发生变化。
- *弹性散射*：若碰撞前后两粒子的内部结构和内部状态均无变化，则称为弹性散射。
- *非弹性散射*：若碰撞前后两粒子的内部状态发生变化，甚至内部结构发生了变化 (反应)，则称为非弹性散射。
这里考虑最简单的弹性散射：两个粒子之间的相互作用可以用相互作用势$V(vb(r))$来描述，
$
  V(vb(r)) = V(vb(r)_1 - vb(r)_2)
$
称为*势散射*。将两体问题约化为单体问题，势散射的图像就成为：入射粒子在*散射中心*的作用下运动状态发生改变。

*势散射的物理图像*：如图所示，O是散射中心和坐标原点，在O附近存在势场$V(vb(r))$。动量为 $vb(p)$、能量为$E = vb(p)^2/(2m)$的粒子源源不断地由从离O很远的狭缝入射，然后实验上在远离O的各个方向上对粒子进行探测。当粒子源源不断地入射时，空间各点的概率密度和概率流密度不随时间变化，从而可以采用定态描述。
#figure(
  image("pic/2025-12-07-02-05-13.png", width: 50%),
  numbering: none,
)
系统的Hamilton量为
$
  hat(H) = hat(vb(p))^2/(2m) + V(vb(r))
$
由于势场$V(vb(r))$不含时，所以系统存在定态。在坐标表象中，定态方程$hat(H) ket(psi) = E ket(psi)$写为
$
  (- hbar^2/(2m) laplacian + V(vb(r))) psi(vb(r)) = E psi(vb(r)), E = vb(p)^2/(2m)
$
我们感兴趣的显然是$E > 0$的*散射态解*，所以很自然地要求Hamilton量在$0 < E < oo$区间存在*连续谱*，所以对$V(vb(r))$最基本的要求是$(r ≡ abs(vb(r)))$
$
  V(vb(r)) -> 0, r -> oo
$

#newpara()

实际上，在量子势散射理论中，对势场的性质还有更高的要求：
1. 势场$V(vb(r))$应当足够光滑；
2. $r -> oo$时势场趋于零，而且趋于零的速度要足够快，Coulomb势不满足这个条件，需要特殊考虑；
3. $r -> 0$时势场的奇异性不能太大。
具体细节由含时散射理论的*渐近条件*决定，此处暂时不做讨论。

=== 定态散射态的波函数与散射振幅

在散射中心O附近取一个区域，其边界为C，在C之外$V(vb(r)) = 0$。由于势场趋于零的速度足够快，所以这是可以做到的。C通常在宏观上不大，但是在微观上足够大，以致于其边界处就可以认为对于O点是无穷远，在边界外面势场$V(vb(r))$不起作用。所以，在C之外(即 $r -> oo$)，定态方程成为
$
  - hbar^2/(2m) laplacian psi(vb(r)) = E psi(vb(r)), E = vb(p)^2/(2m)
$
令$vb(p) = hbar vb(k)$，则$r -> oo$处的定态方程为
$
  laplacian psi(vb(r)) + k^2 psi(vb(r)) = 0
$
不失一般性，取粒子入射方向为$z$方向，则方程的一个解为
$
  psi(vb(r)) prop e^(i vb(k) dot vb(r)) = e^(i k z)
$
这正是入射平面波。由于角动量守恒，另外的解可取为
$
  psi(r, theta, phi) = (u_l (r))/r Y_(l m) (theta, phi)
$
其中球谐函数$Y_(l m) (theta, phi)$是轨道角动量算符$hat(vb(L))^2$和$hat(L)_z$的共同本征函数。代入定态方程，得到径向方程
$
  dv(, r, 2) u_l (r) + (k^2 - l(l+1)/r^2) u_l (r) = 0
$
当$r -> oo$时，$u_l (r)$具有下列形式(与$l$无关)
$
  u_l (r) prop e^(plus.minus i k r)
$
显然，符合散射问题的解应取正号。所以，在$r -> oo$处定态方程的一般解为
$
  psi(vb(r)) & = A e^(i k z) + sum_(l m) c_(l m) e^(i k r)/r Y_(l m) (theta, phi) \
             & = A(e^(i k z) + f(theta, phi) e^(i k r)/r)
$
其中$A$为常数(与采用的归一化有关)，$f(theta, phi)$称为*散射振幅*，它是散射理论计算的目标。要计算散射振幅，通常需要求解边界C以内的定态方程，并利用边界条件。

=== 散射截面

*散射振幅如何与物理观测量联系起来？*上式的物理图像是：从无穷远处沿$z$方向入射的平面波，被势场散射后，在各个方向的无穷远处产生了出射球面波。若粒子的波函数为$psi(vb(r))$，则对应的*概率密度*为
$
  rho = psi^* psi
$
*概率流密度*为
$
  vb(J) = - (- i hbar)/(2m) (psi^* grad psi - psi grad psi^*) = hbar/m Im(psi^* grad psi)
$
满足连续性方程
$
  pdv(rho, t) + div vb(J) = 0
$
根据前面的讨论，在边界C以外，波函数可以写为
$
  psi(vb(r)) & = psi_I (vb(r)) + psi_S (vb(r)) \
$
其中$psi_I$和$psi_S$分别代表入射波和出射波
$
  psi_I (vb(r)) = A e^(i k z) \
  psi_S (vb(r)) = A f(theta, phi) e^(i k r)/r \
$
带入概率流密度公式得到
$
  vb(J) = hbar/m Im(psi^*_I grad psi_I) + hbar/m Im(psi^*_S grad psi_S) + hbar/m Im(psi^*_I grad psi_S) + hbar/m Im(psi^*_S grad psi_I)
$
后两项代表入射波与出射波的干涉项。

可以证明，当$r -> oo$时，除了非常靠近入射方向$(θ approx 0)$之外，干涉效应可以忽略。这样，我们可以把概率流密度写为入射流密度和出射流密度之和
$
  vb(J) = vb(J)_I + vb(J)_S
$
其中，入射流密度可计算为
$
  vb(J)_I = hbar/m Im(psi^*_I grad psi_I) = abs(A)^2 (hbar k)/m vu(e)_z
$
出射流密度可计算为
$
  vb(J)_S & = hbar/m Im(psi^*_S grad psi_S) \
          & = abs(A)^2 (hbar k)/m abs(f(theta, phi))^2/r^2 vu(e)_r + O(1/r^3) vu(e)_theta + O(1/r^3) vu(e)_phi \
$
在$(theta, phi)$方向上的*微分散射截面*$dv(sigma, Omega)$的定义是：在远处$(theta, phi)$方向上单位时间进入单位立体角中的粒子数$dv(N, Omega)$与入射粒子流密度之比。

根据出射流密度计算结果，在$(θ, ϕ)$方向上的立体角元$dd(Ω)$中单位时间出射的粒子数为$(dd(vb(S)) = r^2dd(Ω)vu(e)_r)$
$
  dd(N) = vb(J)_S dot dd(vb(S)) = abs(A)^2 (hbar k)/m abs(f(theta, phi))^2 dd(Ω)
$
微分散射截面为
$
  dv(sigma, Omega) = 1/(vb(J)_I) dv(N, Omega) = abs(f(theta, phi))^2
$
对立体角积分可得到总散射截面$sigma$。

总结：与以往不同，在散射问题中，我们要与Hamilton量的连续谱即散射态打交道。在*定态散射理论*中，我们实际上是求解*以入射粒子能量为能量本征值*的能量本征态
$
  (hat(vb(p))^2/(2m) + V(vb(r))) ket(psi) = E ket(psi), E = vb(p)^2/(2m) > 0
$
在物理上，这意味着粒子的能量没有变化，即*弹性散射*。
- *散射态问题*：在能量本征值$E$已知的情况下，求解本征态$ket(psi)$。这里的能量本征值在Hamilton量的连续谱范围内。对于满足条件的势场$V(vb(r))$，即 $0 < E < oo$。
- *束缚态问题*：与之对比，束缚态问题则是要同时求解离散的能量本征值$E_n$和本征态$ket(n)$。对于满足条件的势场$V(vb(r))$，仍然可能存在束缚态，即 $E_n < 0$。

=== 定态散射理论

考虑粒子在势场$V(vb(r))$中运动，Hamilton量为
$
  hat(H) = hat(H)_0 + hat(V)(vb(r)), hat(H)_0 = hat(vb(p))^2/(2m)
$
在前面对散射问题的描述中，我们将*入射态取为动量本征态$ket(vb(p))$，它也是$hat(H)_0$的本征态*
$
  hat(H)_0 ket(vb(p)) = vb(p)^2/(2m) ket(vb(p))
$
其*能量本征值就是入射能量*。在散射理论中，通常令$vb(p) = hbar vb(k)$，将动量本征态$ket(vb(p))$改写为$ket(vb(k))$，即
$
  hat(H)_0 ket(vb(k)) = (hbar^2 vb(k)^2)/(2m) ket(vb(k))
$
态矢$ket(vb(k))$可定义为(不唯一)
$
  ket(vb(k)) = hbar^(3/2) ket(vb(p)) \
$
态矢$ket(vb(k))$的归一化
$
  braket(vb(p)', vb(p)) = delta(vb(p)' - vb(p)) => braket(vb(k)', vb(k)) = delta(vb(k)' - vb(k))
$
态矢$ket(vb(k))$在坐标表象的波函数
$
  braket(vb(r), vb(p)) = 1/(2 pi hbar)^(3/2) e^(i/hbar vb(p) dot vb(r)) => braket(vb(r), vb(k)) = 1/(2 pi)^(3/2) e^(i vb(k) dot vb(r))
$
完备性条件
$
  integral dd(vb(p), 3) ketbra(vb(p)) = 1 => integral dd(vb(k), 3) ketbra(vb(k)) = 1
$
#newpara()
我们要求解的定态方程为
$
  (hat(H)_0 + hat(V)(vb(r))) ket(psi) = E_vb(k) ket(psi), E_vb(k) = (hbar^2 vb(k)^2)/(2m)\
  => (E_vb(k) - hat(H)_0) ket(psi) = hat(V)(vb(r)) ket(psi) \
$
再利用
$
  (E_vb(k) - hat(H)_0) ket(vb(k)) = 0
$
两式相减得到
$
  (E_vb(k) - hat(H)_0)(ket(psi) - ket(vb(k))) = hat(V) ket(psi)
$
两边从左边作用$(E_vb(k) - hat(H)_0)^(-1)$
$
  ket(psi) = ket(vb(k)) + 1/(E_vb(k) - hat(H)_0) hat(V) ket(psi)
$
这里出现的一个问题是逆算符$(E_vb(k) - hat(H)_0)^(-1)$具有奇异性，它作用在$ket(vb(k))$得到无穷大。解决的办法是将能量$E_vb(k)$扩大到复数域，改写为
$
  ket(psi^(plus.minus)) = ket(vb(k)) + 1/(E_vb(k) plus.minus i epsilon - hat(H)_0) hat(V) ket(psi^(plus.minus))
$
其中$epsilon$是一个无穷小正实数，其来源在含时散射理论中更清楚。这个方程就是*Lippmann-Schwinger方程*。定义算符
$
  hat(G)_0^(plus.minus) (E) = 1/(E - hat(H)_0 plus.minus i epsilon)
$
称为*自由Green算符*(传播子)。利用$ket(vb(k))$的完备性得到
$
  hat(G)_0^(plus.minus) (E) = integral dd(vb(k)', 3) ketbra(vb(k)', vb(k))/(E - E_(vb(k)) plus.minus i epsilon)
$
以上式子中正负号$plus.minus$的意义尚不明确。为此，让Lippmann-Schwinger方程进入坐标表象，得到
$
  braket(vb(r), psi^(plus.minus)) & = braket(vb(r), vb(k)) + integral dd(vb(r)', 3) braket(vb(r), hat(G)_0^(plus.minus) (E_vb(k)), vb(r)') braket(vb(r)', hat(V), psi^(plus.minus)) \
$
计算自由格林算符的矩阵元(为方便乘以常数$hbar^2/(2m)$
$
  cal(G)_0^(plus.minus) (vb(r), vb(r)') & = hbar^2/(2m) braket(vb(r), hat(G)_0^(plus.minus) (E_vb(k)), vb(r)') \
  & = hbar^2/(2m) integral dd(vb(k)', 3) integral dd(vb(k)'', 3) braket(vb(r), vb(k)') braket(vb(k)', hat(G)_0^(plus.minus) (E_vb(k)), vb(k)'') braket(vb(k)'', vb(r)') \
$
其中易证
$
  braket(vb(k)', hat(G)_0^(plus.minus) (E), vb(k)'') &= integral dd(vb(k)''', 3) braket(vb(k)', vb(k)''') 1/(E - E_(vb(k)''') plus.minus i epsilon) braket(vb(k)''', vb(k)'') \
  & = delta(vb(k)' - vb(k)'')/(E - E_(vb(k)') plus.minus i epsilon) \
$
再将
$
  braket(vb(r), vb(k)') = 1/(2 pi)^(3/2) e^(i vb(k)' dot vb(r)), \
  braket(vb(k)'', vb(r)') = 1/(2 pi)^(3/2) e^(- i vb(k)'' dot vb(r)')
$
带入得到$eta = (2 m epsilon)/(hbar^2)$
$
  cal(G)_0^(plus.minus) (vb(r), vb(r)') & = hbar^2/(2m) integral dd(vb(k)', 3)/(2 pi)^3 e^(i vb(k)' dot (vb(r) - vb(r)'))/(E_vb(k) - E_(vb(k)') plus.minus i epsilon) \
  &= integral dd(vb(k)', 3)/(2 pi)^3 e^(i vb(k)' dot (vb(r) - vb(r)'))/(vb(k)^2 - vb(k)'^2 plus.minus i eta) \
$
利用球坐标系，先将角度积掉，得到
$
  cal(G)_0^(plus.minus) (vb(r), vb(r)') & = i/(8 pi^2) 1/abs(vb(r) - vb(r)') integral_(-oo)^(oo) dd(k)' k' (e^(- i k' abs(vb(r) - vb(r)')) - e^(i k' abs(vb(r) - vb(r)')))/(k^2 - k'^2 plus.minus i eta) \
$
利用留数定理计算积分
$
  cal(I) = integral_(-oo)^(oo) dd(k)' k' e^(plus.minus i k' abs(vb(r) - vb(r)'))/(k^2 - k'^2 plus.minus i eta) \
$
#figure(
  image("pic/2025-12-07-15-30-12.png", width: 30%),
  numbering: none,
)
被积函数的奇点位于
$
  k'^2 = k^2 plus.minus i eta => k' = k plus.minus i delta, k' = - k minus.plus i delta
$
第一项含有$e^(-k' abs(vb(r) - vb(r)'))$。此项积分等于*下半平面*半圆围道上的积分，因为下半平面无穷远处的半圆弧对积分贡献为零(Jordan引理)。这样，可以利用留数定理计算这个围道积分，结果为
$
  (-1)2 pi i (minus.plus) (e^(-i(minus.plus k) abs(vb(r) - vb(r)'))/(plus.minus 2 k)) = i pi e^(plus.minus i k abs(vb(r) - vb(r)'))
$
第二项含有$e^(i k' abs(vb(r) - vb(r)'))$。此项积分等于*上半平面*半圆围道上的积分，因为上半平面无穷远处的半圆弧对积分贡献为零。结果与第一项相同 (前面的负号也包含在内)。
$
  cal(I) = i pi e^(plus.minus i k abs(vb(r) - vb(r)'))
$
最终计算得到
$
  cal(G)_0^(plus.minus) (vb(r), vb(r)') = - 1/(4 pi) e^(plus.minus i k abs(vb(r) - vb(r)'))/abs(vb(r) - vb(r)')
$
它实际上是Helmholtz方程的Green函数，即
$
  (laplacian + k^2) cal(G)_0^(plus.minus) (vb(r), vb(r)') = delta(vb(r) - vb(r)')
$
带入Lippmann-Schwinger方程
$
  braket(vb(r), psi^(plus.minus)) & = braket(vb(r), vb(k)) - (2m)/(hbar^2) integral dd(vb(r)', 3) e^(plus.minus i k abs(vb(r) - vb(r)'))/(4 pi abs(vb(r) - vb(r)')) braket(vb(r)', hat(V), psi^(plus.minus)) \
$
我们预期，右边第二项当$abs(vb(r)) -> oo$时给出散射波的渐近形式。

由于$hat(V) = V(hat(vb(r)))$是定域势，所以
$
  braket(vb(r)', hat(V), vb(r)'') = braket(vb(r)', V(hat(vb(r))), vb(r)'') = V(vb(r)') delta(vb(r)' - vb(r)'')
$
得到
$
  braket(vb(r), hat(V), psi^(plus.minus)) = integral dd(vb(r)'', 3) braket(vb(r), ', hat(V), vb(r)'') braket(vb(r)'', psi^(plus.minus)) = V(vb(r)') braket(vb(r)', psi^(plus.minus))
$
这样
$
  braket(vb(r), psi^(plus.minus)) & = braket(vb(r), vb(k)) - (2m)/(hbar^2) integral dd(vb(r)', 3) e^(plus.minus i k abs(vb(r) - vb(r)'))/(4 pi abs(vb(r) - vb(r)')) V(vb(r)') braket(vb(r)', psi^(plus.minus)) \
$
这是关于波函数
$
  braket(vb(r), psi^(plus.minus)) = psi^(plus.minus)(vb(r)) \
$
的*积分方程*。

对于在远处迅速趋于零的*短程势*，右边第二项对$r'$的积分只集中在有限区域(即边界C以内)。考虑$abs(vb(r)) -> oo$，可以认为
$
  abs(vb(r)) >> abs(vb(r)') \
$
设$vb(r)$和$vb(r)'$夹角为$alpha$。当$r>>r'$有
$
  abs(vb(r) - vb(r)') & = sqrt(r^2 + r'^2 - 2 r r' cos alpha) \
                      & = r (1 - (2r')/r cos alpha + (r'^2)/(r^2)) \
                      & approx r - vu(r) dot vb(r)' \
$
其中$vu(r) = vb(r)/r$为$r$方向的单位矢量。定义$vb(k)' = k vu(r)$则
$
  e^(plus.minus i k abs(vb(r) - vb(r)')) & = e^(plus.minus i k r) e^(minus.plus i vb(k)' dot vb(r)') \
$
最终我们得到$r -> oo$时波函数的渐进行为
$
  braket(vb(r), psi^(plus.minus)) -> braket(vb(r), vb(k)) - 1/(4 pi) (2m)/(hbar^2) e^(plus.minus i k r)/(r) integral dd(vb(r)', 3) e^(minus.plus i vb(k)' dot vb(r)') V(vb(r)') braket(vb(r)', psi^(plus.minus)) \
$
显然，对于真实的散射问题，应取$ket(psi^(+))$，因为它对应的是入射波加上出射波的物理图像。所以，散射态波函数的渐进行为为
$
  psi^(+)(vb(r)) &-> braket(vb(r), vb(k)) - 1/(4 pi) (2m)/(hbar^2) e^(i k r)/(r) integral dd(vb(r)', 3) e^(i vb(k)' dot vb(r)') V(vb(r)') braket(vb(r)', psi^+) \
  &= 1/(2 pi)^(3/2) (e^(i vb(k) dot vb(r)) + f(vb(k)', vb(k)) e^(i k r)/r) \
$
散射振幅为
$
  f(vb(k)', vb(k)) = - 1/(4 pi) (2m)/(hbar^2) (2 pi)^3 integral dd(vb(r)', 3) e^(i vb(k)' dot vb(r)')/(2 pi)^(3/2) V(vb(r)') braket(vb(r)', psi^+) \
$
利用
$
  1/(2 pi)^(3/2) e^(- i vb(k)' dot vb(r)') = braket(vb(k)', vb(r)')
$
可将*散射振幅*写为
$
  f(vb(k)', vb(k)) &= - 1/(4pi) (2m)/(hbar^2) (2 pi)^(3/2) integral dd(vb(r)', 3) braket(vb(k)', vb(r)') V(vb(r)') braket(vb(r)', psi^+) \
  &= - (4 pi^2 m)/(hbar^2) integral dd(vb(r)', 3) integral dd(vb(k)'', 3) braket(vb(k)', vb(r)') braket(vb(k)', V, vb(k)'') braket(vb(k)'', psi^+) \
  &= - (4 pi^2 m)/(hbar^2) braket(vb(k)', V, psi^+)
$
利用此式计算散射振幅，需要知道波函数$psi^+(vb(r))$或本征态$ket(psi^+)$的解，*因此实际上只有形式上的意义*。

进一步把散射振幅写成某个算符的矩阵元。为此我们定义*T算符*
$
  hat(T) ket(vb(k)) = hat(V) ket(psi^+)
$
由于$ket(vb(k))$是一组完备的基矢，此式完全定义了T算符。这样，散射振幅可以写为
$
  f(vb(k)', vb(k)) = - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(T), vb(k))
$
将Lippmann-Schwinger方程两边同时作用$hat(V)$得到
$
  hat(V) ket(psi^+) = hat(V) ket(vb(k)) + hat(V) 1/(E_vb(k) + i epsilon - hat(H)_0) hat(V) ket(psi^+)\
  => hat(T) ket(vb(k)) = hat(V) ket(vb(k)) + hat(V) 1/(E_vb(k) + i epsilon - hat(H)_0) hat(T) ket(vb(k))
$
由于$ket(vb(k))$是一组完备的基矢，因此T算符满足方程
$
  hat(T) = hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(T)
$
利用迭代法，可将T算符按照$hat(V)$展开
$
  hat(T) &= hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(V) + ...\
  &= hat(V) sum_(n=0)^oo (hat(G)_0^+ (E_vb(k)) hat(V))^n \
$
此即*T算符的Born级数形式*。利用这个级数形式，可以将散射振幅按照$V$的强度进行逐级近似计算。

一阶近似 (Born 近似) 的结果计算如下
$
  f^((1)) (vb(k)', vb(k)) & = - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(V), vb(k)) \
  & = - (4 pi^2 m)/(hbar^2) integral dd(vb(r)', 3) integral dd(vb(r)'', 3) braket(vb(k)', vb(r)') braket(vb(r)', hat(V), vb(r)'') braket(vb(r)'', vb(k)) \
  & = - (4 pi^2 m)/(hbar^2) integral dd(vb(r)', 3) (e^(- i vb(k)' dot vb(r)'))/(2 pi)^(3/2) V(vb(r)') (e^(i vb(k) dot vb(r)'))/(2 pi)^(3/2) \
  & = - (m)/(2 pi hbar^2) integral dd(vb(r)', 3) e^(i (vb(k) - vb(k)') dot vb(r)') V(vb(r)') \
$
更高阶近似的结果，按照的级数展开继续计算即可。二阶近似的结果为
$
  f^((2)) (vb(k)', vb(k)) & = - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(V) 1/(E_vb(k) + i epsilon - hat(H)_0) hat(V), vb(k)) \
  & = - (4 pi^2 m)/(hbar^2) integral dd(vb(r)', 3) integral dd(vb(r)'', 3) integral dd(vb(r)_1, 3) integral dd(vb(r)_2, 3) \
  & braket(vb(k)', vb(r)') braket(vb(r)', hat(V), vb(r)_1) braket(vb(r)_1, 1/(E_vb(k) + i epsilon - hat(H)_0), vb(r)_2) braket(vb(r)_2, hat(V), vb(r)'') braket(vb(r)'', vb(k)) \
  & = - (m)/(2 pi hbar^2) integral dd(vb(r)', 3) integral dd(vb(r)'', 3) e^(i vb(k) dot vb(r)) V(vb(r)') (2m)/hbar^2 cal(G)_0^+ (vb(r)', vb(r)'') V(vb(r)'') e^(- i vb(k) dot vb(r)'') \
$

=== 定态散射理论的另一种理论框架

要求解的定态方程为
$
  hat(H) ket(psi) = E_vb(k) ket(psi), E_vb(k) = (hbar^2 vb(k)^2)/(2m)\
  => (E_vb(k) - hat(H)) ket(psi) = 0 \
$
同时，对于动量本征态$ket(vb(k))$有
$
  (E_vb(k) - hat(H)_0) ket(vb(k)) = 0 \
  => (E_vb(k) - hat(H)) ket(vb(k)) + hat(V) ket(vb(k)) = 0 \
$
两式相减得到
$
  (E_vb(k) - hat(H))(ket(psi) - ket(vb(k))) - hat(V) ket(vb(k)) = 0\
$
两边同时作用$(E_vb(k) - hat(H))^(-1)$，得到
$
  ket(psi) = ket(vb(k)) - 1/(E_vb(k) - hat(H)) hat(V) ket(vb(k)) \
$
同样，引入$plus.minus i epsilon$，上式成为
$
  ket(psi^(plus.minus)) = ket(vb(k)) - 1/(E_vb(k) - hat(H) plus.minus i epsilon) hat(V) ket(vb(k)) \
$
定义算符
$
  hat(G)^(plus.minus) (E) = 1/(E - hat(H) plus.minus i epsilon)
$
称为*完全Green算符*(传播子)。上式虽然在形式上给出了$ket(psi^(plus.minus))$的解，但是实际上是把难度扔给了完全Green算符。可以证明，完全Green算符和自由Green算符之间的关系为
#theorem(subname: [Dyson-Schwinger方程])[
  Green函数满足的 Dyson-Schwinger 方程
  $
    hat(G)^(plus.minus) (E) = hat(G)_0^(plus.minus) (E) + hat(G)_0^(plus.minus) (E) hat(V) hat(G)^(plus.minus) (E)
  $
]
#proof[
  $
       & E - hat(H) plus.minus i epsilon = E - hat(H)_0 - hat(V) plus.minus i epsilon \
    => & (hat(G)^(plus.minus) (E))^(-1) = (hat(G)_0^(plus.minus) (E))^(-1) - hat(V) \
    => & (hat(G)^(plus.minus) (E))^(-1) = (hat(G)_0^(plus.minus) (E))^(-1) (1 - hat(G)_0^(plus.minus) (E) hat(V)) \
    => & hat(G)^(plus.minus) (E) = (1 - hat(G)_0^(plus.minus) (E) hat(V))^(-1) hat(G)_0^(plus.minus) (E) \
    => & (1 - hat(G)_0^(plus.minus) (E) hat(V)) hat(G)^(plus.minus) (E) = hat(G)_0^(plus.minus) (E) \
    => & hat(G)^(plus.minus) (E) = hat(G)_0^(plus.minus) (E) hat(V) hat(G)^(plus.minus) (E) + hat(G)_0^(plus.minus) (E) \
  $
]
两边同时用$hat(V)$作用，得到
$
  hat(T) = hat(V) + hat(V) hat(G)^(plus) (E_vb(k)) hat(V)
$
利用 Dyson-Schwinger 方程进行迭代，就得到同样的 Born 级数。
