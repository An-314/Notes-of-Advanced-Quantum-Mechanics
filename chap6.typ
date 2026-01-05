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
- *非弹性散射*：若碰撞前后两粒子的内部状态发生变化，甚至内部结构发生了变化(反应)，则称为非弹性散射。
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
2. $r -> oo$时势场趋于零，而且趋于零的速度要足够快，*Coulomb势不满足这个条件*，需要特殊考虑；
3. $r -> 0$时势场的奇异性不能太大。
具体细节由含时散射理论的*渐近条件*决定，此处暂时不做讨论。

=== 定态散射态的波函数与散射振幅

在散射中心O附近取一个区域，其边界为C，在C之外$V(vb(r)) = 0$。由于势场趋于零的速度足够快，所以这是可以做到的。C通常在宏观上不大，但是在微观上足够大，以致于其边界处就可以认为对于O点是无穷远，在边界外面势场$V(vb(r))$不起作用。所以，在C之外(即$r -> oo$)，定态方程成为
$
  - hbar^2/(2m) laplacian psi(vb(r)) = E psi(vb(r)), E = vb(p)^2/(2m)
$
令$vb(p) = hbar vb(k)$，则$r -> oo$处的定态方程为
$
  laplacian psi(vb(r)) + vb(k)^2 psi(vb(r)) = 0
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
其中$A$为常数(与采用的归一化有关)，$f(theta, phi)$称为*散射振幅*，它是*散射理论计算的目标*。要计算散射振幅，通常需要求解边界C以内的定态方程，并利用边界条件。

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

#proposition(subname: [散射截面与散射振幅])[
  在量子势散射理论中，散射振幅描述了入射粒子波函数被势场散射后的出射波函数的角分布特征
  $
    psi(vb(r)) = A(e^(i k z) + f(theta, phi) e^(i k r)/r)
  $
  散射振幅$f(theta, phi)$与微分散射截面$dv(sigma, Omega)$的关系为
  $
    dv(sigma, Omega) = abs(f(theta, phi))^2
  $
]
#newpara()

总结：与以往不同，在散射问题中，我们要与Hamilton量的连续谱即散射态打交道。在*定态散射理论*中，我们实际上是求解*以入射粒子能量为能量本征值*的能量本征态
$
  (hat(vb(p))^2/(2m) + V(vb(r))) ket(psi) = E ket(psi), E = vb(p)^2/(2m) > 0
$
在物理上，这意味着粒子的能量没有变化，即*弹性散射*。
- *散射态问题*：在能量本征值$E$已知的情况下，求解本征态$ket(psi)$。这里的能量本征值在Hamilton量的连续谱范围内。对于满足条件的势场$V(vb(r))$，即 $0 < E < oo$。
- *束缚态问题*：与之对比，束缚态问题则是要同时求解离散的能量本征值$E_n$和本征态$ket(n)$。对于满足条件的势场$V(vb(r))$，仍然可能存在束缚态，即 $E_n < 0$。

=== 定态散射理论

*考虑粒子在势场$V(vb(r))$中运动，Hamilton量为*
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
*态矢$ket(vb(k))$可定义*为(不唯一)
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
我们要求解的*定态方程*为
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
其中$epsilon$是一个无穷小正实数，其来源在含时散射理论中更清楚。这个方程就是*Lippmann-Schwinger方程*。

#theorem(subname: [Lippmann-Schwinger方程])[
  散射态$ket(psi^(plus.minus))$满足Lippmann-Schwinger方程
  $
    ket(psi^(plus.minus)) = ket(vb(k)) + 1/(E_vb(k) plus.minus i epsilon - hat(H)_0) hat(V) ket(psi^(plus.minus))
  $
  其中$epsilon$是一个无穷小正实数
  $
    ket(vb(k)) = hbar^(3/2) ket(vb(p)), hat(H)_0 ket(vb(k)) = (hbar^2 vb(k)^2)/(2m) ket(vb(k))
  $
  是入射态，$E_vb(k) = (hbar^2 vb(k)^2)/(2m)$为入射能量
]
#newpara()

定义算符
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
计算自由Green算符的矩阵元(为方便乘以常数$hbar^2/(2m)$)
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
  braket(vb(r), vb(k)') = 1/(2 pi)^(3/2) e^(i vb(k)' dot vb(r)), braket(vb(k)'', vb(r)') = 1/(2 pi)^(3/2) e^(- i vb(k)'' dot vb(r)')
$
带入得到($eta = (2 m epsilon)/(hbar^2)$)
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
  braket(vb(r), hat(V), psi^(plus.minus)) = integral dd(vb(r)'', 3) braket(vb(r)', hat(V), vb(r)'') braket(vb(r)'', psi^(plus.minus)) = V(vb(r)') braket(vb(r)', psi^(plus.minus))
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

#proposition(subname: [散射态波函数的积分方程])[
  散射态波函数$psi^(plus.minus)(vb(r))$满足Lippmann-Schwinger积分方程
  $
    psi^(plus.minus)(vb(r)) = braket(vb(r), vb(k)) - (2m)/(hbar^2) integral dd(vb(r)', 3) e^(plus.minus i k abs(vb(r) - vb(r)'))/(4 pi abs(vb(r) - vb(r)')) V(vb(r)') psi^(plus.minus)(vb(r)')
  $
]
#newpara()

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
显然，对于真实的散射问题，应取$ket(psi^(+))$，因为它对应的是入射波加上出射波的物理图像。所以，*散射态波函数的渐进行为*为
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
  &= - (4 pi^2 m)/(hbar^2) integral dd(vb(r)', 3) integral dd(vb(r)'', 3) braket(vb(k)', vb(r)') braket(vb(r)', hat(V), vb(r)'') braket(vb(r)'', psi^+) \
  &= - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(V), psi^+)
$
利用此式计算散射振幅，需要知道波函数$psi^+(vb(r))$或本征态$ket(psi^+)$的解，*因此实际上只有形式上的意义*。

#proposition(subname: [散射振幅的表达式])[
  散射振幅$f(vb(k)', vb(k))$可表示为
  $
    f(vb(k)', vb(k)) = - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(V), psi^+)
  $
  其中$ket(psi^+)$为Lippmann-Schwinger方程的解。
]
#newpara()

=== T算符与Born级数展开

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

#proposition(subname: [T算符的Born级数形式])[
  T算符定义为
  $
    hat(T) ket(vb(k)) = hat(V) ket(psi^+) => hat(T) = hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(T)
  $
  T算符的Born级数形式为
  $
    hat(T) = hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(V) + ...\
    = hat(V) sum_(n=0)^oo (hat(G)_0^+ (E_vb(k)) hat(V))^n
  $
]
#newpara()

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

=== 光学定理

#theorem(subname: [光学定理])[
  总散射截面$sigma_t$与向前散射振幅的虚部存在如下关系
  $
    sigma_t = (4 pi)/k Im f(vb(k), vb(k)) = (4 pi)/k Im f(0)
  $
]

#proof[
  将 Lippmann-Schwinger 方程的共轭方程
  $
    bra(vb(k)) = bra(psi^+) - bra(psi^+) hat(V) 1/(E_vb(k) - hat(H)_0 - i epsilon)
  $
  与$hat(V) ket(psi^+)$ 做内积得到
  $
    braket(vb(k), hat(V), psi^+) = braket(psi^+, hat(V), psi^+) - braket(psi^+, hat(V) 1/(E_vb(k) - hat(H)_0 - i epsilon) hat(V), psi^+)
  $
  散射振幅为
  $
    f(vb(k)', vb(k)) = - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(V), psi^+) = - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(T), vb(k))
  $
  考虑$f(vb(k), vb(k))$的虚部
  $
    Im braket(vb(k), hat(T), vb(k)) = Im braket(vb(k), hat(V), psi^+) = Im braket(psi^+, hat(V), psi^+) - Im braket(psi^+, hat(V) 1/(E_vb(k) - hat(H)_0 - i epsilon) hat(V), psi^+)
  $
  由于$hat(V)$是Hermitian算符，所以第一项为零
  $
    Im braket(vb(k), hat(T), vb(k)) &= - Im braket(psi^+, hat(V) 1/(E_vb(k) - hat(H)_0 - i epsilon) hat(V), psi^+)\
    &= - Im braket(vb(k), hat(T)^dagger 1/(E_vb(k) - hat(H)_0 - i epsilon) hat(T), vb(k))\
    &= - Im integral dd(vb(q), 3) braket(vb(k), hat(T)^dagger 1/(E_vb(k) - hat(H)_0 - i epsilon), vb(q)) braket(vb(q), hat(T), vb(k)) \
    &= - Im integral dd(vb(q), 3) braket(vb(k), hat(T)^dagger, vb(q)) braket(vb(q), hat(T), vb(k)) 1/(E_vb(k) - E_vb(q) - i epsilon)\
  $
  利用Dirac恒等式
  $
    1/(x plus.minus i epsilon) = "P" 1/x minus.plus i pi delta(x)
  $
  得到
  $
    Im braket(vb(k), hat(T), vb(k)) &= - pi integral dd(vb(q), 3) abs(braket(vb(k), hat(T), vb(q)))^2 delta(E_vb(k) - E_vb(q))
  $
  采用球坐标，利用
  $
    delta(E_vb(k) - E_vb(q)) = m/(hbar^2 k) delta(k - q), dd(vb(q), 3) = q^2 dd(q) dd(Omega_vb(q))
  $
  完成对$q$的积分，得到$(abs(vb(k)') = abs(vb(k)) = k)$
  $
    Im braket(vb(k), hat(T), vb(k)) & = - (m pi k)/(hbar^2) integral dd(Omega_vb(k)) abs(braket(vb(k)', hat(T), vb(k)))^2 \
  $
  根据散射振幅与$T$矩阵元的关系
  $
    f(vb(k)', vb(k)) = - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(T), vb(k))
  $
  得到
  $
    Im f(vb(k), vb(k)) & = - (4 pi^2 m)/hbar^2 Im braket(vb(k), hat(T), vb(k)) \
    & = - (4 pi^2 m)/hbar^2 (- (pi m k)/(hbar^2) integral dd(Omega_vb(k)) abs(braket(vb(k)', hat(T), vb(k)))^2) \
    & = (k)/(4 pi) integral dd(Omega_vb(k)) abs(f(vb(k)', vb(k)))^2 \
    & = k/(4 pi) sigma_t
  $
]

== 含时散射理论

喀兴林《高等量子力学》第六章 28 小节，或者孙昌璞《量子力学现代教程》 6.4 小节

== 分波展开

如果势场$V(vb(r))$是球对称的，即$V(vb(r)) = V(r)$，根据方程
$
  hat(T) = hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(T)
$
可知$T$算符与$hat(vb(L))^2$和$hat(L)_z$对易。对于自由粒子部分$hat(H)_0$，角动量亦是守恒量，可将力学量完备集取为${hat(H)_0, hat(vb(L))^2, hat(L)_z}$，共同本征态记为$ket(E l m)$，即
$
      hat(H)_0 ket(E l m) & = E ket(E l m), \
  hat(vb(L))^2 ket(E l m) & = l(l+1) hbar^2 ket(E l m), \
      hat(L)_z ket(E l m) & = m hbar ket(E l m)
$
这些本征态的归一化采用如下的 convention
$
  braket(E' l' m', E l m) = delta(E' - E) delta_(l l') delta_(m m')
$
完备性条件为
$
  sum_(l=0)^oo sum_(m=-l)^(l) integral dd(E) ketbra(E l m) = 1
$
可将$ket(vb(k))$态展开为
$
  ket(vb(k)) = sum_(l=0)^oo sum_(m=-l)^(l) integral dd(E) ket(E l m) braket(E l m, vb(k)) \
$
与$vb(k)$表象和$vb(r)$表象之间的变换矩阵元为 (Sakurai, pp.407)
$
  braket(vb(k), E l m) & = hbar/sqrt(m k) delta(E - E_vb(k)) Y_(l m) (vu(k)), E_vb(k) = (hbar^2 vb(k)^2)/(2m) \
  braket(vb(r), E l m) & = i^l/hbar sqrt((2 m k)/pi) j_l (k r) Y_(l m) (vu(r)) \
$
考虑散射振幅，利用$E l m$表象的完备性得到
$
  f(vb(k)', vb(k)) &= - (4 pi^2 m)/(hbar^2) braket(vb(k)', hat(T), vb(k)) \
  &= - (4 pi^2 m)/(hbar^2) sum_l sum_m sum_l' sum_m' integral dd(E) integral dd(E') braket(vb(k)', E' l' m') braket(E' l' m', hat(T), E l m) braket(E l m, vb(k)) \
  &= - (4 pi^2)/(k) sum_l sum_m sum_l' sum_m' braket(E_vb(k) l' m', hat(T), E_vb(k) l m) Y_(l' m') (vu(k)) Y_(l m)^* (vu(k)) \
$
由于$hat(T)$算符与$hat(vb(L))^2$和$hat(L)_z$对易，因此
$
  braket(E_vb(k) l' m', hat(T), E_vb(k) l m) = T_l (E_vb(k)) delta_(l' l) delta_(m' m)
$
$T_l (E_vb(k))$称为*在壳的*$l$-分波$T$矩阵元。于是得到
$
  f(vb(k)', vb(k)) & = - (4 pi^2)/(k) sum_l sum_m T_l (E_vb(k)) Y_(l m) (vu(k)) Y_(l m)^* (vu(k)) \
$
取$vb(k)$的方向为$z$方向
$
  Y_(l m) (vu(k)) = sqrt((2 l + 1)/(4 pi)) delta_(m 0) = Y_(l m)^* (vu(k))
$
这样，对$m$求和只有$m = 0$一项。设$vb(k)'$和$vb(k)$的夹角为$theta$，则
$
  Y_(l 0) (vu(k)') = sqrt((2 l + 1)/(4 pi)) P_l (cos theta) \
$
其中$P_l (cos theta)$是Legendre多项式。

定义*分波散射振幅*
$
  f_l (k) = - pi/k T_l (E_vb(k))
$
则*总散射振幅*可以写为
$
  f(vb(k)', vb(k)) = f(theta, k) = sum_(l=0)^oo (2 l + 1) f_l (k) P_l (cos theta)
$
回顾波函数在无穷远处的渐近行为
$
  braket(vb(r), psi^+) -> 1/(2 pi)^(3/2) (e^(i k z) + f(theta, k) e^(i k r)/r) \
$
利用
$
  e^(i k z) & = sum_(l=0)^oo i^l (2 l + 1) j_l (k r) P_l (cos theta) \
  j_l (k r) & -> 1/(k r) sin(k r - (l pi)/2) 当 r -> oo \
$
可将波函数渐进行为写为
$
  braket(vb(r), psi^+) & -> 1/(2 pi)^(3/2) sum_(l=0)^oo (2 l + 1) (P_l (cos theta))/(2 i k) ((1 + 2 i k f_l (k)) (e^(i k r))/r - e^(-i(k r - l pi))/r)
$

#newpara()

根据*概率守恒*，概率流密度满足
$
  pdv(rho, t) + div vb(J) = 0 => integral vb(J) dot dd(vb(S)) = 0
$
即在任意球面上入射和出射的粒子通量应相等。由于角动量守恒，每个分波都应满足这个要求。定义$l$-分波的$S$矩阵元
$
  vb(J) = (hbar)/(2 m i) (psi^* grad psi - psi grad psi^*)
$
$
  S_l (k) = 1 + 2 i k f_l (k) = 1 - 2 pi i T_l (k)
$
概率守恒要求
$
  abs(S_l (k)) = 1
$
此即幺正性关系。也可以仔细计算概率流密度$vb(J)$来证明这个结论。
#note[
  在含时散射中，$hat(S)$是一个幺正算符，与这里的$S$矩阵元对应。
]
所以*在无穷远处出射波只有相位的改变*。定义
$
  S_l (k) = e^(2 i delta_l)
$
$delta_l = delta_l (k)$被称为$l$-分波的*散射相移*
$
  f_l (k) = (e^(2 i delta_l) - 1)/(2 i k) = (e^(i delta_l) sin delta_l)/k = 1/(k cot delta_l - i k)
$
总散射振幅可以用各分波相移表达为
$
  f(theta) & = sum_(l=0)^oo (2 l + 1) (e^(2 i delta_l) - 1)/(2 i k) P_l (cos theta) \
           & = 1/k sum_(l=0)^oo (2 l + 1) e^(i delta_l) sin delta_l P_l (cos theta)
$
总散射截面可以直接计算为
$
  sigma_t &= integral abs(f(theta))^2 dd(Omega)\
  &= 1/k^2 sum_l sum_l' (2l+1)(2l'+1) e^(-i delta_l) sin delta_l e^(i delta_l') sin delta_l' integral_0^(2pi) dd(phi) integral_(-1)^1 dd(cos theta) P_l (cos theta) P_(l') (cos theta) \
  &= (4pi)/k sum_l (2l+1) sin^2 delta_l \
$
可以验证光学定理
$
  sigma_t = (4 pi)/k Im f(0)
$
#newpara()

按照分波展开，散射振幅的计算归结为散射相移的计算。波函数的一般形式可写为
$
  braket(vb(r), psi^+) = 1/(2 pi)^(3/2) sum_(l=0)^oo (2 l + 1) i^' A_l (r) P_l (cos theta) \
$
径向函数$A_l (r)$满足的方程为
$
  1/r^2 dv(, r) (r^2 dv(, r) A_l (r)) + (k^2 - (2m)/hbar^2 V(r) - l(l+1)/r^2) A_l (r) = 0
$
定义$u_l (r) ≡ r A_l (r)$，则$u_l (r)$满足
$
  dv(u_l, r, 2) + (k^2 - (2m)/hbar^2 V(r) - l(l+1)/r^2) u_l (r) = 0
$
按照散射理论的精神，在边界$C$以外，势场为零。假设当$r > R$时，$V = 0$($R$称为力程)，则在$r > R$区域$A_l (r)$的一般解为
$
  A_l (r) = a_l j_l (k r) + b_l n_l (k r)
$
其中$j_l (x)$和$n_l (x)$分别为第一类和第二类球Bessel函数。或者利用球Hankel函数
$
  h_l^((1)) = j_l + i n_l, h_l^((2)) = j_l - i n_l
$
将一般解写为
$
  A_l (r) = c_l^((1)) h_l^((1)) (k r) + c_l^((2)) h_l^((2)) (k r)
$
利用球Bessel函数或者球Hankel函数在$r -> oo$时的渐近行为，可以确定叠加系数与散射相移 $δ_l$的关系。

当$r -> oo$时，球Hankel函数的渐进行为是
$
  h_l^((1)) (k r) & -> e^(i (k r - (l pi)/2))/(k r) , h_l^((2)) (k r) & -> e^(-i (k r - (l pi)/2))/(k r)
$
与渐近行为
$
  braket(vb(r), psi^+) & -> 1/(2 pi)^(3/2) sum_(l=0)^oo (2 l + 1) (P_l (cos theta))/(2 i k) ((1 + 2 i k f_l (k)) (e^(i k r))/r - e^(-i(k r - l pi))/r)
$
比较可定出系数
$
  c_l^((1)) = 1/2 e^(2 i delta_l), c_l^((2)) = 1/2
$
因此，在$r > R$区域，径向函数的解为
$
  A_l (r) & = 1/2 e^(2 i delta_l) h_l^((1)) (k r) + 1/2 h_l^((2)) (k r) \
          & = e^(i delta_l) (j_l (k r) cos delta_l - n_l (k r) sin delta_l) \
          & -> 1/(k r) sin(k r - (l pi)/2 + delta_l), r -> oo
$
显然，散射相移$δ_l$取决于势场$V(r)$和入射能量$E$。因此需要求解相互作用区域即$r < R$区域的径向方程为（$u_l = r A_l$）
$
  dv(u_l, r, 2) + (k^2 - (2m)/hbar^2 V(r) - l(l+1)/r^2) u_l (r) = 0
$
在$r -> 0$时满足边界条件
$
  u_l (r=0) =0
$
另一个边界条件为径向函数$u_l (r)$及其导数在$r = R$处连续。

#example(subname: [硬球散射])[
  $
    V(r) = cases(oo&r<R, 0&r>R)
  $
  这个问题无需计算$β_l$。利用$r = R$处$A_l (r = R) = 0$得到
  $
    j_l (k R) cos delta_l = n_l (k R) sin delta_l => tan delta_l = (j_l (k R))/(n_l (k R))
  $
  对于$s$波散射$(l = 0)$
  $
    tan delta_0 = (j_0 (k R))/(n_0 (k R)) = (sin(k R)/(k R)) / ((-cos(k R))/(k R)) = -tan(k R)
  $
  得到
  $
    delta_0 = - k R
  $
  $r > R$区域径向波函数为(除去相因子$e^(i δ_0)$
  $
    A_(l=0) (r) & prop (sin k r)/(k r) cos delta_0 + (cos k r)/(k r) sin delta_0 \
                & = 1/(k r) sin(k r + delta_0)
  $
]
#figure(
  image("pic/2025-12-14-14-14-47.png", width: 80%),
  numbering: none,
)

#example(subname: [低能散射与束缚态])[
  在径向方程$(u_l ≡ r A_l)$
  $
    dv(u_l, r, 2) + (k^2 - (2m)/hbar^2 V(r) - l(l+1)/r^2) u_l (r) = 0
  $
  中，取$V = 0$得到边界外的解
  $
    dv(v_l, r, 2) + (k^2 - l(l+1)/r^2) v_l (r) = 0
  $
  以上两式分别乘以$v_l$和$u_l$，相减后得到
  $
    v_l dv(u_l, r, 2) - u_l dv(v_l, r, 2) = - (2m)/hbar^2 v_l u_l V(r)\
    => dv(, r) (v_l dv(u_l, r) - u_l dv(v_l, r)) = - (2m)/hbar^2 v_l u_l V(r)
  $
  两边对$r$从$0$到$oo$积分得到
  $
    eval(v_l dv(u_l, r) - u_l dv(v_l, r))_(r->oo) = - (2m)/hbar^2 integral_0^(oo) v_l u_l V(r) dd(r)
  $
  $v_l$的解即为$v_l (r) = r j_l (k r)$。利用$r -> oo$的渐进行为，得到
  $
    (e^(i delta_l) sin delta_l)/k = - (2m)/hbar^2 integral_0^(oo) A_l (r) j_l (k r) V(r) r^2 dd(r)
  $
  在低能极限下，$1/k$远大于力程，我们预期$A_l (r)$与$j_l (k r)$区别不大。因此，右边$tilde k^(2l)$。这意味着$δ_l$也是很小的，因此左边$tilde δ_l/k$。所以，当$k -> 0$时，相移的行为是
  $
    delta_l tilde k^(2l + 1)
  $
  这意味着，对于*短程势的低能散射问题，只有$s$波是重要的*。
]

#newpara()

考虑在如下势垒$(V_0 > 0)$或势阱$(V_0 < 0)$中的$s$波散射。
$
  V(r) = cases(V_0&r<R, 0&r>R)
$
只考虑$l = 0$。在$r > R$区域，波函数为
$
  A(r) & = e^(i delta_0) (j_0 (k r) cos delta_0 - n_0 (k r) sin delta_0) \
       & = e^(i delta_0) (sin(k r + delta_0))/(k r) \
$
这实际上就是径向方程$u=r A$
$
  dv(u, r, 2) + k^2 u = 0
$
的一般解(差一个整体相因子$e^(i delta_0)$)。

在$r < R$区域，径向方程为
$
  dv(u, r, 2) + (2m)/hbar^2 (E - V_0) u = 0, E = (hbar^2 k^2)/(2m)
$
满足边界条件$u(r = 0) = 0$的解为
$
  u(r) = alpha sin(k' r), k' = sqrt((2m(E - V_0))/hbar^2)
$
对于势垒且$E < V_0$的情形，此式也成立，利用
$
  sin(i kappa r) = i sinh(kappa r), k' = i kappa, kappa = sqrt((2m(V_0 - E))/hbar^2)
$
接下来利用$r = R$处的连续性条件。如果只需要定出散射相移$δ_0$，可利用$u'/u$在$r = R$处连续。得到
$
  tan(k R + delta_0)/k = tan(k' R)/k'\
  => 1/k (tan(k R) + tan delta_0)/(1 - tan(k R) tan delta_0) = tan(k' R)/k'
$
解得
$
  tan delta_0 = (k tan(k' R) - k' tan(k R))/(k' + k tan(k R)tan(k' R))
$
对于势垒且$E < V_0$的情形，利用
$
  tan delta_0 = (k tanh(kappa R) - kappa tan(k R))/(kappa + k tan(k R)tanh(kappa R))
$
#figure(
  image("pic/2025-12-14-16-07-11.png", width: 80%),
  numbering: none,
)
- 势阱是吸引势，相移$δ_0/k$体现为左移
- 势垒是排斥势，相移$δ_0/k$体现为右移
考虑*势阱*的情形$(V_0 < 0)$。定义代表势阱深度的参数
$
  xi = sqrt((2 m abs(V_0))/hbar^2) R => k' R = sqrt((k R)^2 + xi^2)
$
则散射相移为$(x ≡ k R)$
$
  tan delta_0 = (x tan(sqrt(x^2 + xi^2)) - sqrt(x^2 + xi^2) tan x)/(sqrt(x^2 + xi^2) + x tan x tan(sqrt(x^2 + xi^2)))
$
由此可计算参数$xi$取不同值时，相移$δ_0$与$x = k R$的关系。我们发现，当势阱中出现*第一个束缚态*，即
$
  xi = pi/2
$
时，*散射相移的行为发生突变*。
#figure(
  image("pic/2025-12-14-16-11-51.png", width: 80%),
  numbering: none,
)
$
  xi = pi/2 => sqrt((2 m abs(V_0))/hbar^2) R = pi/2
$
这意味着势阱深度到一定程度才有束缚态。

一般性地，考虑散射能量$E -> 0$即$k -> 0$的情形。在$r$大于力程的区域，$l = 0$的径向方程变成
$
  dv(u, r, 2)= 0
$
其一般解为
$
  u(r) = c(r - a), A(r) = c(1-a/r)
$
另一方面，对于有限散射能量，力程以外区域的波函数为
$
  A(r) = e^(i delta_0) (sin(k r + delta_0))/(k r) = e^(i delta_0) (sin(k r) cos delta_0 + cos(k r) sin delta_0)/(k r)
$
取$k -> 0$的极限，得到
$
  lim_(k->0) (sin(k r) cos delta_0 + cos(k r) sin delta_0)/(k r) = (lim_(k->0) cos delta_0) (1 + 1/r lim_(k->0) (tan delta_0)/k)
$
注意$delta_0 = delta_0(k) prop k$，因此常数$a$由下式定出
$
  a = - lim_(k->0) (tan delta_0)/k, lim_(k->0) k cot delta_0 = -1/a
$
这个常数$a$被称为*散射长度*。利用光学定理可以得到$k -> 0$时的散射截面为
$
  sigma_t = lim_(k->0) (4 pi)/k Im 1/(k cot delta_0 - i k) = 4 pi a^2
$
#figure(
  image("pic/2025-12-14-16-23-38.png", width: 50%),
  numbering: none,
)
- 对于势垒：$a>0$，散射长度是正的
- 对于势阱：在没有束缚态时，$a<0$，散射长度是负的
- 对于势阱：在出现束缚态时，$a>0$，散射长度是正的
对于*排斥势*，散射长度$a > 0$，大致与力程同量级。例如，对于势垒$(V_0 > 0)$，可计算出散射长度为
$
  a/R = 1 - (tanh xi)/xi, xi = sqrt((2 m V_0)/hbar^2) R
$
对于硬球散射，$a = R$。
#figure(
  image("pic/2025-12-14-16-24-54.png", width: 30%),
  numbering: none,
)
对于*吸引势*，情形则完全不同。散射长度$a$可正可负，而且可以远大于力程。对于势阱$(V_0 < 0)$的散射，可以计算出散射长度为
$
  a/R = 1 - (tan xi)/xi, xi = sqrt((2 m abs(V_0))/hbar^2) R
$
#figure(
  image("pic/2025-12-14-16-30-41.png", width: 30%),
  numbering: none,
)
可以看到，散射长度在
$
  xi = sqrt((2 m abs(V_0))/hbar^2) R = ((2n + 1) pi)/2, n = 1, 2, 3, ...
$
处发散，这些点是*势阱束缚态的散射共振*。考虑势阱中的*束缚态*问题$(E < 0)$。径向波函数的解为
$
  u(r) = cases(
    alpha sin(sqrt((2m(E-V_0))/hbar^2) r) " " & r < R,
    beta exp(- sqrt(-2m E)/hbar^2 r) " " & r > R
  )
$
利用$u'/u$的连续性可得到$E$满足的方程
$
  cot(sqrt((2m(E-V_0))/hbar^2) R) = - sqrt((-E)/(E - V_0))
$
可见，满足上述条件时，势阱刚好要容纳一个新的束缚态$(E -> 0^-)$，同时，散射长度发散，并发生符号改变。这种现象被称为*散射共振*。在$xi$略大于共振点的地方，散射长度$a$很大，且是正的，对应束缚态能级$E -> 0^-$。在$r > R$的区域，束缚态$E -> 0^-$和散射态$E -> 0^+$的波函数分别为
$
  u(r) = cases(
    alpha (r - a) " " & E-> 0^+,
    beta exp(- kappa r) " " & E-> 0^-
  )
$
令两者的$u'/u$相等，得到$(R ≪ a)$
$
  kappa tilde.eq 1/a
$
因此，接近共振点处的束缚态能量为
$
  E = - (hbar^2 kappa^2)/(2m) tilde.eq - (hbar^2)/(2m a^2)
$
#newpara()
上述结果是普遍的，对于一般性的短程势也成立的。我们还可以证明，*散射振幅或者$S$矩阵的奇点正好就是束缚态能量*。

考虑$l = 0$的情形。对于散射态，径向波函数 $A_(l=0) (r)$在$r -> oo$时的渐进行为是
$
  A_(l = 0) (r) tilde S_(l = 0) (k) e^(i k r)/r - e^(-i k r)/r, E = (hbar^2 k^2)/(2m) > 0
$
而对于束缚态，径向波函数的渐进行为是
$
  A_(l = 0) (r) tilde e^(- kappa r)/r , E = - (hbar^2 kappa^2)/(2m) < 0
$
需要注意的是，对于束缚态，$kappa$只能取某些特定的值。我们预期，如果*将$k$延拓到复数域，那么在虚轴上对应的就是束缚态的情形*。将$kappa$用$E < 0$表示为$(epsilon -> 0^+)$
$
  hbar kappa = sqrt(- 2 m(E + i epsilon))
$
若将此式延拓到$E > 0$的情形，则得到
$
  hbar kappa = - i sqrt(2 m E) = - i hbar k => k = - i kappa
$
在物理上，这表示出射波对应束缚态，入射波则无对应。两者系数的比值是$S_(l=0) (k)$。对于束缚态来说，没有入射波对应意味着两个系数的比值是$oo$。所以，如果认为将$k$延拓到复数域，虚轴上对应束缚态，那么$S_(l=0) (k)$在虚轴上将有一些特定的奇点$k = i kappa$，这些$kappa$对应的就是束缚态能量。

对于给定的势$V(r)$，若求得散射振幅
$
  f_(l=0) (k) = 1/(k cot delta_0 - i k)
$
则$S_(l=0) (k)$可计算为
$
  S_(l=0) (k) = 1 + 2 i k f_(l=0) (k) = (cot delta_0 + i)/(cot delta_0 - i)
$
令
$
  k = i kappa = i sqrt(- 2 m E)
$
求出$E < 0$的奇点，即为束缚态能量。

以有限深方势阱为例。前面已求得
$
  cot delta_0 = (k' + k tan(k R) tan(k' R))/(k tan(k R) - k' tan(k' R))
$
令$k = i kappa$，得到
$
  cot delta_0 = (k' - i kappa tanh(kappa R) tan(k' R))/(i kappa tan(k' R) + i k' tanh(kappa R))
$
因此，决定奇点的方程为
$
  cot delta_0 - i = 0 => cot(k' R) = - kappa/k'
$
与直接求解束缚态问题得到的方程一致。

== 对称性考虑

最后，讨论一下空间反演和时间反演对称性对散射振幅的限制。如果体系具有某种对称性，且对称变换算符$hat(U)$是幺正的，则
$
  hat(U) hat(H)_0 hat(U)^dagger = hat(H)_0\
  hat(U) hat(V) hat(U)^dagger = hat(V)
$
根据$T$算符满足的方程
$
  hat(T) = hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(T)
$
可知$T$算符在对称变换下也不变
$
  hat(U) hat(T) hat(U)^dagger = hat(T)
$
定义
$
  ket(tilde(vb(k))) = hat(U) ket(vb(k)), ket(tilde(vb(k)')) = hat(U) ket(vb(k)')
$
则
$
  braket(tilde(vb(k)'), hat(T), tilde(vb(k))) = braket(vb(k)', hat(U)^dagger hat(U) hat(T) hat(U)^dagger hat(U), vb(k)) = braket(vb(k)', hat(T), vb(k))
$
若体系具有宇称反演对称性$(hat(U) = hat(P))$，根据
$
  hat(P) ket(k) = - ket(k)
$
得到
$
  braket(- tilde(vb(k)'), hat(T), - tilde(vb(k))) = braket(vb(k)', hat(T), vb(k))
$
#newpara()

时间反演变换$hat(Theta)$是反幺正变换。根据$T$算符满足的方程
$
  hat(T) = hat(V) + hat(V) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(T)
$
若体系在时间反演变换下不变，即
$
  hat(Theta) hat(H)_0 hat(Theta)^(-1) = hat(H)_0\
  hat(Theta) hat(V) hat(Theta)^(-1) = hat(V)
$
利用时间反演算符的性质得到
$
  hat(Theta) 1/(E_vb(k) - hat(H)_0 + i epsilon) hat(Theta)^(-1) = 1/(E_vb(k) - hat(H)_0 - i epsilon)
$
因此
$
  hat(Theta) hat(T) hat(Theta)^(-1) = hat(T)^dagger
$
对于任意量子态$ket(alpha)$和$ket(beta)$，定义时间反演态
$
  ket(tilde(alpha)) = hat(Theta) ket(alpha), ket(tilde(beta)) = hat(Theta) ket(beta)
$
反幺正变换要求
$
  braket(tilde(alpha), tilde(beta)) = braket(beta, alpha)
$
令
$
  ket(alpha) = hat(T) ket(vb(k)), ket(beta) = ket(vb(k)')
$
得
$
  ket(tilde(alpha)) &= hat(Theta) hat(T) ket(vb(k)) = hat(Theta) hat(T) hat(Theta)^(-1) hat(Theta) ket(vb(k)) = hat(T)^dagger hat(Theta) ket(vb(k)) = hat(T)^dagger ket(- vb(k))\
  ket(tilde(beta)) &= hat(Theta) ket(vb(k)') = ket(- vb(k)')
$
因此，若体系具有时间反演不变性，则有
$
  braket(- vb(k)', hat(T), - vb(k)) = braket(vb(k)', hat(T), vb(k))
$
#newpara()

若体系同时具有空间反演和时间反演不变性，根据以上结果，得到
$
  braket(vb(k)', hat(T), vb(k)) = braket(- vb(k), hat(T), - vb(k)') = braket(vb(k), hat(T), vb(k)')
$
第一个等号利用时间反演不变性，第二个等号利用空间反演不变性。所以
$
  f(vb(k)', vb(k)) = f(vb(k), vb(k)') => dv(sigma, Omega) (vb(k)' -> vb(k)) = dv(sigma, Omega) (vb(k) -> vb(k)')
$
这个结果被称为*细致平衡*。
