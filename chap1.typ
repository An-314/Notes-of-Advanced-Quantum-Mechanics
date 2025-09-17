#import "@preview/scripst:1.1.1": *
#import "@preview/tablex:0.0.9": colspanx, gridx, hlinex, rowspanx, tablex, vlinex
#import "@preview/fletcher:0.5.5" as fletcher: diagram, edge, node

= 量子力学的抽象形式和表象理论

== 引言

物理系统的动力学理论，通常需要包含如下几个要素：
- 系统的*状态*用什么量来描述？
  - 经典力学的Hamilton形式：正则坐标$q$和正则动量$p$
- 系统的*可观测力学量*是什么？
  - 经典力学的Hamilton形式：力学量是$q$和$p$的函数$f(q, p)$
- 系统状态的*时间演化方程*或者说*运动方程*是什么？
  - 经典力学的Hamilton形式：Hamilton正则方程
    $
      dot(q) = pdv(cal(H), p), dot(p) = - pdv(cal(H), q)
    $
- 对系统的力学量进行*测量*，得到的结果是什么？
  - 经典力学的Hamilton形式: 结果就是力学量$f(q, p)$
虽然本科学习的量子力学已经足够烧脑，但还是很具象的：
- 粒子的状态由位形空间的波函数$psi(vb(r), t)$来描述
- 粒子的力学量也是位形空间的函数或微分算符，例如动量算符$hat(p)_x = - hbar dv(, x)$
- 粒子状态的时间演化满足Schrödinger方程
  $
    i hbar pdv(, t) psi(vb(r), t) = (- hbar^2/(2m) laplacian + V(vb(r))) psi(vb(r), t)
  $
  其定态是典型的微分方程的本征值问题

#note[
  这是非相对论的，单粒子的，位形空间（坐标表象）的量子力学，我们后面将要推广到更一般的情形。

  实际上，量子系统存在一个抽象的“原神” 。尤其是，在将上述理论推广到其他更为复杂的量子体系 (量子多体、量子场论等系统)的时候，这样一套抽象的形式是很方便的。
]

#figure(
  three-line-table[
    | 动力学理论 | 经典力学的Hamilton形式 | 本科量子力学 |
    | --- | --- | --- |
    | 状态 | $(q, p)$ | $psi(vb(r), t)$ 波函数 |
    | 力学量 | $f(q, p)$ | $hat(x) = x, hat(p)_x = - i hbar dv(, x), [hat(x), hat(p)_x] = i hbar$ |
    | 时间演化 | Hamilton方程 \ $dot(q) = pdv(cal(H), p), dot(p) = - pdv(cal(H), q)$ | Schrödinger方程 \ $i hbar pdv(, t) psi(vb(r), t) = (- hbar^2/(2m) laplacian + V(vb(r))) psi(vb(r), t)$ |
    | 测量 | $f(q, p)$ | 几率诠释 |
  ],
  kind: table,
  caption: [经典力学和本科量子力学的动力学理论比较],
)

== 量子态

*量子力学的第一个基本假设*是关于如何描述量子系统的状态。在本科量子力学中，是这样表述的：

#definition(subname: [量子力学的第一个基本假设])[
  一个微观粒子的状态总可以用一个*波函数*$ψ(vb(r), t)$来完全描述。波函数是粒子坐标和时间的复值函数，模平方$abs(ψ)^2 ≡ ψ^* ψ$称为概率密度，就是说，在波函数分布区域的小体积元$dd(V)$中找到粒子的概率由
  $
    dd(P) = abs(ψ(vb(r), t))^2 dd(V)
  $
  给出。如果$ψ_1$和$ψ_2$是粒子可能的状态的波函数，那么任意的复系数线性叠加$ψ = c_1 ψ_1 + c_2 ψ_2$也是粒子可能的状态。
]

这种表述，实际上是微观粒子的“原神”在*坐标空间*(位形空间)的具体表现形式。它就是系统的*Hilbert空间*，简单地说，它是定义在复数域上的完备的内积空间。找到了“原神”，我们就可以把第一个基本假设用一种更为抽象的方式表述：

#definition(subname: [量子力学的第一个基本假设（抽象形式）])[
  一个量子系统的状态总可以用Hilbert空间$cal(H)$中的一个*态矢量*$ket(psi)$来完全描述。
]
通常说系统的量子态，就是指系统的*态矢*。Dirac给态矢设计了一套符号体系，称为*Dirac符号*。

态矢用符号$ket(psi)$来表示。若$ket(psi_1)$和$ket(psi_2)$是系统可能的两个量子态，那么任意的复系数线性叠加$ket(psi) = c_1 ket(psi_1) + c_2 ket(psi_2)$也是系统可能的状态。

波函数的概率诠释在这里成了不必要的表述，在后面的测量公设中可以看到，相应的概率密度是$abs(braket(vb(r), psi))^2$。

=== 量子系统的Hilbert空间

考虑一个*复数线性空间*，其中的矢量我们用Dirac符号表示为$ket(psi)$，称为右矢(ket)，这个空间称为*右矢空间*(ket space)。

1. 若$c$为任意复数，则
  $
    c ket(psi) = ket(psi) c
  $
2. 当$c=0$时，上述右矢为零矢量 (null ket)。
  #note[
    注意$ket(0),0$不一致，前者是基态，后者是零矢量。
  ]

引入右矢空间的*对偶*(dual) 空间即*左矢空间*(bra space)，其中存在$bra(psi)$的对偶矢量$bra(psi)$，称为左矢(bra)。这个对偶关系可以写为
$
                         ket(psi) & <-> bra(psi) \
                       c ket(psi) & <-> c^* bra(psi) \
  c_1 ket(psi_1) + c_2 ket(psi_2) & <-> c_1^* bra(psi_1) + c_2^* bra(psi_2)
$

#newpara()

接下来我们定义态矢的*内积*。将左矢$bra(phi)$和右矢$ket(psi)$的内积写为（有时也称为$ket(psi)$和$bra(phi)$的内积或交叠 overlap）
$
  bra(phi) dot ket(psi) ≡ braket(phi, psi)
$
可以形象地理解为$ket(alpha)$向$bra(beta)$的投影。我们有如下假设：
+ 内积一般是一个复数，并且有
  $
    braket(phi, psi) = braket(psi, phi)^*
  $
+ 左矢和与自己对偶的右矢的内积是半正定的
  $
    braket(psi, psi) ≥ 0
  $
  当且仅当$ket(psi)$为零矢量时，等号成立。这是量子力学的概率诠释要求的。
+ 若$ket(psi) = c_1 ket(psi_1) + c_2 ket(psi_2)$，则
  $
    braket(phi, psi) = c_1 braket(phi, psi_1) + c_2 braket(phi, psi_2)
  $

之后引入性质
- 正交性：态矢$phi$和$psi$正交，当且仅当
  $
    braket(phi, psi) = 0 "  " ("implies" braket(psi, phi) = 0)
  $
- 归一性：$sqrt(braket(psi, psi))$称为态矢$ket(psi)$的模。考虑非零态矢$ket(psi)$，如果$braket(psi)=1$，则称态矢$ket(psi)$为归一的。如果$braket(psi)!=1$但为有限值，则其总可以归一化为
  $
    ket(tilde(psi)) = 1/sqrt(braket(psi)) ket(psi)
  $
  注意：量子力学中仍然“非法”使用不能归一化的态矢，例如坐标本征态$ket(x)$和动量本征态$ket(p)$。
