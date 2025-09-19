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

*量子力学的第一个基本假设*是关于如何*描述量子系统的状态*。在本科量子力学中，是这样表述的：

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

引入右矢空间的*对偶*(dual) 空间即*左矢空间*(bra space)，其中存在$bra(psi)$的对偶矢量$bra(psi)$，称为*左矢*(bra)。这个对偶关系可以写为
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
- *正交性*：态矢$ket(phi)$和$ket(psi)$正交，当且仅当
  $
    braket(phi, psi) = 0 "  " ("implies" braket(psi, phi) = 0)
  $
- *归一性*：$sqrt(braket(psi, psi))$称为态矢$ket(psi)$的模。考虑非零态矢$ket(psi)$，如果$braket(psi)=1$，则称态矢$ket(psi)$为归一的。如果$braket(psi)!=1$但为有限值，则其总可以归一化为
  $
    ket(tilde(psi)) = 1/sqrt(braket(psi)) ket(psi)
  $
  注意：量子力学中仍然“非法”使用不能归一化的态矢，例如坐标本征态$ket(x)$和动量本征态$ket(p)$。

  那么对于
  $
    ket(alpha) = sum_n a_n ket(n)\
    ket(beta) = sum_n b_n ket(n) <-> bra(beta) = sum_n b_n^* bra(n)\
    braket(m, n) = delta_(m n)
  $
  就有
  $
    braket(alpha, beta) = sum_n b^*_n a_n
  $

定义了内积以后，复数线性空间便升级为*内积空间*。在此空间中，如果存在一组矢量${ket(n)}$，使得空间中任意矢量$ket(psi)$都可以用这组矢量展开
$
  ket(psi) = sum_n c_n ket(n)
$
则称此空间是*完备*的，这组矢量称为*基矢量*。通常选取正交归一的基矢量，即满足
$
  braket(m, n) = delta_(m n)
$
此时，展开系数$c_n$可以很方便地求出为
$
  c_n = braket(n, psi)
$
因此我们得到
$
  ket(psi) = sum_n braket(n, psi) ket(n) = sum_n ket(n) braket(n, psi)
$
由于$ket(psi)$是任意矢量，所以 ($hat(I)$为单位算符)
$
  sum_n ketbra(n) = hat(I)
$
这称为内积空间的*完备性条件*。满足此条件的空间称为*Hilbert空间*。

不过，量子力学中仍然“非法”使用不可归一的基矢，此时基矢指标是*连续*的，应将求和改为积分
$
  integral dd(xi) ketbra(xi) = hat(I)
$
若基矢既有离散部分也有连续部分，则写为
$
  sum_n ketbra(n) + integral dd(xi) ketbra(xi) = hat(I)
$

#newpara()

在上面论述完备性的过程中，我们不可避免地遇到了类似$ketbra(alpha, beta)$的东西，称为态矢的*外积*，即
$
  ket(beta) dot bra(alpha) ≡ ketbra(beta, alpha)
$
按照 Dirac 的结合律假设(associative axiom)，外积实际上是Hilbert空间中的算符。根据此假设，有如下的结果
$
  (ketbra(beta, alpha)) dot ket(gamma) & = ket(beta) (braket(alpha, gamma)) = braket(alpha, gamma) ket(beta) \
  bra(delta) dot (ketbra(beta, alpha)) & = braket(delta, beta) ket(alpha)
$
由于$braket(alpha, gamma)$和$braket(delta, beta)$都是普通的复数，因此$ketbra(beta, alpha)$作用在任意右矢和左矢的结果都是Hilbert空间中的矢量，即$ketbra(beta, alpha)$是Hilbert空间中的算符。

== 力学量

*量子力学的第二个基本假设*是关于*系统的可观测量*，即*力学量*。

#definition(subname: [量子力学的第二个基本假设])[
  任一可观测力学量$A$可以用相应的线性Hermite算符$hat(A)$来表示。这些算符作用于态的波函数，在这种由力学量$A$到算符$hat(A)$的众多对应规则中，基本的规则是坐标$x$和动量$p$向它们的算符$hat(x)$和$hat(p)$对应。这种对应要求
  $
    hat(x) hat(p) - hat(p) hat(x) = i hbar
  $
]
说得更简单一些，在我们常用的坐标空间中，坐标$x$仍然是$x$，而动量变成了微分算符
$
  hat(p)_x = - i hbar dv(, x)
$
#newpara()
若采用抽象的“原神”，则可表述为：
#definition(subname: [量子力学的第二个基本假设（抽象形式）])[
  *量子系统的可观测力学量*是定义在希尔伯特空间中的*线性Hermite变换*，即作用在态矢上的线性Hermite算符。

  在有经典对应的时候，遵循的基本规则是*正则量子化条件*：若经典系统存在一系列正则坐标${q_n}$和正则动量${p_n}$，则量子化后需要满足正则对易关系
  $
    [hat(q)_m, hat(p)_n] = i hbar delta_(m,n)
  $
]
在量子场论等系统中这种对易关系可以是无穷多个，并且也可能是反对易关系(费米子)。一般来说，经典力学中的力学量$A(q, p)$可以通过替换得到对应的量子力学算符$hat(A)(hat(q), hat(p))$，不过可能会遇到不同的排序产生的问题。
#note[
  但例如$q^2p^2 -> hat(q)^2 hat(p)^2 ->^"Hermite" 1/2(hat(q)^2 hat(p)^2 + hat(p)^2 hat(q)^2)$就需要将其Hermite化。
]

=== Hilbert空间中的算符

希尔伯特空间中的算符是对其中矢量的一种变换。考虑定义在某右矢空间中的算符$hat(X)$，它作用在任意右矢$ket(alpha)$上的结果是将其变换为空间中的另一个右矢$ket(alpha')$，即
$
  hat(X) dot ket(alpha) ≡ hat(X) ket(alpha) = ket(alpha')
$
如果两个算符$hat(X)$和$hat(Y)$作用在空间中任意矢量$ket(alpha)$上，结果相同，即
$
  hat(X) ket(alpha) = hat(Y) ket(alpha)
$
则称$hat(X)$和$hat(Y)$是*相等的*，记
$
  hat(X) = hat(Y)
$
#newpara()
如果算符$hat(X)$作用在空间中任意矢量$ket(alpha)$的结果都等于$ket(alpha)$
$
  hat(X) ket(alpha) = ket(alpha)
$
则称$hat(X)$为*单位算符*(identity operator)，记为$hat(I)$。

如果算符$hat(X)$作用在空间中任意矢量$ket(alpha)$的结果都等于零
$
  hat(X) ket(alpha) = ket(0)
$
则称$hat(X)$为*零算符*。

实际上，若$ket(alpha)$是空间的一组完备基矢，上述结论就成立。

*算符的加法*定义为
$
  (hat(X) + hat(Y)) ket(alpha) = hat(X) ket(alpha) + hat(Y) ket(alpha)
$
容易验证，加法满足交换律和结合律
$
  hat(X) + hat(Y) = hat(Y) + hat(X)\
  (hat(X) + hat(Y)) + hat(Z) = hat(X) + (hat(Y) + hat(Z))
$
#newpara()
*线性算符*：在量子力学中遇到的算符绝大部分是线性算符。若$hat(X)$是线性算符，则对空间中的任意态矢$ket(α), ket(β)$和任意复数$a, b, c$有
$
  hat(X) (a ket(alpha) + b ket(beta)) = a hat(X) ket(alpha) + b hat(X) ket(beta) \
  hat(X) (c ket(alpha)) = c hat(X) ket(alpha)
$
#newpara()
算符$hat(X)$也可以作用在左矢空间。它作用在任意左矢$bra(alpha)$上的结果是将其变换为空间中的另一个左矢$bra(alpha'')$，即
$
  bra(alpha) dot hat(X) ≡ bra(alpha) hat(X) = bra(alpha'')
$
需要注意的是右矢 $hat(X) ket(alpha)$ 和左矢$bra(alpha) hat(X)$ 一般并不是互相对偶的。

定义算符$hat(X)^dagger$ ，使得右矢$hat(X) ket(alpha)$和左矢 $bra(alpha) hat(X)^dagger$ 互为对偶，即
$
  hat(X) ket(alpha) & <-> bra(alpha) hat(X)^dagger \
$
算符$hat(X)^dagger$称为算符$hat(X)$的*Hermite共轭*(Hermitian adjoint)或*伴随*。如果
$
  hat(X) = hat(X)^dagger
$
则称算符$hat(X)$为*Hermite算符*(Hermitian operator)。

*算符的乘法*定义为
$
  (hat(X) hat(Y)) ket(alpha) = hat(X) (hat(Y) ket(alpha))
$
容易验证，乘法满足*结合律*
$
  (hat(X) hat(Y)) hat(Z) = hat(X) (hat(Y) hat(Z)) = hat(X) hat(Y) hat(Z)
$
*非对易性*：一般来说，算符的乘法不满足交换律，即不可交换，
$
  hat(X) hat(Y) != hat(Y) hat(X)
$
*算符乘积的Hermite共轭*：算符$hat(X) hat(Y)$的Hermite共轭为
$
  (hat(X) hat(Y))^dagger = hat(Y)^dagger hat(X)^dagger
$

=== Dirac 结合律假设

定义了左矢、右矢和算符以后，它们之间就存在“合法乘积”和“非法乘积”。典型的“非法乘积”列举如下：
$
  ket(alpha) hat(X), hat(X) bra(alpha), ket(alpha) ket(beta), bra(alpha) bra(beta)
$
对于“合法乘积”，Dirac 假设它们都满足结合律。

第一个应用就是前面提到的*外积*
$
  (ketbra(beta, alpha)) ket(gamma) = ket(beta) (braket(alpha, gamma)) = braket(alpha, gamma) ket(beta)
$
第二个应用是力学量的*矩阵元*。根据结合律假设，我们有
$
  bra(beta) (hat(X) ket(alpha)) = (bra(beta) hat(X)) ket(alpha)
$
由于以上的等式，我们将其写成更为紧凑的形式
$
  braket(beta, hat(X), alpha)
$
结果是一个复数，称为力学量在态矢$ket(beta)$和$ket(alpha)$之间的*矩阵元*。

也就是说，在考虑矩阵元$braket(beta, hat(X), alpha)$的时候，既可以理解成 $hat(X)$先向右作用，得到右矢$hat(X) ket(alpha)$，再与左矢$bra(beta)$做内积，也可以理解成$hat(X)$先向左作用，得到左矢$bra(beta) hat(X)$，再与右矢$ket(alpha)$做内积。

#note[
  需要注意的是，这样的理解只有当*$hat(X)$是线性算符*的时候才成立。当遇到*反线性算符*的时候，这样的理解会“翻车”。最著名的例子就是*时间反演算符*$hat(T)$。它是一个典型的反线性算符，满足
  $
    hat(T) (a ket(alpha) + b ket(beta)) = a^* hat(T) ket(alpha) + b^* hat(T) ket(beta)
  $
]
#newpara()

进一步考虑线性算符$hat(X)$的矩阵元$braket(beta, hat(X), alpha)$。由于$bra(alpha) hat(X)^dagger$是$hat(X) ket(alpha)$的对偶，所以
$
  braket(beta, hat(X), alpha) & = bra(beta) dot (hat(X) ket(alpha)) \
                              & = ((bra(alpha) hat(X)^dagger) dot ket(beta))^* \
                              & = (bra(alpha) dot (hat(X)^dagger ket(beta)))^* \
$
即
$
  braket(beta, hat(X), alpha) = braket(alpha, hat(X)^dagger, beta)^*
$
其中，用到了内积的基本假设$braket(alpha, beta) = braket(beta, alpha)^*$。如果$hat(X)$是Hermite算符，则
$
  braket(beta, hat(X), alpha) = braket(alpha, hat(X), beta)^*
$

#newpara()

考虑某个力学量$hat(Omega)$，根据基本假设，它必定为Hermite算符。对于任意态矢$ket(psi), ket(phi)$，有
$
  braket(phi, hat(Omega), psi) = braket(psi, hat(Omega), phi)^*
$
特别是，当两个态矢相同时，我们得到
$
  expval(hat(Omega))_psi = braket(psi, hat(Omega), psi) = braket(psi, hat(Omega), psi)^*
$
这说明在任意态中，*Hermite算符的平均值都为实数*。

考虑Hermite算符$hat(Omega)$的*本征方程*
$
  hat(Omega) ket(omega) = omega_n ket(omega)
$
即$hat(Omega)$的*本征值*为${omega_n}$，*本征态*为${ket(omega_n)}$。

考虑矩阵元$braket(m, hat(Omega), n)$，其中$ket(m)$和$ket(n)$分别是本征值$omega_m$和$omega_n$对应的本征矢。根据 Dirac
结合律假设，可以理解为$hat(Omega)$先向右作用，得到
$
  braket(m, hat(Omega), n) = omega_n braket(m, n)
$
也可以理解为$hat(Omega)$先向左作用，得到
$
  braket(m, hat(Omega), n) = omega_m^* braket(m, n)
$
两种结果应相等，所以有
$
  (omega_n - omega_m^*) braket(m, n) = 0
$
- 当$m = n$时，由于$braket(n, n) != 0$，所以$omega_n = omega_n^*$，即*Hermite算符的本征值为实数*。
- 当$m != n$时，分两种情况讨论。
  - 若 $omega_m != omega_n$，即本征态$ket(m)$和$ket(n)$是非简并的，则$braket(m, n) = 0$，即*Hermite算符属于不同本征值的本征态正交*。
  - 若 $omega_m = omega_n$，即本征态$ket(m)$和$ket(n)$是简并的，那么总可以通过将这两个态矢线性叠加，构造出两个线性无关的态矢$ket(m')$和$ket(n')$，使得它们正交。同理，若有多个态是简并的，那么总可以通过将这这些态矢线性叠加，构造出相同数目的互相正交的态矢。更为物理的做法是引入力*学量完备集*。

对于本征值离散的力学量，$braket(n)$是有限的，其本征态总可以归一化。对于本征值连续的力学量，通常可以将其本征态归一化到$δ$函数。因此，对于力学量，总可以找到一组*正交归一*的本征态。进一步地，我们通常假设这组本征态是*完备*的。因此，*力学量或者力学量完备集就提供了Hilbert空间的一组基矢量*。

== 测量

*量子力学的第三个基本假设*是关于*测量*。

#definition(subname: [量子力学的第三个基本假设])[
  一个微观粒子处于波函数$ψ(vb(r))$的状态，若对它测量可观测力学量$hat(A)$的数值，则所得的$hat(A)$的期望值为
  $
    expval(hat(A))_psi = (integral ψ^*(vb(r)) hat(A) ψ(vb(r)) dd(vb(r)) )/ (integral psi^*(vb(r)) psi(vb(r))dd(vb(r)))
  $
  若$psi(vb(r))$是归一的，则为
  $
    expval(hat(A))_psi = integral ψ^*(vb(r)) hat(A) ψ(vb(r)) dd(vb(r))
  $
]
由于波函数公设中已经给出了Born几率诠释，所以不再说明单次测量的结果，而只给出对量子系综多次测量的平均结果(期望值)。另外，这个表述也只针对单粒子系统，采用坐标空间进行表述。

如果采用抽象的 Dirac 符号，则测量公设可以表述为：

#definition(subname: [量子力学的第三个基本假设（抽象形式）])[
  若一个量子系统处于状态$ket(psi)$，在此状态下测量某力学量$hat(A)$，得到的结果只能是$hat(A)$的某个本征值。也就是说，如果$hat(A) ket(n) = a_n ket(n)$，可将状态$ket(psi)$展开为
  $
    ket(psi) = sum_n c_n ket(n) => c_n = braket(n, psi)
  $
  若$ket(psi)$是归一化的，则测量力学量$hat(A)$的结果为$a_n$的概率为
  $
    P(a_n) = abs(c_n)^2 = abs(braket(n, psi))^2
  $
  #newpara()
  这里假定了本征态无简并，对于有简并的情形，只需简单推广。例
  如：假设$ket(m), ket(n),...$是简并的，即$a_m = a_n = ...$，则测量结果为$a_n$的概率为$abs(c_m)^2 + abs(c_n)^2 + ...$。
]

#newpara()

根据上述假设，可计算出$hat(A)$的*期望值*为
$
  expval(hat(A))_psi = sum_n a_n P(a_n) = sum_n a_n abs(braket(n, psi))^2
$
容易证明
$
  braket(psi, hat(A), psi) = sum_(n,m) c_m^* c_n braket(m, hat(A), n) = sum_(n,m) c_m^* c_n a_n braket(m, n) = sum_n a_n abs(c_n)^2 = expval(hat(A))_psi
$
因此期望值$expval(hat(A))_psi$实际上是*量子力学平均值*(与混合态对比)。

需要注意的是，测量公设中既有随机性，也有确定性。在单次测量中，结果是哪个本征值完全是*随机*的。但是，测量得到某个本征值的*概率是确定的*。当$ket(psi)$是$hat(A)$的*本征态*时，测量结果是完全确定的，称$hat(A)$具有*确定的取值*。

=== 多个力学量同时具有确定值的条件

#theorem[
  若力学量$hat(A)$和$hat(B)$有完备的共同本征态，则它们一定对易。
]

#proof[
  设完备的共同本征态为${ket(n)}$，即
  $
    hat(A) ket(n) = a_n ket(n), hat(B) ket(n) = b_n ket(n)
  $
  由于共同本征态是完备的，所以任意态矢$ket(psi)$都可以展开为
  $
    ket(psi) = sum_n ket(n) braket(n, psi)
  $
  利用此展开式
  $
    [hat(A), hat(B)] ket(psi) & = (hat(A) hat(B) - hat(B) hat(A)) ket(psi) \
                              & = sum_n (hat(A) hat(B) - hat(B) hat(A)) ket(n) braket(n, psi) \
                              & = sum_n (a_n b_n - b_n a_n) ket(n) braket(n, psi) = 0
  $
  由于$ket(psi)$是任意态矢，所以
  $
    [hat(A), hat(B)] = 0
  $
]

#theorem[
  若力学量$hat(A)$和$hat(B)$对易，则它们有完备的共同本征态。
]

#proof[
  设$hat(A)$的完备本征态为${ket(n)}$，即
  $
    hat(A) ket(n) = a_n ket(n)
  $
  由$[hat(A), hat(B)] = 0$，有
  $
    hat(A) hat(B) ket(n) = hat(B) hat(A) ket(n) = a_n hat(B) ket(n)
  $
  这说明$hat(B) ket(n)$也是$hat(A)$对应本征值为$a_n$的本征态。若$a_n$无简并，则$hat(B) ket(n)$与$ket(n)$是同一个态，只能相差一个常数，即
  $
    hat(B) ket(n) = b_n ket(n)
  $
  由于$hat(B)$是力学量，$b_n$必为实数。

  所以$ket(n)$也是$hat(B)$的本征态，即$hat(A), hat(B)$有完备的共同本征态。

  若$a_n$有简并，则可以用Schmidt方法来证明有相同的结果。
]

如果两个力学量对易，它们可以有共同的完备本征态，这两个力学量称为*相容的*(compatible)。当系统处于共同本征态时，它们同时具有确定的取值，即：*存在这样的状态，测量之前两个力学量在客观上都各自具有确定的数值*。因此，有时也称两个力学量“*可以同时测量*”。反之，如果两个力学量不对易，则不存在这样的状态，测量之前两个力学量在客观上都各自具有确定的数值，或者“*不能同时测量*”。

对于多个力学量，同样可以证明：如果*两两对易*，它们可以有共同的完备本征态。当体系处于共同本征态时，它们同时有确定的取值。

#note[
  如果$[hat(A), hat(B)]$且$[hat(B), hat(C)] = 0$，不一定有$[hat(A), hat(C)] = 0$。

  例如角动量$hat(vb(J))$，有
  $
    [hat(vb(J))^2,J_x] = [hat(vb(J))^2,J_y] = [hat(vb(J))^2,J_z] = 0 \
    [J_x,J_y] = i hbar J_z, [J_y,J_z] = i hbar J_x, [J_z,J_x] = i hbar J_y
  $
  因此，$hat(vb(J))^2$与$J_x, J_y, J_z$两两对易，但$J_x, J_y, J_z$彼此不对易。$hat(vb(J))^2, J_x$可以同时有确定的取值，但$J_x, J_y$一般不能同时有确定值。
]

=== 力学量完备集

若力学量$hat(A)$的本征态有*简并*，则本征方程可写为
$
  hat(A) ket(n\,i) = a_n ket(n\,i), i = 1, 2, ...
$
在数学上，可以用Schmidt方法来构造简并子空间中正交归一的本征态。物理的做法是引入力学量完备集来实现正交归一，其它力学量的本征值或*量子数*自动区分这些简并的本征态。(请注意：是区分简并态，不是解除简并)

力学量$hat(A)$的本征态有简并，说明其对应的量子数不足以完全标记这些本征态，必须引入新的力学量即新的量子数来区分这些简并的本征态。而这通常又与*对称性*有关。

引入另一个与$hat(A)$对易的力学量$hat(B)$，$[hat(A), hat(B)] = 0$，则$hat(A), hat(B)$具有完备的共同本征态。如果$hat(A)$的简并本征态对于$hat(B)$来说是不简并的，即
$
  hat(A) ket(n\,i) = a_n ket(n\,i), hat(B) ket(n\,i) = b_i ket(n\,i)
$
其中，对于不同的$i$，$b_i$是不同的。则对于力学量集合${hat(A), hat(B)}$，本征值${a_n , b_i}$仅对应唯一的本征态$ket(n\, i)$，那么本征值集合${a_n , b_i}$或者量子数集合${n, i}$标记的共同本征态自动是正交归一的
$
  braket(m\,j, n\,i) = delta_(m n) delta_(j i)
$
我们称${hat(A), hat(B)}$构成*力学量完备集*。

如果力学量$hat(B)$的本征值仍不足以区分$hat(A)$的简并本征态，则可进一步引入与$hat(A)$和$hat(B)$都对易的力学量$hat(C)$，直到简并态被完全区分，这些两两对易的力学量就是力学量完备集。


#example(subname: [粒子轨道角动量])[
  粒子轨道角动量$hat(vb(L))$的平方满足本征方程
  $
    hat(vb(L))^2 ket(l\,m) = hbar^2 l(l+1) ket(l\,m), m = -l, -l+1, ..., l
  $
  简并度为$2l + 1$。但是这$2l + 1$个正交归一本征态的构造方法不是唯一的。

  引入新的力学量$hat(L)_z$，它与$hat(vb(L))^2$对易，有共同本征态
  $
    hat(L)_z ket(l\,m) = hbar m ket(l\,m)\
    hat(vb(L))^2 ket(l\,m) = hbar^2 l(l+1) ket(l\,m)
  $
  这说明，$hat(L)_z$的量子数$m$完全区分了这$2l + 1$个简并态，${hat(vb(L))^2, hat(L)_z}$构成力学量完备集。
]

=== 不确定度关系

=== 波包塌缩假设

前面实际上是Born几率诠释的抽象版本及其演绎出来的结果，但是并未回答：*测量结束以后，系统的状态是什么*。如果要追问这个问题，则需要附加波包塌缩假设：
#definition(subname: [波包塌缩假设])[

  测量得到结果后，状态$ket(psi)$受到*严重干扰或者破坏*，系统的状*态突*变为测量得到的那个本征值对应的本征态。即若测量结果为$a_n$，则测量结束以后系统状态突变为对应的本征态$ket(n)$。这个过程称为*波包塌缩*。除非$ket(psi)$本身就是力学量$hat(A)$的某个本征态，否则在单次测量后系统塌缩到哪个本征态，就和测量得到哪个本征值一样，是完全随机的、不可事先预测的。
]

在这样的学说(Copenhagen诠释)中，测量及其伴随的塌缩不同于动力学演化(Schrödinger方程)，它是一种*非幺正演化*，这种演化是随机的、不可逆的、破坏相干性的。需要注意的是，波包塌缩假设对量子动力学演绎出的结果(“计算量子力学”)不产生任何影响，因此才有那句名言：*Shut up and calculate.*
