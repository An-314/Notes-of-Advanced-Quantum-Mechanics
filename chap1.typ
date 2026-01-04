#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

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
  *量子系统的可观测力学量*是定义在Hilbert空间中的*线性Hermite变换*，即作用在态矢上的线性Hermite算符。

  在有经典对应的时候，遵循的基本规则是*正则量子化条件*：若经典系统存在一系列正则坐标${q_n}$和正则动量${p_n}$，则量子化后需要满足正则对易关系
  $
    [hat(q)_m, hat(p)_n] = i hbar delta_(m,n)
  $
]
在量子场论等系统中这种对易关系可以是无穷多个，并且也可能是反对易关系(Fermi子)。一般来说，经典力学中的力学量$A(q, p)$可以通过替换得到对应的量子力学算符$hat(A)(hat(q), hat(p))$，不过可能会遇到不同的排序产生的问题。
#note[
  但例如$q^2p^2 -> hat(q)^2 hat(p)^2 ->^"Hermite" 1/2(hat(q)^2 hat(p)^2 + hat(p)^2 hat(q)^2)$就需要将其Hermite化。
]

=== Hilbert空间中的算符

Hilbert空间中的算符是对其中矢量的一种变换。考虑定义在某右矢空间中的算符$hat(X)$，它作用在任意右矢$ket(alpha)$上的结果是将其变换为空间中的另一个右矢$ket(alpha')$，即
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
#proposition(subname: [Hermite算符的性质])[
  对于任意Hermite算符$hat(Omega)$，有如下性质：
  1. *Hermite算符的平均值为实数*。
  2. *Hermite算符的本征值为实数*。
  3. *Hermite算符属于不同本征值的本征态正交*。
]
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
  - 若 $omega_m = omega_n$，即本征态$ket(m)$和$ket(n)$是简并的，那么总可以通过将这两个态矢线性叠加，构造出两个线性无关的态矢$ket(m')$和$ket(n')$，使得它们正交。同理，若有多个态是简并的，那么总可以通过将这这些态矢线性叠加，构造出相同数目的互相正交的态矢。更为物理的做法是引入*力学量完备集*。

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
  这里假定了本征态无简并，对于有简并的情形，只需简单推广。例如：假设$ket(m), ket(n),...$是简并的，即$a_m = a_n = ...$，则测量结果为$a_n$的概率为$abs(c_m)^2 + abs(c_n)^2 + ...$。
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

若$[hat(A), hat(B)] != 0$，则$hat(A)$和$hat(B)$*不能同时具有确定的取值*。

考虑任意状态$ket(psi)$，可以用$hat(A)$的本征态$ket(n)$展开为
$
  ket(psi) = sum_n ket(n) braket(n, psi)
$
力学量$hat(A)$取值为$a_n$的概率为$P(n) = abs(braket(n, psi))^2$。定义力学量的涨落
$
  Delta hat(A) = hat(A) - expval(hat(A))
$
*力学量的不确定度*的平方定义为方差
$
  expval((Delta hat(A))^2) = expval(hat(A)^2) - expval(hat(A))^2
$
可以证明如下关系成立：
#theorem(subname: [不确定度关系])[
  对于任意两个力学量$hat(A)$和$hat(B)$，有不确定度关系
  $
    expval((Delta hat(A))^2) expval((Delta hat(B))^2) >= 1/4 abs(expval([hat(A), hat(B)]))^2
  $
]

#proof[
  对于任意态$ket(alpha), ket(beta)$和任意复数$lambda$，线性组合态$ket(alpha) + lambda ket(beta)$的内积恒为非负，即
  $
    (bra(alpha) + lambda^* bra(beta)) dot (ket(alpha) + lambda ket(beta)) >= 0
  $
  取$lambda = - (braket(beta, alpha))/(braket(beta, beta))$，则有
  $
    braket(alpha, alpha) braket(beta, beta) >= abs(braket(beta, alpha))^2
  $
  此即*Schwarz不等式*。令$ket(alpha) = Delta hat(A) ket(psi)$, $ket(beta) = Delta hat(B) ket(psi)$，上式变为
  $
    expval((Delta hat(A))^2) expval((Delta hat(B))^2) >= abs(braket(psi, Delta hat(A) Delta hat(B), psi))^2
  $
  可将右边写为
  $
    Delta hat(A) Delta hat(B) = 1/2 [Delta hat(A), Delta hat(B)] + 1/2 {Delta hat(A), Delta hat(B)}
  $
  易证，对易子$[∆ hat(A), ∆ hat(B)]$是反Hermite算符，平均值为纯虚数，反对易子${∆ hat(A), ∆ hat(B)}$是Hermite算符，平均值为实数。因此
  $
    abs(expval(Delta hat(A) Delta hat(B)))^2 = (1/4) abs(expval([Delta hat(A), Delta hat(B)]))^2 + (1/4) abs(expval({Delta hat(A), Delta hat(B)}))^2
  $
  利用
  $
    [Delta hat(A), Delta hat(B)] = [hat(A) - expval(hat(A)), hat(B) - expval(hat(B))] = [hat(A), hat(B)]
  $
  得到
  $
    expval((Delta hat(A))^2) expval((Delta hat(B))^2) >= 1/4 abs(expval([hat(A), hat(B)]))^2
  $
]
要使等号成立，即满足*最小不确定度*，量子态$psi$必须满足条件：
- Schwarz 不等式取等号
- $expval({Delta hat(A), Delta hat(B)}) = 0$
两个条件分别要求：
$
  Delta hat(B) ket(psi) = c Delta hat(A) ket(psi)\
  braket(psi, Delta hat(A) Delta hat(B), psi) + braket(psi, Delta hat(B) Delta hat(A), psi) = 0
$
将第一个式子代入第二个式子得到
$
  c braket(psi, (Delta hat(A))^2, psi) + c^* braket(psi, (Delta hat(A))^2, psi) = 0
$
所以$c$必为纯虚数，即$c = i a$，$a$为实数。故最小不确定度对应的量子态$ket(psi)$满足
$
  (hat(B) - expval(hat(B))) ket(psi) = i a (hat(A) - expval(hat(A))) ket(psi)
$

#example(subname: [位置-动量不确定度关系])[
  取$hat(A) = hat(x)$，$hat(B) = hat(p)$，则得到位置-动量不确定度关系
  $
    Delta x Delta p >= hbar / 2
  $
  最小不确定度对应的量子态满足方程
  $
    (hat(p) - expval(hat(p))) ket(psi) = i a (hat(x) - expval(hat(x))) ket(psi)
  $
  在坐标空间(表象)，写为$(expval(p) = p_0, expval(x) = x_0)$
  $
    (-i hbar dv(, x) - p_0) psi(x) = i a (x - x_0) psi(x)
  $
  其解为
  $
    psi(x) = cal(N) e^(- a (x - x_0)^2 / (2 hbar)) e^((i p_0 x) / hbar)
  $
]

=== 波包塌缩假设

前面实际上是Born几率诠释的抽象版本及其演绎出来的结果，但是并未回答：*测量结束以后，系统的状态是什么*。如果要追问这个问题，则需要附加波包塌缩假设：
#definition(subname: [波包塌缩假设])[

  测量得到结果后，状态$ket(psi)$受到*严重干扰或者破坏*，系统的状*态突*变为测量得到的那个本征值对应的本征态。即若测量结果为$a_n$，则测量结束以后系统状态突变为对应的本征态$ket(n)$。这个过程称为*波包塌缩*。除非$ket(psi)$本身就是力学量$hat(A)$的某个本征态，否则在单次测量后系统塌缩到哪个本征态，就和测量得到哪个本征值一样，是完全随机的、不可事先预测的。
]

在这样的学说(Copenhagen诠释)中，测量及其伴随的塌缩不同于动力学演化(Schrödinger方程)，它是一种*非幺正演化*，这种演化是随机的、不可逆的、破坏相干性的。需要注意的是，波包塌缩假设对量子动力学演绎出的结果(“计算量子力学”)不产生任何影响，因此才有那句名言：*Shut up and calculate.*

=== 量子系综

由于波包*塌缩*，我们不可能对同一个量子系统进行反复测量，获得力学量的期望值。为了理解测量得到某个本征值的概率和力学量的期望值，可以引入量子系综的概念。所谓*量子系综*，就是假设存在*大量相同的系统，即都处于相同的状态*$ket(psi)$。

对系综中的任何一个系统测量的时候，得到的结果都是随机的。对系综中的所有系统都进行测量后，将测量结果进行统计，那么就会发现，测量得到某个本征值$a_n$的概率就等于$P(a_n) = abs(braket(n, psi))^2 = abs(c_n)^2$。力学量$hat(A)$的平均值就是按照概率分布$P(a_n)$加权的本征值的期望值
$
  expval(hat(A))_psi = sum_n a_n abs(c_n)^2
$

#example(subname: [量子测量])[
  对量子系统的三个力学量$hat(A), hat(B), hat(C)$分别进行如下两种测量。

  *第一种：*先测$hat(A)$，再测$hat(B)$，最后测$hat(C)$。

  设体系初态为$ket(psi)$，对$hat(A)$进行测量，将以$abs(braket(a, psi))^2$ 的概率得到$hat(A)$的某个本征值$a$（玻恩几率诠释），测量后体系塌缩为$ket(a)$（波包塌缩）。

  然后体系处于状态$ket(a)$，对$hat(B)$进行测量，将以$abs(braket(b, a))^2$的概率得到$hat(B)$的某个本征值$b$，测量后体系塌缩为$ket(b)$。

  之后体系处于状态$ket(b)$，对$hat(C)$进行测量，将以$abs(braket(c, b))^2$的概率得到$hat(C)$的某个本征值$c$，测量后体系塌缩为$ket(c)$。

  因此，三次测量后，得到结果$c$的概率是
  $
    P(c) = sum_(a,b) abs(braket(c, b))^2 abs(braket(b, a))^2 abs(braket(a, psi))^2
  $
  由于前两次测量得到各种$a, b$的可能性都存在，因此需要对$a, b$求和。

  *第二种：*直接测$hat(A)$以后直接测$hat(C)$。

  解答：设体系初态为$ket(psi)$，对$hat(A)$进行测量，将以$abs(braket(a, psi))^2$ 的概率得到$hat(A)$的某个本征值$a$（玻恩几率诠释），测量后体系塌缩为$ket(a)$（波包塌缩）。然后体系处于状态$ket(a)$，对$hat(C)$进行测量，将以$abs(braket(c, a))^2$的概率得到$hat(C)$的某个本征值$c$，测量后体系塌缩为$ket(c)$。因此，两次测量后，得到结果$c$的概率是
  $
    P'(c) = sum_a abs(braket(c, a))^2 abs(braket(a, psi))^2
  $
  由于前一次测量得到各种$a$的可能性都存在，因此需要对$a$求和。

  显然，这两种测量方案得到结果$c$的概率是不同的，即
  $
    P(c) != P'(c)
  $
  但是，如果$[hat(A), hat(B)] = 0$，那么$hat(A), hat(B)$具有完备的共同本征态
  $
    braket(a, b) = delta_(a b)
  $
  因此
  $
    P(c) & = sum_(a,b) abs(braket(c, b))^2 abs(braket(b, a))^2 abs(braket(a, psi))^2 \
         & = sum_a abs(braket(c, a))^2 abs(braket(a, psi))^2 = P'(c)
  $
  即，当$hat(A), hat(B)$对易时，两种测量方案得到结果$c$的概率相同。

  同样可以证明，如果$[hat(B), hat(C)] = 0$，两种测量方案结果也相同。
]

== 动力学演化

*量子力学的第四个基本假设*是关于量子态的*动力学演化*。

#definition(subname: [量子力学的第四个基本假设])[
  一个微观粒子系统的状态波函数满足如下的*Schrödinger方程*
  $
    i hbar pdv(, t) psi(vb(r), t) = hat(H)(vb(r), hat(vb(p))) psi(vb(r), t) = hat(H)(vb(r), - i hbar grad) psi(vb(r), t)
  $
  这里$hat(H)$是系统的Hamilton算符，称为*Hamilton量*
  $
    hat(H) = hat(T) + hat(V)(vb(r)) = hat(vb(p))^2/2m + V(vb(r)) = - hbar^2/(2m) laplacian + V(vb(r))
  $
]
这样表述显然有很大的局限性。
- 第一，只适用于非相对论性单粒子系统，如果要考虑量子多体和量子场论等系统，则需要抽象的形式；
- 第二，只适用于坐标空间 (位形空间)，不能体现量子系统的“原神”。

#definition(subname: [量子力学的第四个基本假设])[
  采用抽象的 Dirac 符号，将系统在 t 时刻的量子态写为 |ψ(t)⟩，则量子态的动力学演化遵守运动方程
  $
    i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
  $
  这个方程仍然称为*Schrödinger方程*，其中$hat(H)$是系统的Hamilton算符或Hamilton量。
]

Schrödinger方程是*关于时间的一阶微分方程*，这是因为我们假设量子态$ket(psi(t))$完全描述了系统的状态。

Schrödinger方程是*线性*的，容易验证：若$ket(psi_1 (t))$和$ket(psi_2 (t))$都是方程的解，则复系数线性叠加 $c_1 ket(psi_1 (t)) + c_2 ket(psi_2 (t))$ 也是方程的解。只要Hamilton量是Hermite算符，那么量子态的模方$braket(psi(t))$将不随时间变化 (概率守恒)。

需要注意的是，抽象的Schrödinger方程理论意义上更高。是一切量子理论都要遵循的时间演化方程，包括量子多体系统、量子场论等。不过，这些系统的自由度、Hilbert空间、相互作用和Hamilton量将会非常复杂。量子态的模方$braket(psi(t))$也不再解释为粒子出现的概率。

大部分对于量子系统的研究，可以归结为两件事情：*造Hamilton量和解Hamilton量*，即造模型和解模型。

如果Hamilton量$hat(H)$已知，且不显含时间，则Schrödinger方程的求解完全归结为求解*定态方程*
$
  hat(H) ket(psi) = E ket(psi)
$
也就是求解Hamilton量的本征值 (能谱) 和本征态。

材料电子结构的计算和格点 QCD 计算等暴力运算，实际上就是在干这件事。理论上，如果Hamilton量的本征值能计算出来，那么统计物理问题也解决了$cal(Z) = sum_n e^(- beta E_n)$。

如果系统的Hamilton量未知，则需要先构造出系统的Hamilton量。根据对应原理，如果系统有*经典对应*，则可以利用经典系统的Hamilton量得到量子的Hamilton量算符。如果不存在经典对应，则只能靠“猜”了。归根结底，量子系统的Hamilton量是什么形式，应该由实验来检验。即：理论上先“猜”出量子Hamilton量，然后将各种计算结果与实验结果对比，检验Hamilton量是否正确。

#example(subname: [构造Hamilton量])[
  如果经典Hamilton量的形式为
  $
    H(q,p) = T(p) + V(q)
  $
  即$q$和$p$是*分离*的(分别出现的不同的项中)，则量子Hamilton量可以直接通过替换$q -> hat(q), p -> hat(p)$得到
  $
    hat(H) = T(hat(p)) + V(hat(q))
  $
]

== 表象理论

虽然量子理论可以写成抽象的形式，但是针对具体的系统，通常我们需要进入具体的*表象*进行计算。

如何进入一个具体的*表象*？条件是*Hilbert空间一组基矢的完备性*。选取不同的基矢，即不同的表象，量子态和力学量将呈现出不同的形式。

这组*基矢${ket(n)}$通常选为某力学量(完备集)$hat(F)$的正交归一本征态*，即
$
  hat(F) ket(n) = f_n ket(n)
$
如果指标$n$是离散的，完备性表达为
$
  sum_n ketbra(n) = hat(I)
$
如果指标$n$是连续的，则求和改为积分
$
  integral dd(xi) ketbra(xi) = hat(I)
$
如果指标$n$是离散的，正交归一性写为
$
  ketbra(n, n') = delta_(n,n')
$
如果指标$n$是连续的，归一到$δ$函数
$
  ketbra(xi, xi') = delta(xi - xi')
$

对于一个量子系统，这组完备的基矢就构成了一个*表象*，称为*$F$表象*(注意：$F$通常代表力学量完备集)。*利用完备性，我们就可以很方便地进入任何一个表象*。

=== 量子态

*量子态*如何进入表象？假设某个量子态为$ket(psi)$，那么在它前面作用一个单位算符，然后把单位算符用完备性关系带入。

对于*离散表象*
$
  ket(psi) & = hat(I) ket(psi) = sum_n ketbra(n) ket(psi) \
           & = sum_n ket(n) braket(n, psi) \
           & = sum_n psi_n ket(n)
$
我们把$psi_n = braket(n, psi)$称为态矢$ket(psi)$在$F$表象下的*波函数*，或者简称为波函数。对于离散表象，这个波函数实际上可以写成一个*列向量*
$
  psi = mat(psi_1, psi_2, ..., psi_n, ...)^TT
$
对于*连续表象*，同样操作
$
  ket(psi) & = hat(I) ket(psi) = integral dd(xi) ketbra(xi) ket(psi) \
           & = integral dd(xi) ket(xi) braket(xi, psi) \
           & = integral dd(xi) psi(xi) ket(xi)
$
连续表象下的波函数是一个*连续函数*
$
  psi(xi) = braket(xi, psi)
$
例如我们常见的*坐标表象$(ξ → x)$和动量表象$(ξ → p)$*。

如何计算*量子态的内积*？对于*离散表象*
$
  braket(phi, psi) & = braket(phi, hat(I), psi) = sum_n braket(phi, n) braket(n, psi) \
                   & = sum_n phi_n^* psi_n = phi^dagger psi
$
对于*连续表象*
$
  braket(phi, psi) & = braket(phi, hat(I), psi) = integral dd(xi) braket(phi, xi) braket(xi, psi) \
                   & = integral dd(xi) phi(xi)^* psi(xi)
$

#proposition(subname: [量子态进入表象])[
  量子态$ket(psi)$在$F$表象下的波函数为
  - 离散表象：$psi_n = braket(n, psi)$，波函数为列向量$psi = mat(psi_1, psi_2, ..., psi_n, ...)^TT$
  - 连续表象：$psi(xi) = braket(xi, psi)$，波函数为连续函数$psi(xi)$
]

=== 力学量

*力学量*如何进入表象？考虑力学量$hat(A)$，它是Hilbert空间的一个算符，它作用在任意态$ket(psi)$上的结果是另一个态$ket(phi)$，即
$
  hat(A) ket(psi) = ket(phi)
$
#newpara()
利用完备性口诀，灵活地在上式中合适的地方插入单位算符。对于*离散表象*，我们这么操作
$
  hat(A) ket(psi) & = hat(I) hat(A) hat(I) ket(psi) \
                  & = (sum_m ketbra(m)) hat(A) (sum_n ketbra(n)) ket(psi) \
                  & = sum_(m,n) ket(m) braket(m, hat(A), n) braket(n, psi) \
                  & = sum_m (sum_n A_(m n) psi_n) ket(m)
$
其中$A_(m n) = braket(m, hat(A), n)$称为力学量$hat(A)$在$F$表象下的*矩阵元*

另一方面，态$ket(phi)$可以用完备性口诀变为
$
  ket(phi) & = hat(I) ket(phi) = sum_m ketbra(m) ket(phi) \
           & = sum_m ket(m) braket(m, phi) = sum_m phi_m ket(m)
$
利用正交性，即展开系数相等，得到
$
  sum_n A_(m n) psi_n = phi_m
$
这实际上是一个*矩阵方程*
$
  A psi = phi
$
即：方阵$A$乘以列向量$psi$等于列向量$phi$。因此，力学量在离散表象下的波函数和矩阵元分别是列向量和方阵。

因此，在*离散表象下，力学量以矩阵的形式呈现*。若Hilbert空间是$N$维的，那么力学量就是$N × N$的方阵。任意力学量$hat(A)$在此表象下的矩阵$A$的矩阵元为
$
  A_(m n) = braket(m, hat(A), n)
$
上面推导矩阵形式的过程还可以简化为
$
  hat(A) & = hat(I) hat(A) hat(I) = sum_m ketbra(m) hat(A) sum_n ketbra(n) \
         & = sum_(m,n) ket(m) braket(m, hat(A), n) bra(n) \
         & = sum_(m,n) A_(m n) ketbra(m, n)
$
由于$hat(A)$是Hermite算符，容易证明对应矩阵$A$是Hermite矩阵，即$A_(m n) = A_(n m)^*$。

对于*连续表象*，做类似操作
$
  hat(A) ket(psi) & = hat(I) hat(A) hat(I) ket(psi) \
                  & = (integral dd(xi') ketbra(xi')) hat(A) (integral dd(xi) ketbra(xi)) ket(psi) \
                  & = integral dd(xi') dd(xi) ket(xi') braket(xi', hat(A), xi) braket(xi, psi) \
                  & = integral dd(xi') (integral dd(xi) braket(xi', hat(A) psi(xi)) ket(xi')
$
$
  ket(phi) & = hat(I) ket(phi) = integral dd(xi') ketbra(xi') ket(phi) \
           & = integral dd(xi') ket(xi') braket(xi', phi) = integral dd(xi) phi(xi) ket(xi)
$
利用基矢正交性得到
$
  integral dd(xi) braket(xi', hat(A), xi) psi(xi) = phi(xi')
$
这实际上是一个*积分方程*，其中
$
  A(xi', xi) = braket(xi', hat(A), xi)
$
称为力学量$hat(A)$在$F$表象下的*核函数*。由于基矢指标是连续的，此时算符$hat(A)$并不呈现为矩阵形式，而是呈现为*函数和微分算符*的形式。

#proposition(subname: [力学量进入表象])[
  力学量$hat(A)$在$F$表象下的表示为
  - 离散表象：矩阵元$A_(m n) = braket(m, hat(A), n)$，矩阵方程$sum_n A_(m n) psi_n = phi_m$
  - 连续表象：核函数$A(xi', xi) = braket(xi', hat(A), xi)$，积分方程$integral dd(xi) A(xi', xi) psi(xi) = phi(xi')$
]
#newpara()

#example(subname: [坐标表象下的动量算符])[
  以常见的*坐标表象*为例，$ξ → x$。考虑动量算符$hat(A) = hat(p)$，我们需要计算“矩阵元” $braket(x, hat(p), x')$。再次利用口诀，做如下计算：
  $
    braket(x, hat(p), x') & = braket(x, hat(I) hat(p) hat(I), x') \
                          & = braket(x, integral dd(p) ketbra(p) hat(p) integral dd(p') ketbra(p'), x') \
                          & = integral dd(p) integral dd(p') braket(x, p) braket(p, hat(p), p') braket(p', x') \
  $
  因为$braket(x, p)$就是*动量算符本征态在坐标表象的波函数*，所以
  $
    braket(x, p) = 1/sqrt(2 pi hbar) e^(i p x / hbar)\
    braket(p, x) = 1/sqrt(2 pi hbar) e^(- i p x / hbar)
  $
  并利用$braket(p, hat(p), p') = p' delta(p - p')$，得到
  $
    braket(x, hat(p), x') & = integral dd(p) integral dd(p') 1/(2 pi hbar) e^(i p x / hbar) p' delta(p - p') e^(- i p' x' / hbar) \
    & = integral dd(p) 1/(2 pi hbar) p e^(i p (x - x') / hbar) \
    & = - i hbar pdv(, x) delta(x - x')
  $
  其中用到了
  $
    integral dd(p) p e^(i p (x - x') / hbar) = - i hbar pdv(, x) integral dd(p) e^(i p (x - x') / hbar) = - i hbar pdv(, x) (2 pi hbar delta(x - x'))
  $
  从而
  $
    - i hbar pdv(, x) psi(x) = phi(x)
  $
  这就是方程$hat(p) ket(psi) = ket(phi)$在坐标表象的形式，表明动量算符$hat(p)$在坐标表象呈现为一个*微分算符*$- i hbar pdv(, x)$。

  同样，我们可以证明方程
  $
    hat(x) ket(psi) = ket(phi)
  $
  在*动量表象*下的形式
  $
    i hbar pdv(, p) psi(p) = phi(p)
  $
  因此，在动量表象，坐标算符$hat(x)$呈现为*微分算符*$i hbar pdv(, p)$。
]

=== 算符的不确定性

#note[
  前面的推导有循环论证的嫌疑！前面说$braket(x, p)$就是动量算符本征态在坐标表象的波函数，这一点没错，但是要得出
  $
    braket(x, p) = 1/sqrt(2 pi hbar) e^(i p x / hbar)
  $
  实际上需要事先知道动量算符的形式是
  $
    hat(p) = - i hbar pdv(, x)
  $
  只有这样，本征方程
  $
    hat(p) psi_p (x) = p psi_p (x)
  $
  的解才是*平面波*
  $
    psi_p(x) = 1/sqrt(2 pi hbar) e^((i p x) / hbar)
  $
  所以，在推导动量算符形式的过程中，平面波是一种*假设*。事实上，动量算符的形式并不是唯一的，其本征态也不一定必须是平面波形式。

  根据基本公设，正则坐标和正则动量满足正则对易关系
  $
    [hat(x), hat(p)] = i hbar
  $
  如果动量算符是如下形式
  $
    hat(p) = - i hbar pdv(, x) + f(x)
  $<p-hat>
  那么正则对易关系照样成立，其中$f(x)$可以为任意实函数。

  如何理解这样一种不确定性呢？下面我们从物理上探讨一下这个问题。如果动量算符可以取成一般形式，那么动量算符的本征方程变为
  $
    (- i hbar pdv(, x) + f(x)) psi_p (x) = p psi_p (x)
  $
  首先遇到的问题是：*动量本征态$psi_p (x)$不再是平面波*！

  不过，波函数不是直接观测量，我们需要去考察力学量的本征值和对力学量进行测量的概率分布是否发生变化。

  实际上，本征方程 @p-hat 的解可以求出来，结果是
  $
    psi_p (x) = 1/sqrt(2 pi hbar) e^((i (p x - g(x))) / hbar), g(x) = integral_a^x dd(y) f(y)
  $
  我们看到，即使动量算符多出一个尾巴$f(x)$，*其本征值不发生变化，其本征函数与平面波相比仅仅多出了一个相位因子*。

  因此，如果一个微观粒子处于某个量子态$Ψ(x)$，那么对其动量进行测量，得到的结果(包括观测值和概率分布)是不会发生变化的。对其他力学量的观测结果也不会产生影响。

  如果动量算符形式发生变化，那么描述波函数时间演化的Schrödinger方程也发生了变化。如果采用动量算符的一般形式 @p-hat ，Schrödinger方程为
  $
    i hbar pdv(, t) psi(x, t) = ( 1/(2m) (- i hbar pdv(, x) + f(x))^2 + V(x) ) psi(x, t) \
  $
  受到动量本征态讨论的启发，我们可以猜测$f(x) != 0$和$f(x) = 0$时的波函数在任意时刻也只差一个整体相位因子。这样，两种情况在物理上完全没有区别。

  事实上，如果我们定义
  $
    psi(x, t) = phi(x, t) e^((- i g(x)) / hbar), g(x) = integral_a^x dd(y) f(y)
  $
  那么经过简单计算就可以发现$phi(x, t)$满足标准形式的Schrödinger方程
  $
    i hbar pdv(, t) phi(x, t) = ( - hbar^2/(2m) pdv(, x)^2 + V(x) ) phi(x, t)
  $
  这实际上就是取$f(x)=0$时波函数满足的Schrödinger方程。

  因此，*两种情况下的波函数在任意时刻也只差一个整体相位因子，在物理上完全等价*。
]

因此，在量子力学中存在这样的不确定性，但是不影响物理结果。我们事先进行*约定(convention)*即可。

*动量算符本征态在坐标表象的波函数取为*
$
  braket(x, p) = 1/sqrt(2 pi hbar) e^((i p x) / hbar)
$
这实际上是将动量本征态的相位进行固定，即选取$g(x) = 0$。在一些教科书中干脆直接将*动量算符的矩阵元固定*，即约定
$
  braket(x, hat(x), x') = x delta(x - x')\
  braket(x, hat(p), x') = - i hbar pdv(, x) delta(x - x')
$
#newpara()
我们还可以更一般地理解这种不确定性。

量子系统的Hilbert空间是*复数线性空间*。因此，用一组完备基矢构建一个表象时，每个基矢的相位可以随意选择而不影响物理的结果。*不同的相位选取就导致了力学量的矩阵元的不确定性*。

假设力学量$hat(F)$的本征值是离散的，即
$
  hat(F) ket(n) = f_n ket(n)
$
这些正交归一本征态${|n⟩}$可以构建一个表象，称为$F$表象。在此表象下，力学量$hat(A)$呈现为矩阵形式，矩阵元为
$
  A_(m n) = braket(m, hat(A), n)
$
但是，如果$hat(n)$是$hat(F)$的本征态，那么它*乘上任意相位因子*
$e^(i θ_n)$后也是$hat(F)$的本征态。

这就意味着，我们既可以选择${ket(n)}$来构$F$表象，也可以选择${e^(i θ_n) ket(n)}$来构建$F$表象。后者，力学量$hat(A)$的矩阵元为
$
  A_(m n) & = braket(m e^(- i θ_m), hat(A), e^(i θ_n) n) \
          & = e^(- i (θ_m - θ_n)) braket(m, hat(A), n)
$
后者的矩阵元$A_(m n)$与前者相比，相差了一个相位因子
$
  e^(i (θ_n - θ_m))
$
由于这种不确定性不会改变物理结果，所以通常可以选取最简单的相位约定。

#example(subname: [Pauli矩阵])[
  一个熟悉的例子就是自旋 1/2 的 Pauli 矩阵。根据自旋角动量的性质得到，自旋算符在 $sigma_z$ 表象下的矩阵形式为
  $
    sigma_x = mat(0, e^(i alpha); e^(- i alpha), 0), sigma_y = mat(0, - i e^(i alpha); i e^(- i alpha), 0), sigma_z = mat(1, 0; 0, -1)
  $
  其中的相位$alpha$可以取任意实数。这个相位的不确定性可以理解为$sigma_z$表象的基矢$ket(arrow.t)$和$ket(arrow.b)$可以乘上任意的相位。由于这个不确定性不影响物理结果，通常我们选取$alpha = 0$，即常用的*Pauli矩阵*。
]
#newpara()
对于动量算符的不确定性，也可以认为是坐标表象的基矢$ket(x)$被乘上了任意的相位因子，即$ket(x) -> e^(i g(x))/hbar ket(x)$。相应的波函数$phi(x) = braket(x, psi)$作如下变换
$
  phi(x) -> e^(- (i g(x))/hbar) phi(x)
$
上述变换实际上是*幺正变换*，相应地，动量算符做如下变换
$
  - i hbar pdv(, x) -> e^(- (i g(x))/hbar) (- i hbar pdv(, x)) e^((i g(x))/hbar) = - i hbar pdv(, x) + f(x)
$
上式应在作用任意函数$phi(x)$的意义下理解。

=== Schrödinger方程

*Schrödinger方程*如何进入表象？考虑时间演化方程
$
  i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
$
#newpara()
对于*离散表象*，我们这么操作：
$
  i hbar pdv(, t) ket(psi(t)) & = i hbar pdv(, t) hat(I) ket(psi(t)) \
                              & = i hbar pdv(, t) sum_m ketbra(m) ket(psi(t)) \
                              & = sum_m ket(m) i hbar pdv(, t) braket(m, psi(t)) \
                              & = sum_m ket(m) i hbar pdv(, t) psi_m (t)
$
右边操作也是一样的：
$
  hat(H) ket(psi(t)) & = hat(I) hat(H) hat(I) ket(psi(t)) \
                     & = sum_m ketbra(m) hat(H) sum_n ketbra(n) ket(psi(t)) \
                     & = sum_(m,n) ket(m) braket(m, hat(H), n) braket(n, psi(t)) \
                     & = sum_m ket(m) sum_n H_(m n) psi_n (t)
$
利用正交性得到
$
  i hbar pdv(, t) psi_m (t) = sum_n H_(m n) psi_n (t)
$
所以，在*离散表象*中，Schrödinger方程呈现为*矩阵形式*
$
  i hbar pdv(, t) psi(t) = H psi(t)\
  i hbar pdv(, t) mat(psi_1(t); psi_2(t); ...) = mat(H_(11), H_(12), ...; H_(21), H_(22), ...; ..., ..., ...) mat(psi_1(t); psi_2(t); ...)
$
定态Schrödinger方程即为求解Hamilton量矩阵的本征值问题
$
  H phi = E phi\
  mat(H_(11), H_(12), ...; H_(21), H_(22), ...; ..., ..., ...) mat(phi_1; phi_2; ...) = E mat(phi_1; phi_2; ...)
$
在实际问题中，Hamilton量矩阵通常是无穷维的。不过如果我们只对较低的能级感兴趣，通常可以对Hilbert空间做一个截断，只考虑其中某个有限维的子空间。


对于*连续表象*，为了简单起见，考虑一维势阱中的粒子，Hamilton量为
$
  hat(H) = hat(p)^2/(2m) + V(hat(x))
$
要进入*坐标表象*，可如下操作：
$
  "RHS" & = hat(I) hat(H) hat(I) ket(psi(t)) \
        & = integral dd(x) ketbra(x) (hat(p)^2/(2m) + V(hat(x))) integral dd(x') ketbra(x') ket(psi(t)) \
        & = integral dd(x) dd(x') ket(x) braket(x, hat(H), x') braket(x', psi(t)) \
        & = integral dd(x) dd(x') ket(x) braket(x, hat(H), x') psi(x', t)
$
我们需要计算“矩阵元”$braket(x, hat(H), x')$。利用口诀，做如下计算
$
  braket(x, V(hat(x)), x') = V(x) delta(x - x')\
$
$
  braket(x, hat(p)^2/(2m), x') &= 1/(2m) braket(x, hat(p)^2, x') \
  & = 1/(2m) braket(x, hat(I) hat(p)^2 hat(I), x') \
  & = 1/(2m) braket(x, integral dd(p) ketbra(p) hat(p)^2 integral dd(p') ketbra(p'), x') \
  & = 1/(2m) integral dd(p) integral dd(p') braket(x, p) braket(p, hat(p)^2, p') braket(p', x') \
  & = 1/(2m) integral dd(p) integral dd(p') 1/(2 pi hbar) e^((i p x) / hbar) p'^2 delta(p - p') e^(- (i p' x') / hbar) \
  & = - 1/(2m) pdv(, x, 2) integral dd(p) 1/(2 pi hbar) e^((i p (x - x')) / hbar) \
  & = - hbar^2/(2m) pdv(, x, 2) delta(x - x')
$
Schrödinger方程左边为
$
  "LHS" & = i hbar pdv(, t) ket(psi(t)) = i hbar pdv(, t) integral dd(x) ketbra(x) ket(psi(t)) \
        & = integral dd(x) ket(x) i hbar pdv(, t) braket(x, psi(t)) = integral dd(x) ket(x) i hbar pdv(, t) psi(x, t)
$
因此在坐标表象中的Schrödinger方程为
$
  i hbar pdv(, t) psi(x, t) = ( - hbar^2/(2m) pdv(, x, 2) + V(x) ) psi(x, t)
$
定态Schrödinger方程呈现为微分方程本征值问题
$
  ( - hbar^2/(2m) pdv(, x, 2) + V(x) ) phi(x) = E phi(x)
$
要进入*动量表象*，可如下操作：
$
  "RHS" & = hat(I) hat(H) hat(I) ket(psi(t)) \
        & = integral dd(p) ketbra(p) (hat(p)^2/(2m) + V(hat(x))) integral dd(p') ketbra(p') ket(psi(t)) \
        & = integral dd(p) dd(p') ket(p) braket(p, hat(H), p') braket(p', psi(t)) \
        & = integral dd(p) dd(p') ket(p) braket(p, hat(H), p') psi(p', t)
$
我们需要计算“矩阵元”$braket(p, hat(H), p')$。利用口诀，做如下计算
$
  braket(p, hat(p)^2/(2m), p') = p^2/(2m) delta(p - p')\
$
$
  braket(p, V(hat(x)), p') & = braket(p, hat(I) V(hat(x)) hat(I), p') \
  & = braket(p, integral dd(x) ketbra(x) V(hat(x)) integral dd(x') ketbra(x'), p') \
  & = integral dd(x) integral dd(x') braket(p, x) braket(x, V(hat(x)), x') braket(x', p') \
  & = integral dd(x) integral dd(x') 1/(2 pi hbar) e^(- (i p x) / hbar) V(x) delta(x - x') e^((i p' x') / hbar) \
  & = 1/(2 pi hbar) integral dd(x) V(x) e^(- (i (p - p') x) / hbar)\
  &= V(i hbar pdv(, p)) integral dd(x) 1/(2 pi hbar) e^(- (i (p - p') x) / hbar) \
  &= V(i hbar pdv(, p)) delta(p - p')
$
因此，在*动量表象*中的Schrödinger方程为
$
  i hbar pdv(, t) psi(p, t) = ( p^2/(2m) + V(i hbar pdv(, p)) ) psi(p, t)
$
定态Schrödinger方程呈现为微分方程本征值问题
$
  ( p^2/(2m) + V(i hbar pdv(, p)) ) phi(p) = E phi(p)
$
可见，对于势阱中的量子力学问题，如果势函数比较特殊，有时采用动量表象计算更为简单。例如：线性势、$δ$函数势。

#proposition(subname: [Schrödinger方程进入表象])[
  Schrödinger方程在不同表象下的形式为
  - 离散表象：矩阵形式$i hbar pdv(, t) psi(t) = H psi(t)$，定态方程$H phi = E phi$
  - 坐标表象：微分形式$i hbar pdv(, t) psi(x, t) = ( - hbar^2/(2m) pdv(, x, 2) + V(x) ) psi(x, t)$，定态方程$( - hbar^2/(2m) pdv(, x, 2) + V(x) ) phi(x) = E phi(x)$
  - 动量表象：微分形式$i hbar pdv(, t) psi(p, t) = ( p^2/(2m) + V(i hbar pdv(, p)) ) psi(p, t)$，定态方程$( p^2/(2m) + V(i hbar pdv(, p)) ) phi(p) = E phi(p)$
]

=== 表象变换

*不同表象之间如何变换？*

考虑两个离散表象，其完备基矢分别用${ket(n)}$和${ket(alpha)}$表示，为方便起见，分别称为$A$表象和$B$表象。

将$A$表象的任意基矢进入$B$表象，得到
$
  ket(n) = sum_alpha ket(alpha) braket(alpha, n) = sum_alpha S_(alpha n) ket(alpha)
$
这里定义了变换矩阵$S$，其矩阵元为
$
  S_(alpha n) = braket(alpha, n)
$
$S$的Hermite共轭矩阵$S^dagger$的矩阵元为
$
  S^dagger_(n alpha) = (S_(alpha n))^* = braket(n, alpha)
$
我们可以计算
$
  (S S^dagger)_(alpha beta) & = sum_n S_(alpha n) S^dagger_(n beta) = sum_n braket(alpha, n) braket(n, beta) \
                            & = braket(alpha, beta) = delta_(alpha beta)
$
从而
$
  S S^dagger = hat(I)
$
同样地，可以证明$S^dagger S = hat(I)$，因此$S$是一个*幺正矩阵*。

考虑不同表象的波函数之间的变换。对于任意态$ket(psi)$，$A$表象的波函数为$psi_n = braket(n, psi)$，$B$表象的波函数为$phi_alpha = braket(alpha, psi)$。不同表象的波函数之间的变换可如下推导：
$
  phi_alpha & = braket(alpha, psi) = braket(alpha, hat(I), psi) = sum_n braket(alpha, n) braket(n, psi) \
            & = sum_n braket(alpha, n) psi_n \
            & = sum_n S_(alpha n) psi_n
$
因此两个表象的波函数通过一个幺正变换联系起来，
$
  phi = S psi\
  psi = S^dagger phi = S^(-1) phi
$
#newpara()

考虑力学量$F$的矩阵形式在不同表象之间的变换。推导如下：
$
  braket(alpha, hat(F), beta) & = braket(alpha, hat(I) hat(F) hat(I), beta) \
                              & = sum_(m,n) braket(alpha, m) braket(m, hat(F), n) braket(n, beta) \
                              & = sum_(m,n) S_(alpha m) braket(m, hat(F), n) S^dagger_(n beta)
$
所以力学量在不同表象的矩阵也通过这个幺正变换联系
$
  F_B = S F_A S^dagger\
  F_A = S^dagger F_B S
$

#newpara()

对于*连续表象*，我们来考虑坐标表象和动量表象的波函数之间的变换。从坐标表象的波函数$ψ(x)$出发，
$
  psi(x) & = braket(x, psi) = braket(x, hat(I), psi) \
         & = integral dd(p) braket(x, p) braket(p, psi) \
         & = integral dd(p) 1/sqrt(2 pi hbar) e^((i p x) / hbar) psi(p)
$
此即*Fourier变换*。逆变换是
$
  psi(p) & = braket(p, psi) = braket(p, hat(I), psi) \
         & = integral dd(x) braket(p, x) braket(x, psi) \
         & = integral dd(x) 1/sqrt(2 pi hbar) e^(- (i p x) / hbar) psi(x)
$
#newpara()

#proposition(subname: [表象变换])[
  不同表象之间的变换由幺正矩阵$S$给出：
  - 基矢变换：$ket(n) = sum_alpha S_(alpha n) ket(alpha), ket(alpha) = sum_n S^dagger_(n alpha) ket(n), S_(alpha n) = braket(alpha, n)$
  - 波函数变换：$phi = S psi$，$psi = S^dagger phi$
  - 力学量矩阵变换：$F_B = S F_A S^dagger$，$F_A = S^dagger F_B S$
]
#newpara()

表象变换有如下性质：
+ *表象变换不改变力学量的本征值*。假设算符$hat(F)$在$A$表象中的本征方程为
  $
    F_A psi_A = f_n^A psi_A
  $
  则在$B$表象中有
  $
    F_B psi_B = S F_A S^dagger S psi_A = S F_A psi_A = f_n^A S psi_A = f_n^A psi_B
  $
  进一步，可以证明*表象变换不改变力学量的平均值*。

+ *表象变换不改变力学量的对易关系*。假设在$A$表象中有对易关系
  $
    [F_A, G_A] = X_A
  $
  则在$B$表象中有
  $
    [F_B, G_B] & = S [F_A, G_A] S^dagger = S X_A S^dagger = X_B
  $
+ *表象变换不改变力学量的迹*。任意算符$hat(F)$的迹定义为它在某个表象中的对角矩阵元的求和，即
  $
    tr(hat(F)) = sum_n braket(n, hat(F), n)
  $
  假设在$A$表象中求迹为$tr F_A$，则在$B$表象中有
  $
    tr F_B = tr(S F_A S^dagger) = tr(F_A S^dagger S) = tr F_A
  $
  由于在任意表象求迹结果都一样，因此也可以在连续表象求迹，此时求和变为积分
  $
    tr(hat(F)) = integral dd(xi) braket(xi, hat(F), xi)
  $

  #example(subname: [平衡态量子统计物理])[
    平衡态量子统计物理中，配分函数的基本公式是
    $
      cal(Z) = tr(e^(- beta hat(H))), beta = 1/(k_B T)
    $
    由于在任意表象求迹结果都一样，我们可以在任何表象中求出算符$e^(- beta hat(H))$的对角矩阵元，然后求和。以谐振子为例，要计算出结果，最简单的方法是采用能量表象或者粒子数表象
    $
      cal(Z) = sum_n braket(n, e^(- beta hat(H)), n) = sum_n e^(- beta E_n)
    $
    但是，对于相互作用多体系统，能量本征值一般无法全部解析求出，这个办法就抓瞎了。一个想法是直接在坐标表象中求迹，即
    $
      cal(Z) = integral dd(x) braket(x, e^(- beta hat(H)), x)
    $
    然后将“虚时”间隔$[0, β]$无限细分，继续采用口诀，就可以推出配分函数的虚时路径积分形式。这样就可以进入*有限温度量子场论*。
  ]
