#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

= 二次量子化

本章讨论处理*全同粒子体系*的普遍方法：二次量子化。它是进一步学习量子场论和量子多体理论（凝聚态场论）的基础。
- 全同粒子系统最主要的特点是*微观粒子的不可分辨性*，即粒子不可以区分，也就意味着对粒子进行编号并无实质意义。
- 在本科量子力学中，我们采用“*一次量子化*”的方法描述全同粒子系统。在这种描述中，我们需要对粒子进行编号，但是系统的态矢或波函数要满足*对称性要求*：对于交换任意两个粒子的操作必须是*对称或反对称*的。这样，虽然我们对粒子编了号，但是这些它们仍然是不可分辨的。
- 在一次量子化的描述中，我们仍然将每个粒子的*坐标当作动力学变量*（自由度，如果有内部自由度也要考虑，比如自旋），*用单粒子的波函数来构造全同粒子体系的波函数*。我们需要*对粒子进行编号，然后进行对称化或反对称化操作*。
- 有没有不对粒子进行编号的描述方法呢？最简单的方法就是，将*每个单粒子态上的粒子数交代清楚，全同粒子系统的状态就完全确定了*。这就是“*二次量子化*”。
- 为何称为“二次量子化”？
  - 若单粒子态采用坐标表象，其图像可以理解为：将单粒子波函数$psi(vb(x))$及其复共轭$psi^*(vb(x))$看成是一个场，即所谓的*Schrödinger场*。
  - 我们将$psi(vb(x))$和$psi^*(vb(x))$看成是系统的正则坐标，求出对应的正则动量，然后采用*正则量子化*方法，引入相应的正则对易关系。
  - 二次量子化以后，波函数成了*场算符*$hat(psi)(vb(x))$和$hat(psi)^dagger (vb(x))$，即所谓的*湮灭算符和产生算符*，坐标$vb(x)$不再是算符，而只是一个参数。
  - 多粒子体系的态矢可以利用产生算符作用在真空态上构造出来，从而无须对粒子编号，自动满足对称性要求。
  - 二次量子化描述中，体系的自由度可以是无穷大。

== 全同粒子系统

考虑一个两粒子系统。假设单粒子的Hilbert空间的基矢取为一组完备力学量的共同本征态，记为${ket(k)}$，$k$代表这组量子数。两粒子系统的Hilbert空间是两个单粒子空间的直积(张量积)。其基矢最简单的取法就是取为两个空间的基矢的直积
$
  ket(k')_1 times.o ket(k'')_2 = ket(k') ket(k'')
$
将两个粒子的状态交换，则成为
$
  ket(k'')_1 times.o ket(k')_2 = ket(k'') ket(k')
$
若$k' != k''$，由于两个基矢是正交的，所以两种状态是可区分的。如果两个粒子是全同粒子，这两种状态显然都不满足对称性要求。下面，讨论全同粒子的情形。

如果我们对这两种两粒子状态进行测量，我们将对其中一个粒子得到$k'$，对另一个粒子得到$k''$。但是，由此我们并不能确定系统的状态到底是$ket(k')ket(k''), ket(k'')ket(k')$，还是两种状态的任意叠加
$
  c_1 ket(k') ket(k'') + c_2 ket(k'') ket(k')
$
这被称为*交换简并*。这样，与单粒子情形不同，我们无法通过一组力学量完备集对应的本征值（量子数）来确定系统的量子态。

定义*置换算符*$hat(P)_12$
$
  hat(P)_12 ket(k') ket(k'') = ket(k'') ket(k')
$
很显然
$
  hat(P)_12 = hat(P)_21 , hat(P)_12^2 = hat(I)
$
#newpara()

考虑最简单的情形，单粒子态由某个力学量$hat(A)$的本征值刻画，则
$
  hat(A)_1 ket(a') ket(a'') = a' ket(a') ket(a'')\
  hat(A)_2 ket(a') ket(a'') = a'' ket(a') ket(a'')
$
第一个式子两边作用$hat(P)_12$，得到
$
  hat(P)_12 hat(A)_1 hat(P)_12^(-1) hat(P)_12 ket(a') ket(a'') = a' hat(P)_12 ket(a') ket(a'')\
  => hat(P)_12 hat(A)_1 hat(P)_12^(-1) ket(a'') ket(a') = a' ket(a'') ket(a')
$
由此我们得到置换算符对力学量的作用是
$
  hat(P)_12 hat(A)_1 hat(P)_12^(-1) = hat(A)_2,\
  hat(P)_12 hat(A)_2 hat(P)_12^(-1) = hat(A)_1
$
即：改变力学量中的粒子编号指标。

考虑两个全同粒子系统的Hamilton量，一般性地可以写为
$
  hat(H) = hat(H)_1 + hat(H)_2 + hat(H)_"int"
$
相互作用满足
$
  hat(P)_12 hat(H)_"int" hat(P)_12^(-1) = hat(H)_"int" => hat(P)_12 hat(H) hat(P)_12^(-1) = hat(H)
$
例如，两粒子系统Hamilton量的具体形式可以写为
$
  hat(H) = hat(vb(p))_1^2/(2m) + hat(U)_"ext" (vb(x)_1) + hat(vb(p))_2^2/(2m) + hat(U)_"ext" (vb(x)_2) + hat(V)_"int" (abs(vb(x)_1 - vb(x)_2))
$
这意味着
$
  [hat(H), hat(P)_12] = 0
$
$hat(P)_12$是守恒量。由$hat(P)_12^2 = 1$可知其本征值为$+1$和$-1$。所以，如果体系在初始时刻处于本征值为$+1$和$-1$的本征态，这种状态将不随时间变化。之前所学的量子力学并没有对全同粒子的态矢做出要求，因此全同粒子体系的态矢可以是$hat(P)_12$的线性组合态。但全同性原理要求，全同粒子的态矢必须是$hat(P)_12$的本征态，要么是对称的，要么是反对称的。

以上讨论可以直接推广到任意$N$粒子系统，置换算符$hat(P)_(i j)$定义为
$
  hat(P)_(i j) ket(k') ket(k'') ... ket(k^(i)) ... ket(k^(j)) ... ket(k^(N)) = ket(k') ket(k'') ... ket(k^(j)) ... ket(k^(i)) ... ket(k^(N))
$
*全同性原理*则对全同粒子体系的量子态有着更高的要求，即对于任意两个粒子的置换操作，要么是对称的，要么是反对称的
$
  hat(P)_(i j) ket(N "identical bosons") = + ket(N "identical bosons")\
  hat(P)_(i j) ket(N "identical fermions") = - ket(N "identical fermions")
$
很显然，直接将单粒子态直积起来构造的基矢不满足对称性要求，多粒子系统的Hilbert空间中存在大量不满足对称性要求的态矢。因此，我们需要构造对称化或反对称化的基矢。

在非相对论量子力学的框架下，*全同性原理是其第五大公设*。在相对论的量子场论中，自旋-统计定理表明，自旋为整数的粒子必须服从对称化原理，自旋为半整数的粒子必须服从反对称化原理。因此，全同性原理与粒子的自旋密切相关。

#example(subname: [两粒子系统])[
  对于两粒子系统，假设单粒子空间基矢为$ket(k')$和$ket(k'')$。若不考虑全同性，则基矢为
  $
    ket(k') ket(k''), ket(k'') ket(k'), ket(k') ket(k'), ket(k'') ket(k'')
  $
  若为全同粒子，则对称化的基矢为
  $
    ket(k') ket(k'), ket(k'') ket(k''), 1/sqrt(2) (ket(k') ket(k'') + ket(k'') ket(k'))
  $
  反对称化的基矢为
  $
    1/sqrt(2) (ket(k') ket(k'') - ket(k'') ket(k'))
  $
]

#example(subname: [三粒子系统])[
  对于三粒子系统，假设单粒子空间基矢为$ket(k')$、$ket(k'')$和$ket(k''')$。若考虑三个不同的状态构造的对称化或反对称化的基矢，则为
  $
    ket(k' k'' k''')_plus.minus = & 1/sqrt(6) (ket(k') ket(k'') ket(k''') plus.minus ket(k') ket(k''') ket(k'') \
                                  & + ket(k'') ket(k''') ket(k') plus.minus ket(k'') ket(k') ket(k''') \
                                  & + ket(k''') ket(k') ket(k'') plus.minus ket(k''') ket(k'') ket(k'))
  $
  其中正号和负号代表对称化和反对称化的基矢。当然，还存在有两个单粒子态相同的基矢，例如
  $
    ket(k' k' k'')_+ = 1/sqrt(3) (ket(k') ket(k') ket(k'') + ket(k') ket(k'') ket(k') + ket(k'') ket(k') ket(k'))
  $
]

== Fock 空间

对于一个$N$粒子系统，原则上可以构造对称化或反对称化的基矢，从而构造出满足全同性原理的Hilbert空间。但是如果粒子数很大，这种基矢的数量必然也很大，彼此之间又没有明显联系，理论上处理起来很不方便。因此前人发展了另一种处理方法，即*二次量子化*。

=== Fock 空间的构造与产生湮灭算符

设单粒子力学量$hat(K)$的本征方程为
$
  hat(K) ket(k_i) = k_i ket(k_i)
$
即: 本征值为${k_i}$，本征态为$ket(k_i)$。*二次量子化方法，就是用单粒子本征态上的粒子数分布来描述全同粒子体系的状态。*定义多粒子系统态矢为
$
  ket(n_1\, n_2\, ... n_i\, ...), n_i = 0, 1, 2, ...
$
$n_i$是处于第$i$个单粒子本征态$ket(k_i)$的粒子数。

由这些态矢构成的矢量空间称为*Fock空间*。它实际上是*所有可能的粒子数的对称化（或反对称化）Hilbert空间的直和*。因此，它既可以处理粒子数不变的系统，也可以处理粒子数可变的系统。
$
  cal(H) = cal(H)_1 times.o cal(H)_2 times.o ... times.o cal(H)_N
$
一次量子化构造的满足对称和反对称的子空间分别为
$
  cal(H)_"B"^((n)), cal(H)_"F"^((n))
$
这都是对于$N$粒子系统而言的。Fock空间的基矢为
$
  ket(n_1\, n_2\, ... n_i\, ...), sum_i n_i = N != "const"
$
这样给出的Hilbert空间不再局限于固定粒子数$N$，而是包含了所有可能的粒子数。从而
$
  cal(H)_"Fock" = sum_n cal(H)^((n))_"B or F"
$

#newpara()
Fock 空间中有两个特殊的态矢。第一个是*真空态*
$
  ket(0\, 0\, ... 0\, ...) = ket(vb(0))
$
即：所有单粒子态上都没有粒子。第二个是*单粒子态*
$
  ket(0\, 0\, ... n_i = 1\, ... 0\, ...) = ket(k_i)
$
即：一个粒子处于本征值为$k_i$的本征态$ket(k_i)$。这两个态都已经归一化，即
$
  braket(vb(0)) = 1, braket(k_i) = 1
$
#newpara()
多粒子态如何满足全同性原理呢？既然全同粒子系统的状态用粒子数分布表示，则可以定义改变粒子数的算符：*产生算符和湮灭算符*（某种意义下也可称为场算符）。

#note[
  这里的产生湮灭算符和当时谐振子中的产生湮灭算符有着本质的区别：谐振子中的产生湮灭算符是改变谐振子能级的，而这里的产生湮灭算符是改变粒子数的。但我们如此命名的方式是其对易关系是一致的。后面我们会说这样的算符也可以被称为*场算符*，因为它们可以看成是场的量子化。
]

定义*产生算符*$hat(a)^dagger_i$，对Fock空间中任意态矢的作用为
$
  hat(a)^dagger_i ket(n_1\, n_2\, ... n_i\, ...) prop ket(n_1\, n_2\, ... n_i + 1\, ...)
$
其中归一化常数待定。假设$hat(a)^dagger_i$作用在真空态上，得到单粒子态
$
  hat(a)^dagger_i ket(vb(0)) = ket(k_i)
$
由此得到
$
  1 = braket(k_i) = braket(vb(0), hat(a)_i hat(a)^dagger_i, vb(0)) = braket(vb(0), hat(a)_i, k_i)\
  => hat(a)_i ket(k_i) = ket(vb(0))
$
所以产生算符$hat(a)^dagger_i$的Hermite共轭是*湮灭算符*，定义其对Fock空间中任意态矢的作用为
$
  hat(a)_i ket(n_1\, n_2\, ... n_i\, ...) prop ket(n_1\, n_2\, ... n_i - 1\, ...)
$
尤其是
$
  hat(a)_i ket(vb(0)) = 0, hat(a)_i ket(k_j) = 0 (i!-j)
$
湮灭算符对单粒子态的作用可合并写为
$
  hat(a)_i ket(k_j) = delta_(i j) ket(vb(0))
$
产生算符和湮灭算符可以很自然地用来表示多粒子系统的力学量。

=== 表象变换

*表象变换*：Fock态中粒子占据的单粒子态，原则上*可以用任意单粒子力学量的完备本征态*。之前用的是单粒子力学量$hat(K)$的本征态，如果换成别的力学量，例如$hat(G)$，在新的 Fock 态中定义的产生湮灭算符与原来的有何关系？Fock的Hilbert空间是一致的，我们认为是其基矢的表示方式不同，也就是表象不同。

假设$hat(G)$的本征方程为
$
  hat(G) ket(g_i) = g_i ket(g_i), i = 1, 2, ...
$
即：本征值为${g_i}$，本征态为${ket(g_i)}$。利用这组单粒子态定义的Fock态为
$
  ket(m_1\, m_2\, ... m_i\, ...), m_i = 0, 1, 2, ...
$
同样，定义产生算符和湮灭算符
$
  hat(b)_i^dagger ket(m_1\, m_2\, ... m_i\, ...) prop ket(m_1\, m_2\, ... m_i + 1\, ...)\
  hat(b)_i ket(m_1\, m_2\, ... m_i\, ...) prop ket(m_1\, m_2\, ... m_i - 1\, ...)
$
#newpara()

*单粒子态的表象变换*：(K 表象$<->$G 表象)
$
  ket(k_i) = sum_m ket(g_m) braket(g_m, k_i) = sum_m S_(m i) ket(g_m)\
  S_(m i) = braket(g_m, k_i)
$
#newpara()
*产生湮灭算符的表象变换*：利用表象变换矩阵的幺正性
$
  (S^dagger S)_(i j) = (S S^dagger)_(i j) = delta_(i j)
$
得到($S_(i j)^dagger = (S^dagger)_(i j) = (S^*)_(j i)$)
$
  hat(a)_i ket(k_j) &= delta_(i j) ket(vb(0)) = sum_m S_(i m)^dagger S_(m j) ket(vb(0)) = sum_(m l) S_(i m)^dagger S_(l j) delta_(m l) ket(vb(0))\
  &= sum_(m l) S_(i m)^dagger S_(l j) hat(b)_m ket(g_l) = sum_m S_(i m)^dagger hat(b)_m (sum_l S_(l j) ket(g_l)) = sum_m S_(i m)^dagger hat(b)_m ket(k_j)\
$
由于$ket(k_j)$是任意单粒子态，所以
$
  hat(a)_i = sum_m braket(k_i, g_m) hat(b)_m\
  hat(a)_i^dagger = sum_m hat(b)_m^dagger braket(g_m, k_i)
$
逆变换为
$
  hat(b)_m = sum_i braket(g_m, k_i) hat(a)_i\
  hat(b)_m^dagger = sum_i hat(a)_i^dagger braket(k_i, g_m)
$
以上结果表明：*利用某个（单粒子）表象定义的产生算符（湮灭算符），可以用别的表象定义的产生算符（湮灭算符）的线性叠加表示出来，叠加系数就是（单粒子）表象变换的矩阵元*。我们还发现*单粒子态和产生湮灭算符的对应关系*
$
  ket(k_i) <-> hat(a)_i^dagger, ket(k_i) <-> hat(a)_i^dagger, ketbra(k_i, k_j) <-> hat(a)_i^dagger hat(a)_j
$
意思是
$
  ket(k_i) = sum_m braket(g_m, k_i) ket(g_m) <-> hat(a)_i^dagger = sum_m braket(g_m, k_i) hat(b)_m^dagger
$

=== 全同性原理

接下来考虑*全同性原理*，即Fock空间中的任意态矢，对于任意两个粒子的置换操作，要么是对称的，要么是反对称的。对于Fock空间中的任意态矢，考虑先在$ket(k_i)$态加一个粒子，再在$ket(k_j)$态加一个粒子。*置换等价于改变这两个操作的顺序*。所以
$
  hat(a)^dagger_i hat(a)^dagger_j ket(n_1\, n_2\, ...\, n_i\, ...) = plus.minus hat(a)^dagger_j hat(a)^dagger_i ket(n_1\, n_2\, ...\, n_i\, ...)
$
其中$+$号对应Bose子， $-$号对应Fermi子。态矢是任意的，所以
$
  hat(a)^dagger_i hat(a)^dagger_j minus.plus hat(a)^dagger_j hat(a)^dagger_i = 0
$
即（反对易子${A, B} ≡ A B + B A$）
$
  hat(a)^dagger_i hat(a)^dagger_j - hat(a)^dagger_j hat(a)^dagger_i = [hat(a)^dagger_i, hat(a)^dagger_j] = 0 "   " & "Bosons"\
  hat(a)^dagger_i hat(a)^dagger_j + hat(a)^dagger_j hat(a)^dagger_i = {hat(a)^dagger_i, hat(a)^dagger_j} = 0 "   " & "Fermions"\
$
此关系的Hermite共轭给出
$
  [hat(a)_i, hat(a)_j] = 0 "   " & "Bosons" \
  {hat(a)_i, hat(a)_j} = 0 "   " & "Fermions" \
$
#newpara()
*Pauli不相容原理*自动包含在了其中：对于Fermi子
$
  hat(a)^dagger_i hat(a)^dagger_i = 0
$
即两个费米子不能占据同一个量子态。

同样，当$i != j$时，由全同性原理得到
$
  hat(a)_i hat(a)^dagger_j ket(n_1\, n_2\, ...\, n_i\, ...) = plus.minus hat(a)^dagger_j hat(a)_i ket(n_1\, n_2\, ...\, n_i\, ...)
$
注意：当$i = j$时，此式不一定成立，例如
$
  hat(a)_i hat(a)^dagger_i ket(vb(0)) = ket(vb(0)) != hat(a)^dagger_i hat(a)_i ket(vb(0)) = 0
$
从而
$
  hat(a)_i hat(a)^dagger_j minus.plus hat(a)^dagger_j hat(a)_i = 0, i != j
$
即
$
  hat(a)_i hat(a)^dagger_j - hat(a)^dagger_j hat(a)_i = [hat(a)_i, hat(a)^dagger_j] = 0 "   " & "Bosons" \
  hat(a)_i hat(a)^dagger_j + hat(a)^dagger_j hat(a)_i = {hat(a)_i, hat(a)^dagger_j} = 0 "   " & "Fermions" \
$
当$i = j$时上述对易（反对易）关系结果如何？我们假设
$
  hat(a)_i hat(a)^dagger_i - hat(a)^dagger_i hat(a)_i = [hat(a)_i, hat(a)^dagger_i] = 1 "   " & "Bosons" \
  hat(a)_i hat(a)^dagger_i + hat(a)^dagger_i hat(a)_i = {hat(a)_i, hat(a)^dagger_i} = 1 "   " & "Fermions" \
$
#newpara()
下面将证明：对于Bose子，这将使得粒子数$n_i$取值为$0, 1, 2, ...$；而对于Fermi子，这将使得粒子数$n_i$只能取值为$0$和$1$。

根据以上讨论，当$i = j$时有
$
  "Bosons": [hat(a)_i, hat(a)^dagger_i] = 1, [hat(a)_i, hat(a)_i] = 0, [hat(a)^dagger_i, hat(a)^dagger_i] = 0\
  "Fermions": {hat(a)_i, hat(a)^dagger_i} = 1, {hat(a)_i, hat(a)_i} = 0, {hat(a)^dagger_i, hat(a)^dagger_i} = 0
$
对于给定的$i$，可定义*粒子数算符*
$
  hat(N)_i = hat(a)^dagger_i hat(a)_i
$
不论对于Bose子还是Fermi子，都有
$
  [hat(N)_i, hat(a)^dagger_i] = hat(a)^dagger_i, [hat(N)_i, hat(a)_i] = - hat(a)_i
$
类似*谐振子的代数解法*，求解$hat(N)_i$的本征方程
$
  hat(N)_i ket(n_i) = n_i ket(n_i)
$
利用对易关系得到
$
  hat(N)_i hat(a)^dagger_i ket(n_i) = ([hat(N)_i, hat(a)^dagger_i] + hat(a)^dagger_i hat(N)_i) ket(n_i) = (n_i + 1) hat(a)^dagger_i ket(n_i)\
  hat(N)_i hat(a)_i ket(n_i) = ([hat(N)_i, hat(a)_i] + hat(a)_i hat(N)_i) ket(n_i) = (n_i - 1) hat(a)_i ket(n_i)
$
所以$hat(a)^dagger_i ket(n_i)$和$hat(a)_i ket(n_i)$也是$hat(N)_i$的本征态，本征值增加或减少1。

*对于Bose子*，$hat(a)_i hat(a)^dagger_i = hat(a)^dagger_i hat(a)_i + 1$，所以
$
  hat(a)_i ket(n_i) = alpha ket(n_i - 1) => braket(n_i, hat(a)_i^dagger hat(a)_i, n_i) = abs(alpha)^2 => n_i = abs(alpha)^2 => hat(a)_i ket(n_i) = sqrt(n_i) ket(n_i - 1)\
  hat(a)^dagger_i ket(n_i) = beta ket(n_i + 1) => braket(n_i, hat(a)_i hat(a)_i^dagger, n_i) = abs(beta)^2 => n_i + 1 = abs(beta)^2 => hat(a)^dagger_i ket(n_i) = sqrt(n_i + 1) ket(n_i + 1)
$
利用湮灭算符$hat(a)_i$可不断使得本征值减小，再利用
$
  n_i = braket(n_i, hat(a)^dagger_i hat(a)_i, n_i) >= 0
$
可知$n_i$最小值为0，所以$n_i$取值为$0, 1, 2, ...$。利用上述结果还可以得到
$
  hat(a_i) ket(0) = 0\
  ket(n_i) = (hat(a)^dagger_i)^(n_i)/sqrt(n_i !) ket(0)
$
这与之前谐振子的结果完全一致，所以又称为*Bose型谐振子*。

*对于Fermi子*，同样
$
  hat(a)_i ket(n_i) = alpha ket(n_i - 1) => braket(n_i, hat(a)_i^dagger hat(a)_i, n_i) = abs(alpha)^2 => n_i = abs(alpha)^2 => hat(a)_i ket(n_i) = sqrt(n_i) ket(n_i - 1)\
  hat(a)^dagger_i ket(n_i) = beta ket(n_i + 1) => braket(n_i, hat(a)_i hat(a)_i^dagger, n_i) = abs(beta)^2 => 1 - n_i = abs(beta)^2 => hat(a)^dagger_i ket(n_i) = sqrt(1 - n_i) ket(n_i + 1)
$
利用湮灭算符$hat(a)_i$不断使得本征值减小，可知$n_i$最小值为0。然后利用产生算符$hat(a)^dagger_i$不断使得本征值增大，得到
$
  hat(a)^dagger_i ket(0) = ket(1) => hat(a)^dagger_i hat(a)^dagger_i ket(0) = hat(a)^dagger_i ket(1) = 0
$
所以，$n_i$最大值为$1$。

因此，对于Fermi子，$n_i$只能取$0$和$1$两个值。上述结果可总结为
$
  hat(a)_i ket(0) = 0, hat(a)_i ket(1) = ket(0)\
  hat(a)^dagger_i ket(0) = ket(1), hat(a)^dagger_i ket(1) = 0\
$
这个系统也被称为*Fermi型谐振子*。对于Fermi子，我们也可以将本征态$ket(n)_i$写为与Bose子一样的形式，利用上述结果还可以得到
$
  ket(n_i) = (hat(a)^dagger_i)^(n_i)/sqrt(n_i !) ket(0) (n_i = 0, 1)
$
推广到Fock空间中的态矢，不论Bose子还是Fermi子，都可写为
$
  ket(n_1\, n_2\, ... n_i\, ...) = product_i (hat(a)^dagger_i)^(n_i)/sqrt(n_i !) ket(vb(0))
$
现在可以确定产生算符和湮灭算符作用在Fock空间中任意态矢得到的结果。

- *对于Bose子*：
  $
    hat(a)^dagger_i ket(n_1\, n_2\, ... n_i\, ...) = sqrt(n_i + 1) ket(n_1\, n_2\, ... n_i + 1\, ...)\
    hat(a)_i ket(n_1\, n_2\, ... n_i\, ...) = sqrt(n_i) ket(n_1\, n_2\, ... n_i - 1\, ...)
  $
- *对于Fermi子*：
  $
    hat(a)^dagger_i ket(n_1\, n_2\, ... n_i\, ...) = xi sqrt(1 - n_i) ket(n_1\, n_2\, ... n_i + 1\, ...)\
    hat(a)_i ket(n_1\, n_2\, ... n_i\, ...) = xi sqrt(n_i) ket(n_1\, n_2\, ... n_i - 1\, ...)
  $
  其中因子$xi$来源于交换产生的$(-1)$的乘积
  $
    xi = (-1)^(sum_(j < i) n_j)
  $
  由于Fermi子的$n_i$只能取$0, 1$，上述关系又可明确写为
  $
    hat(a)^dagger_i ket(n_1\, n_2\, ... n_i\, ...) = xi delta_(n_i, 0) ket(n_1\, n_2\, ... 1_i\, ...)\
    hat(a)_i ket(n_1\, n_2\, ... n_i\, ...) = xi delta_(n_i, 1) ket(n_1\, n_2\, ... 0_i\, ...)
  $
很显然，不论对Bose子还是Fermi子
$
  hat(N)_i ket(n_1\, n_2\, ... n_i\, ...) = n_i ket(n_1\, n_2\, ... n_i\, ...)
$
定义*总粒子数算符*
$
  hat(N) = sum_i hat(N)_i = sum_i hat(a)^dagger_i hat(a)_i
$
得
$
  hat(N) ket(n_1\, n_2\, ... n_i\, ...) = (sum_i n_i) ket(n_1\, n_2\, ... n_i\, ...)
$
#newpara()
产生算符和湮灭算符满足的对易或反对易关系可总结如下：
- *对于Bose子*：
  $
    [hat(a)_i, hat(a)_j] = 0, [hat(a)^dagger_i, hat(a)^dagger_j] = 0,[hat(a)_i, hat(a)^dagger_j] = delta_(i j)
  $
- *对于Fermi子*：
  $
    {hat(a)_i, hat(a)_j} = 0, {hat(a)^dagger_i, hat(a)^dagger_j} = 0, {hat(a)_i, hat(a)^dagger_j} = delta_(i j)
  $
Bose子满足的对易关系等价于要求体系的量子态对于任意两个粒子的置换是对称的，即*Bose-Einstein统计*；Fermi子满足的反对易关系等价于要求体系的量子态对于任意两个粒子的置换是反对称的，即*Fermi-Dirac统计*。

因此，定义了Fock空间和产生湮灭算符，我们就可以描述全同粒子构成的多粒子系统。

可以证明，*表象变换不改变产生湮灭算符的对易关系*。在以上讨论中，单粒子态采用$hat(K)$表象的基矢。若采用$hat(G)$表象，则可计算
$
  &[hat(b)_i, hat(b)_j] = sum_(k l) S_(i k) S_(j l) [hat(a)_k, hat(a)_l] = 0\
  &[hat(b)^dagger_i, hat(b)^dagger_j] = sum_(k l) S_(i k)^dagger S_(j l)^dagger [hat(a)^dagger_k, hat(a)^dagger_l] = 0\
  &[hat(b)_i, hat(b)^dagger_j] = sum_(k l) S_(i k) S_(j l)^dagger [hat(a)_k, hat(a)^dagger_l] = sum_(k l) S_(i k) S_(l j)^dagger delta_(k l) = delta_(i j)
$
此外，还可以证明：*总粒子数算符的形式与表象无关*
$
  hat(N) = sum_i hat(a)^dagger_i hat(a)_i = sum_(i k l) S_(k i) hat(b)^dagger_k S_(i l)^dagger hat(b)_l = sum_l hat(b)^dagger_l hat(b)_l
$

== Fock空间中的力学量

=== Fock空间中的力学量

构造了满足全同性原理的Fock空间以后，我们需要研究如何构造Fock空间中的力学量，即*力学量的二次量子化形式*。一般来说，$N$个粒子构成的全同粒子系统的力学量有如下几种类型：
- *单体算符*：可以写成$N$个单体力学量之和，例如粒子的动能、在外场中的势能；
- *两体算符*：可以写成$N(N - 1)$个两体力学量之和，例如一对粒子之间的相互作用势能；
- *三体算符*：可以写成$N(N - 1)(N - 2)$个三体力学量之和，例如核子系统中的三体相互作用；
- ...
力学量的二次量子化形式，实际上就是求出上述算符在Fock空间中的形式，将其用产生算符和湮灭算符表示出来。

=== 单体算符

设单粒子力学量$hat(K)$的本征方程为
$
  hat(K) ket(k_i) = k_i ket(k_i), i = 1, 2, ...
$
Fock空间的基矢可以用本征态${ket(k_i)}$上的粒子数分布构造
$
  ket(Psi) = ket(n_1\, n_2\, ... n_i\, ...)
$
对于$N$粒子系统，与单粒子力学量$hat(K)$对应的单体算符为
$
  hat(cal(K)) = sum_(alpha = 1)^N hat(K)_alpha
$
容易看到，$ket(Psi)$是$hat(K)$的本征态，本征值为
$
  K_"total" = sum_i n_i k_i
$
所以*单体算符$hat(cal(K))$可以用产生算符和湮灭算符表示为*
$
  hat(cal(K)) = sum_i k_i hat(N)_i = sum_i k_i hat(a)^dagger_i hat(a)_i
$

#example(subname: [自由粒子])[
  自由粒子，单粒子Hamilton量为
  $
    hat(H) = hat(vb(p))^2/(2m), hat(H) ket(vb(p)) = vb(p)^2/(2m) ket(vb(p))
  $
  二次量子化的Hamilton量为
  $
    hat(cal(H)) = sum_vb(p) vb(p)^2/(2m) hat(a)^dagger_vb(p) hat(a)_vb(p)
  $
  这里通常采用箱归一化，将体系看成边长为$L$的立方体，从而动量取值是离散化的。计算完再令 $L -> oo$。
]

#example(subname: [外势场中的粒子])[
  外势场中的粒子，例如冷原子实验中囚禁在磁光阱中的原子、光晶格中的原子、固体晶格中的电子等。单粒子Hamilton量为
  $
    hat(H) = hat(vb(p))^2/(2m) + hat(U) (hat(vb(r)))
  $
  假设单粒子Hamilton量的本征值问题可以求解
  $
    hat(H) ket(phi_i) = epsilon_i ket(phi_i)
  $
  则多粒子系统二次量子化的Hamilton量为
  $
    hat(cal(H)) = sum_i epsilon_i hat(a)^dagger_i hat(a)_i
  $
]

以上例子表明：*在利用单粒子力学量自身表象构造的Fock空间中，该单粒子力学量对应的单体算符的二次量子化形式是对角的*。

如果单粒子态变了，那么单体算符的形式也会发生变化。考虑另一个单粒子力学量$hat(G)$，本征方程为
$
  hat(G) ket(g_i) = g_i ket(g_i), i = 1, 2, ...
$
Fock空间的基矢也可以用本征态$ket(g_i)$上的粒子数分布构造
$
  ket(Psi') = ket(m_1\, m_2\, ... m_i\, ...)
$
对应的产生(湮灭)算符为$hat(b)_i, hat(b)^dagger_i$。问题：采用这组基矢对应的产生湮灭算符，与单粒子力学量$hat(K)$对应的单体算符$hat(cal(K))$的形式如何？

根据之前的讨论，新旧产生和湮灭算符之间存在表象变换
$
  hat(a)_i = sum_j braket(k_i, g_j) hat(b)_j\
  hat(a)_i^dagger = sum_j hat(b)_j^dagger braket(g_j, k_i)
$
带入得到
$
  hat(cal(K)) & = sum_i k_i sum_(m n) hat(b)_m^dagger braket(g_m, k_i) braket(k_i, g_n) hat(b)_n \
              & = sum_(m n) hat(b)_m^dagger (sum_i k_i braket(g_m, k_i) braket(k_i, g_n)) hat(b)_n \
              & = sum_(m n) hat(b)_m^dagger bra(g_m) (hat(K) sum_i ketbra(k_i)) ket(g_n) hat(b)_n \
$
即
$
  hat(cal(K)) = sum_(m n) hat(b)_m^dagger bra(g_m) hat(K) ket(g_n) hat(b)_n = sum_(m n) K_(m n) hat(b)_m^dagger hat(b)_n
$
这就是单体算符最*一般性的二次量子化形式*。

#example(subname: [外势场中的粒子])[
  外势场中的粒子，单粒子Hamilton量为
  $
    hat(H) = hat(vb(p))^2/(2m) + hat(U) (hat(vb(r)))
  $
  若单粒子态仍取为动量本征态$ket(vb(p))$，则二次量子化的Hamilton量为
  $
    hat(cal(H)) &= sum_(vb(p) vb(p)') bra(vb(p)) hat(H) ket(vb(p)') hat(a)^dagger_vb(p) hat(a)_vb(p')\
    &= sum_(vb(p)) vb(p)^2/(2m) hat(a)^dagger_vb(p) hat(a)_vb(p) + sum_(vb(p) vb(p)') bra(vb(p)) hat(U) ket(vb(p)') hat(a)^dagger_vb(p) hat(a)_vb(p')
  $
  可见，由于$ket(vb(p))$不是$hat(H)$的本征态，单体哈密顿量中出现了非对角项。
]

=== 两体算符

在多粒子系统中，还存在两体算符($α, β$：粒子编号)
$
  hat(cal(V)) = sum_(alpha < beta) hat(V)_(alpha beta) = 1/2 sum_(alpha != beta) hat(V)_(alpha beta)
$
假设对任意两粒子的$K ⊗ K$表象，有如下本征方程
$
  hat(V)_(alpha beta) ket(k_i)_alpha ket(k_j)_beta = V_(i j) ket(k_i)_alpha ket(k_j)_beta
$
由于$hat(V)_(alpha beta)$是Hermite算符，所以$V_(i j)$是实数，且$V_(i j) = V_(j i)$。所以，在$K ⊗ K$表象中，$hat(V)_(alpha beta)$是对角的。

两体算符可以用产生湮灭算符表示为
$
  hat(cal(V)) = 1/2 sum_(i != j) V_(i j) hat(N)_i hat(N)_j + 1/2 sum_i V_(i i) hat(N)_i (hat(N)_i - 1) \
$

#proof[
  利用单粒子的$K$表象构造Fock空间的基矢
  $
    ket(Psi) = ket(n_1\, n_2\, ... n_i\, ...)
  $
  考虑两体算符$hat(cal(V))$对 Fock 空间基矢的作用。易知，$ket(Psi)$是$hat(cal(V))$的本征态，本征值为
  $
    V_"total" = 1/2 sum_(i != j) V_(i j) n_i n_j + 1/2 sum_i V_(i i) n_i (n_i - 1)
  $
  另一方面，考虑原始定义的两体算符对$ket(Psi)$的作用。在$ket(Psi)$定义的Fock态中，有 $n_i$个粒子处于$ket(k_i)$ 态，$n_j$个粒子处于$ket(k_j)$ 态。利用假设可知，$ket(Psi)$是$hat(cal(V))$的本征态，本征值与上式相同。因此，两体算符的二次量子化形式即为上述形式。
]
#newpara()
可以将$i != j$和$i = j$的两项合并写为
$
  hat(cal(V)) = 1/2 sum_(i j) V_(i j)(hat(N)_i hat(N)_j - delta_(i j) hat(N)_i) = 1/2 sum_(i j) V_(i j) hat(Pi)_(i j)
$
其中
$
  hat(Pi)_(i j) = hat(N)_i hat(N)_j - delta_(i j) hat(N)_i
$
称为*对分布算符*。利用对易或反对易关系，可计算
$
  hat(Pi)_(i j) &= hat(N)_i hat(N)_j - delta_(i j) hat(N)_i = hat(a)^dagger_i hat(a)_i hat(a)^dagger_j hat(a)_j - delta_(i j) hat(a)^dagger_i hat(a)_i \
  &= hat(a)^dagger_i (delta_(i j) plus.minus hat(a)^dagger_j hat(a)_i) hat(a)_j - delta_(i j) hat(a)^dagger_i hat(a)_i \
  &= plus.minus hat(a)^dagger_i hat(a)^dagger_j hat(a)_i hat(a)_j\
  &= (plus.minus) (plus.minus) hat(a)^dagger_i hat(a)^dagger_j hat(a)_j hat(a)_i\
  &= hat(a)^dagger_i hat(a)^dagger_j hat(a)_j hat(a)_i
$
所以，两体算符可以写为更为简洁的形式
$
  hat(cal(V)) = 1/2 sum_(i j) V_(i j) hat(a)^dagger_i hat(a)^dagger_j hat(a)_j hat(a)_i
$
可知
$
  V_(i j) = braket(k_i\, k_j, hat(V), k_i\, k_j) = braket(k_i, braket(k_j, hat(V), k_i), k_j)
$
通常简记为
$
  V_(i j) = braket(i j, hat(V), i j)
$
则
$
  hat(cal(V)) = 1/2 sum_(i j) braket(i j, hat(V), i j) hat(a)^dagger_i hat(a)^dagger_j hat(a)_j hat(a)_i
$
#newpara()

在以上形式中，选取了特殊的$K$表象，其中$hat(V)$是对角的。一般表象下两体算符的形式，可以通过表象变换得到。考虑另一个表象，即$hat(G)$表象，利用表象变换可以计算
$
  hat(cal(V)) & = 1/2 sum_(i j) V_(i j) hat(a)^dagger_i hat(a)^dagger_j hat(a)_j hat(a)_i \
  & = 1/2 sum_(i j) V_(i j) sum_(m n p q) (S_(m i) hat(b)^dagger_m) (S_(n j) hat(b)^dagger_n) (S_(j p)^dagger hat(b)_p) (S_(i q)^dagger hat(b)_q) \
  & = 1/2 sum_(i j) sum_(m n p q) V_(i j) braket(g_m, k_i) braket(g_n, k_j) braket(k_j, g_p) braket(k_i, g_q) hat(b)^dagger_m hat(b)^dagger_n hat(b)_p hat(b)_q \
  & = 1/2 sum_(m n p q) braket(m n, hat(V), p q) hat(b)^dagger_m hat(b)^dagger_n hat(b)_p hat(b)_q\
$
其中
$
    & sum_(i j) V_(i j) braket(g_m, k_i) braket(g_n, k_j) braket(k_j, g_p) braket(k_i, g_q) = braket(m n, hat(V), p q) \
  = & sum_(i j) V_(i j) (bra(g_m) bra(h_n))(ket(k_i) ket(k_j)) (bra(k_j) bra(k_i))(ket(g_p) ket(g_q)) \
  = & sum_(i j) bra(g_m) bra(h_n) hat(V) (ket(k_i) ket(k_j) bra(k_j) bra(k_i))ket(g_p) ket(g_q) \
  = & braket(g_m, braket(g_n, hat(V), g_p), g_q) \
  = & braket(m n, hat(V), p q)
$

#example(subname: [粒子之间的两体相互作用])[
  粒子之间的两体相互作用。单粒子取坐标表象，有
  $
    hat(V)_(alpha beta) ket(vb(x))_alpha ket(vb(x)')_beta = V(vb(x), vb(x)') ket(vb(x))_alpha ket(vb(x)')_beta
  $
  其中$V(vb(x), vb(x)') = V(abs(vb(x) - vb(x)'))$。由于$hat(V)$在坐标表象是对角的，所以两体相互作用的二次量子化形式为
  $
    hat(cal(V)) = 1/2 integral dd(vb(x), 3) integral dd(vb(x)', 3) V(abs(vb(x) - vb(x)')) hat(a)^dagger (vb(x)) hat(a)^dagger (vb(x)') hat(a)(vb(x)') hat(a)(vb(x))
  $
  下面将看到，这里的产生湮灭算符，实际上就是场算符，即
  $
    hat(cal(V)) = 1/2 integral dd(vb(x), 3) integral dd(vb(x)', 3) V(abs(vb(x) - vb(x)')) hat(psi)^dagger (vb(x)) hat(psi)^dagger (vb(x)') hat(psi)(vb(x)') hat(psi)(vb(x))
  $
  若采用动量表象，则$hat(V)$不是对角的，需要采用一般的形式。将在后面的实例中进行推导。
]

== 实例：电子气

作为一个计算实例，考虑金属中的电子气。做为简化模型，取边长为$L$的立方体$(V = L^3)$，$N$个电子在其中运动。由于整体是电中性的（正离子），将正电荷看成是在立方体中*均匀分布*的背景。

只考虑Coulomb相互作用，系统的Hamilton量为
$
  hat(H) = hat(H)_e + hat(H)_"b" + hat(H)_(e"b")\
$
其中，$hat(H)_e$包含电子的动能和电子之间的Coulomb相互作用
$
  hat(H)_e = sum_i hat(vb(p))_i^2/(2m) + 1/2 e^2 sum_(i != j) e^(- mu abs(hat(vb(x))_i - hat(vb(x))_j))/abs(hat(vb(x))_i - hat(vb(x))_j)
$
这里，我们假设Coulomb力被屏蔽了，引入了屏蔽质量$mu$，计算的最后令$mu -> 0$。

$hat(H)_"b"$为正电背景的相互作用能
$
  hat(H)_"b" = 1/2 integral dd(vb(x), 3) integral dd(vb(x)', 3) rho(vb(x)) rho(vb(x)') e^(- mu abs(vb(x) - vb(x)'))/(abs(vb(x) - vb(x)'))
$
根据模型假设，$ρ(x) = N/V$，计算得到
$
  hat(H)_"b" = 1/2 e^2 (N/V)^2 integral dd(vb(x)', 3) integral dd(vb(x), 3) e^(- mu abs(vb(x)))/(abs(vb(x))) = 1/2 e^2 N^2/V (4pi)/(mu^2) \
$
$hat(H)_(e"b")$为电子与正电背景的相互作用能
$
  hat(H)_(e"b") = - e^2 sum_i integral dd(vb(x), 3) rho(vb(x)) e^(- mu abs(hat(vb(x))_i - vb(x)))/(abs(hat(vb(x))_i - vb(x))) = - e^2 N^2/V (4pi)/(mu^2) \
$
注意：要用国际单位制，做代换$e^2 -> e^2/(4 pi epsilon_0)$。

所以，体系的Hamilton顿量可写为
$
  hat(H) = sum_i hat(vb(p))_i^2/(2m) + 1/2 e^2 sum_(i != j) e^(- mu abs(hat(vb(x))_i - hat(vb(x))_j))/abs(hat(vb(x))_i - hat(vb(x))_j) - 1/2 e^2 N^2/V (4pi)/(mu^2)
$
虽然最后一项当$µ -> 0$时发散，但是后面会看到，它将被精确抵消。接下来，将此Hamilton量写成二次量子化的形式。

将单电子态取为动量$hat(vb(p))$和自旋$sigma_z$的共同本征态$ket(vb(k) lambda), lambda=plus.minus$。考虑由于周期性边界条件，$vb(k)$的取值是离散的($vb(k) = vb(p)/hbar$)
$
  vb(k) = (2pi)/L vb(n), vb(n) in ZZ^3
$
动量本征态在坐标表象的波函数为
$
  braket(vb(x), vb(k) lambda) = 1/sqrt(V) e^(i vb(k) dot vb(x)) chi_lambda
$
正交归一和完备性条件
$
  braket(vb(k)' lambda', vb(k) lambda) = delta_(vb(k)', vb(k)) delta_(lambda', lambda)\
  sum_(vb(k) lambda) ketbra(vb(k) lambda) = 1
$
#newpara()
*电子动能项*：这项对应单体算符，根据前面的公式，得到
$
  sum_i hat(vb(p))_i^2/(2m) -> hat(cal(K))& = sum_(vb(k) lambda) sum_(vb(k)' lambda') braket(vb(k) lambda, hat(vb(p))^2/(2m), vb(k)' lambda') hat(a)^dagger_(vb(k) lambda) hat(a)_(vb(k)' lambda') \
  & = sum_(vb(k) lambda) (hbar^2 vb(k)^2)/(2m) hat(a)^dagger_(vb(k) lambda) hat(a)_(vb(k) lambda)
$
#newpara()
*电子间相互作用项*：这项对应两体算符，根据前面的公式，得到
$
  hat(cal(V)) = sum_(vb(k)_1 lambda_1) sum_(vb(k)_2 lambda_2) sum_(vb(k)_3 lambda_3) sum_(vb(k)_4 lambda_4) braket(vb(k)_1 lambda_1\, vb(k)_2 lambda_2, hat(V), vb(k)_3 lambda_3\, vb(k)_4 lambda_4) times hat(a)^dagger_(vb(k)_1 lambda_1) hat(a)^dagger_(vb(k)_2 lambda_2) hat(a)_(vb(k)_4 lambda_4) hat(a)_(vb(k)_3 lambda_3)
$
其中的矩阵元可如下计算
$
  &braket(vb(k)_1 lambda_1\, vb(k)_2 lambda_2, hat(V), vb(k)_3 lambda_3\, vb(k)_4 lambda_4) \
  = &integral dd(vb(x)_1, 3) integral dd(vb(x)_2, 3) integral dd(vb(x)_3, 3) integral dd(vb(x)_4, 3) braket(vb(k)_1 lambda_1, vb(x)_1) braket(vb(k)_2 lambda_2, vb(x)_2) bra(vb(x)_1\, vb(x)_2) hat(V) ket(vb(x)_3\, vb(x)_4) braket(vb(x)_3, vb(k)_3 lambda_3) braket(vb(x)_4, vb(k)_4 lambda_4) \
  = &integral dd(vb(x), 3) integral dd(vb(x)', 3) braket(vb(k)_1 lambda_1, vb(x)) braket(vb(k)_2 lambda_2, vb(x)') hat(V)(vb(x), vb(x)') braket(vb(x), vb(k)_3 lambda_3) braket(vb(x)', vb(k)_4 lambda_4)\
  = &e^2/V^2 integral dd(vb(x), 3) integral dd(vb(x)', 3) e^(i (vb(k)_3 - vb(k)_1) dot vb(x)) e^(i (vb(k)_4 - vb(k)_2) dot vb(x)') e^(- mu abs(vb(x) - vb(x)'))/(abs(vb(x) - vb(x)')) chi_(lambda_1)^dagger chi_(lambda_3) chi_(lambda_2)^dagger chi_(lambda_4)\
  = &e^2/V^2 delta_(lambda_1, lambda_3) delta_(lambda_2, lambda_4) integral dd(vb(y), 3) e^(- i (vb(k)_1 + vb(k)_2 - vb(k)_3 - vb(k)_4) dot vb(y)) integral dd(vb(z), 3) e^(- mu abs(vb(z)))/(abs(vb(z))) e^(- i (vb(k)_1 - vb(k)_3) dot vb(z))\
$
第二步中做了变量代换：$vb(y) = vb(x)'$,$vb(z) = vb(x) - vb(x)'$。定义动量转移 $vb(q) = vb(k)_1 - vb(k)_3$，完成积分得到
$
  &braket(vb(k)_1 lambda_1\, vb(k)_2 lambda_2, hat(V), vb(k)_3 lambda_3\, vb(k)_4 lambda_4) \
  = &e^2/V delta_(lambda_1, lambda_3) delta_(lambda_2, lambda_4) delta_(vb(k)_1 + vb(k)_2, vb(k)_3 + vb(k)_4) (4pi)/(mu^2 + vb(q)^2)
$
所以，电子间相互作用项的二次量子化形式为
$
  hat(cal(V)) = e^2/(2V) sum_(lambda_1 lambda_2) sum_(vb(k)_1 vb(k)_2 vb(k)_3 vb(k)_4) delta_(vb(k)_1 + vb(k)_2, vb(k)_3 + vb(k)_4) (4pi)/(mu^2 + (vb(k)_1 - vb(k)_3)^2) times hat(a)^dagger_(vb(k)_1 lambda_1) hat(a)^dagger_(vb(k)_2 lambda_2) hat(a)_(vb(k)_4 lambda_2) hat(a)_(vb(k)_3 lambda_1)
$
令$vb(k)_3 = vb(k)$, $vb(k)_4 = vb(p)$，则$vb(k)_1 = vb(k) + vb(q)$。当$vb(k)_1 + vb(k)_2 = vb(k)_3 + vb(k)_4$时$vb(k)_2 = vb(p) - vb(q)$。所以，上式又可以写为
$
  hat(cal(V)) = e^2/(2V) sum_(lambda_1 lambda_2) sum_(vb(k) vb(p) vb(q)) (4pi)/(mu^2 + vb(q)^2) times hat(a)^dagger_(vb(k) + vb(q) lambda_1) hat(a)^dagger_(vb(p) - vb(q) lambda_2) hat(a)_(vb(p) lambda_2) hat(a)_(vb(k) lambda_1)
$
其中动量转移为零$vb(q) = 0$的项为
$
  e^2/(2V) sum_(lambda_1 lambda_2) sum_(vb(k) vb(p)) (4pi)/(mu^2) hat(a)^dagger_(vb(k) lambda_1) hat(a)^dagger_(vb(p) lambda_2) hat(a)_(vb(p) lambda_2) hat(a)_(vb(k) lambda_1) &= e^2/(2V) (4pi)/(mu^2) sum_(vb(k) lambda_1) sum_(vb(p) lambda_2) hat(N)_(vb(k) lambda_1) hat(a)^dagger_(vb(k) lambda_1) hat(a)_(vb(k) lambda_1) (hat(a)^dagger_(vb(p) lambda_2) hat(a)_(vb(p) lambda_2) - delta_(vb(k), vb(p)) delta_(lambda_1, lambda_2))\
  & = e^2/(2V) (4pi)/(mu^2) (hat(N)^2 - hat(N))
$
由于讨论固定粒子数的系统， $hat(N)$可用总粒子数$N$代替。因此，第一项与
$
  hat(H) = sum_i hat(vb(p))_i^2/(2m) + 1/2 e^2 sum_(i != j) e^(- mu abs(hat(vb(x))_i - hat(vb(x))_j))/abs(hat(vb(x))_i - hat(vb(x))_j) - 1/2 e^2 N^2/V (4pi)/(mu^2)
$
中的背景项精确抵消。第二项表示平均每个粒子的能量为$- (2 π e^2)/(µ^2 V)$，当体积$V → ∞$时这项消失。所以，$vb(q) = 0$的项没有贡献(正好抵消背景带来的发散项)。

去掉$vb(q) = 0$的项后，可以令屏蔽质量$mu = 0$而不出现奇异性。最后得到二次量子化的Hamilton量为
$
  hat(cal(H)) = hat(cal(H))_0 + hat(cal(H))_"int"
$
其中
$
  hat(cal(H))_0 &= sum_(vb(k) lambda) (hbar^2 vb(k)^2)/(2m) hat(a)^dagger_(vb(k) lambda) hat(a)_(vb(k) lambda)\
  hat(cal(H))_"int" &= e^2/(2V) sum_(lambda_1 lambda_2) sum_(vb(k) vb(p) (vb(q) != 0)) (4pi)/(vb(q)^2) times hat(a)^dagger_(vb(k) + vb(q) lambda_1) hat(a)^dagger_(vb(p) - vb(q) lambda_2) hat(a)_(vb(p) lambda_2) hat(a)_(vb(k) lambda_1)
$
#newpara()

求解上页给出的Hamilton量的本征值问题是非常困难的。可以用*微扰论*计算其基态能量，即将$hat(cal(H))_"int"$看成微扰，体系基态的零级近似就是*自由Fermi气体*。为了进一步计算基态能量的微扰修正，定义两个长度。
- 第一个是电子间的平均距离$r_0$ ，定义为
  $
    n = N/V = ((4 pi)/3 r_0^3)^(-1) => r_0 = (3/(4 pi n))^(1/3)
  $
- 第二个是*Bohr半径* $a_0$，定义为
  $
    a_0 = (hbar^2)/(m e^2)
  $
可以用Bohr半径作为长度标度，定义无量纲参数
$
  r_s = r_0/a_0
$
#newpara()
零级近似：在零级近似下，体系的基态就是$hat(cal(H))_0$的基态，即*自由Fermi气体*，其能量就是基态能量的零级近似$E^((0))$。自由Fermi气体基态可以表示为
$
  ket(G) = product_(abs(vb(k)) < k_F, lambda) hat(a)^dagger_(vb(k) lambda) ket(0)
$
其中$k_F$为Fermi动量（波矢），总粒子数$N$可如下计算
$
  N = sum_(vb(k) lambda) braket(G, hat(a)^dagger_(vb(k) lambda) hat(a)_(vb(k) lambda), G) = 2 sum_(abs(vb(k)) < k_F) theta(k_F - k) = 2 integral dd(vb(k), 3)/(2 pi)^3 theta(k_F - k) = (V)/(3 pi^2) k_F^3\
  => k_F = (3 pi^2 N/2)^(1/3) = ((9 pi)/4)^(1/3) 1/r_0
$
基态能量的零级近似$E^((0))$可如下计算
$
  E^((0)) &= sum_(vb(k) lambda) (hbar^2 vb(k)^2)/(2m) braket(G, hat(a)^dagger_(vb(k) lambda) hat(a)_(vb(k) lambda), G) \
  &= 2 sum_(abs(vb(k)) < k_F) (hbar^2 vb(k)^2)/(2m) = 2 integral dd(vb(k), 3)/(2 pi)^3 (hbar^2 vb(k)^2)/(2m) theta(k_F - k) \
  &= V/(pi^2) (hbar^2)/(2m) integral_0^(k_F) dd(k) k^4 = V/(pi^2) (hbar^2)/(2m) (k_F^5)/5 \
  &= 3/5 (hbar^2 k_F^2)/(2m) N = e^2/(2 a_0) N 3/5 ((9 pi)/4)^(2/3) 1/r_s^2
$

一级近似：基态能量的一阶修正为
$
  E^((1)) &= braket(G, hat(cal(H))_"int", G) \
  &= e^2/(2V) sum_(lambda_1 lambda_2) sum_(vb(k) vb(p) (vb(q) != 0)) (4pi)/(vb(q)^2) times braket(G, hat(a)^dagger_(vb(k) + vb(q) lambda_1) hat(a)^dagger_(vb(p) - vb(q) lambda_2) hat(a)_(vb(p) lambda_2) hat(a)_(vb(k) lambda_1), G)\
$
对角矩阵元中的 4 个算符必须两两配对，矩阵元才可能不为零。只有两种情况
- $vb(k) + vb(q) = vb(p)$, $vb(p) - vb(q) = vb(l)$, $lambda_1 = lambda_2$
- $vb(k) + vb(q) = vb(k)$, $vb(p) - vb(q) = vb(p)$
第二种情况要求$vb(q) = 0$，不在求和范围内，因此不用考虑。

对于第一种情况，矩阵元为
$
  delta_(vb(k) + vb(q), vb(p))delta_(lambda_1, lambda_2) times braket(G, hat(a)^dagger_(vb(k)+vb(q) lambda_1) hat(a)^dagger_(vb(k) lambda_1) hat(a)_(vb(k)+vb(q) lambda_1) hat(a)_(vb(k) lambda_1), G) &= - delta_(vb(k) + vb(q), vb(p)) delta_(lambda_1, lambda_2) times braket(G, hat(N)_(vb(k)+vb(q) lambda_1) hat(N)_(vb(k) lambda_2), G)\
  &= - delta_(vb(k) + vb(q), vb(p)) delta_(lambda_1, lambda_2) theta(k_F - abs(vb(k) + vb(q))) theta(k_F - abs(vb(k)))
$
带入计算基态能量的一阶修正
$
  E^((1)) &= - e^2/(2V) sum_(lambda_1) sum_(vb(k) vb(q)) (4pi)/(vb(q)^2) theta(k_F - abs(vb(k) + vb(q))) theta(k_F - abs(vb(k))) \
  &= - e^2 V integral dd(vb(k), 3)/(2 pi)^3 integral dd(vb(q), 3)/(2 pi)^3 (4pi)/(vb(q)^2) theta(k_F - abs(vb(k) + vb(q))) theta(k_F - abs(vb(k)))\
$
做积分变量代换
$
  vb(k) -> vb(P) = vb(k) + 1/2 vb(q)
$
得到
$
  E^((1)) = - e^2 V integral dd(vb(P), 3)/(2 pi)^3 integral dd(vb(q), 3)/(2 pi)^3 (4pi)/(vb(q)^2) theta(k_F - abs(vb(P) + 1/2 vb(q))) theta(k_F - abs(vb(P) - 1/2 vb(q)))
$
完成积分后得到
$
  E^((1)) = - e^2/(2 a_0) N (3)/(2 pi) ((9 pi)/4)^(1/3) 1/r_s
$
所以，计算到一级近似的基态能量为
$
  E/N = e^2/(2 a_0) [3/5 ((9 pi)/4)^(2/3) 1/r_s^2 - (3)/(2 pi) ((9 pi)/4)^(1/3) 1/r_s]
$

== Schrödinger场的量子化

采用任意的单粒子表象都可以建立多粒子系统的Fock空间，不同的表象之间通过表象变换联系起来。

=== 单粒子表象与场算符

设单粒子力学量$hat(K)$的本征方程为($K$表象)
$
  hat(K) ket(k_i) = k_i ket(k_i), i = 1, 2, ...
$
Fock空间的基矢可以用本征态${ket(k_i)}$上的粒子数分布构造
$
  ket(Psi) = ket(n_1\, n_2\, ... n_i\, ...)
$
产生湮灭算符为$hat(a)_i, hat(a)^dagger_i$。当然，我们也可以采用单粒子的*坐标表象*构造Fock态$(i → vb(x)_i)$
$
  ket(Psi') = ket(n(vb(x)_1)\, n(vb(x)_2)\, ... n(vb(x)_i)\, ...)
$
将坐标表象中的产生湮灭算符记为$hat(psi)^dagger (vb(x))$和$hat(psi)(vb(x))$。根据表象变换关系得到
$
  hat(a)^dagger_i = integral dd(vb(x), 3) hat(psi)^dagger (vb(x)) braket(vb(x), k_i)\
  hat(a)_i = integral dd(vb(x), 3) braket(k_i, vb(x)) hat(psi) (vb(x))
$
其逆变换为
$
  hat(psi)^dagger (vb(x)) = sum_i hat(a)^dagger_i braket(k_i, vb(x))\
  hat(psi) (vb(x)) = sum_i braket(vb(x), k_i) hat(a)_i
$
上述变换中的内积$braket(vb(x), k_i)$即为$ket(k_i)$态在坐标表象中的波函数
$
  phi_i (vb(x)) = braket(vb(x), k_i)
$
所以，上述变换又可写为
$
  hat(psi)^dagger (vb(x)) = sum_i hat(a)^dagger_i phi_i^* (vb(x))\
  hat(psi) (vb(x)) = sum_i hat(a)_i phi_i (vb(x))
$
*在形式上，这与将任意波函数用完备本征函数进行展开完全一致*。

利用$hat(a)^dagger_i$和$hat(a)_i$满足的对易或反对易关系，可以计算$hat(psi)^dagger (vb(x))$和$hat(psi) (vb(x))$满足的对易或反对易关系
$
  [hat(psi) (vb(x)), hat(psi) (vb(x)')]_(plus.minus) = sum_(i j) braket(vb(x), k_i) braket(vb(x)', k_j) [hat(a)_i, hat(a)_j]_(plus.minus) = 0\
  [hat(psi)^dagger (vb(x)), hat(psi)^dagger (vb(x)')]_(plus.minus) = sum_(i j) braket(k_i, vb(x)) braket(k_j, vb(x)') [hat(a)^dagger_i, hat(a)^dagger_j]_(plus.minus) = 0\
$
其中下标正号($[..., ...]_+$)代表对易关系(Bose子)，负号($[..., ...]_-$)代表反对易关系(Fermi子)。最后一个对易或反对易关系可计算为
$
  [hat(psi) (vb(x)), hat(psi)^dagger (vb(x)')]_(plus.minus) = sum_(i j) braket(vb(x), k_i) braket(k_j, vb(x)') [hat(a)_i, hat(a)^dagger_j]_(plus.minus) = sum_i phi_i (vb(x)) phi_i^* (vb(x)') = delta(vb(x) - vb(x)')
$
#newpara()
采用坐标表象，多粒子系统的单体算符$hat(cal(K))$的形式为
$
  hat(cal(K)) = integral dd(vb(x), 3) integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)) bra(vb(x)) hat(K) ket(vb(x)') hat(psi) (vb(x)')\
$
对于*动能算符*$hat(T) = hat(vb(p))^2/(2m)$，有
$
  hat(cal(T)) &= integral dd(vb(x), 3) integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)) bra(vb(x)) hat(vb(p))^2/(2m) ket(vb(x)') hat(psi) (vb(x)')\
  &= - (hbar^2)/(2m) integral dd(vb(x), 3) hat(psi)^dagger (vb(x)) nabla^2 hat(psi) (vb(x))\
  &= (hbar^2)/(2m) integral dd(vb(x), 3) nabla hat(psi)^dagger (vb(x)) dot nabla hat(psi) (vb(x))\
$
*粒子在外场中的势能*$hat(U)(x)$，可进一步计算为
$
  hat(cal(U)) &= integral dd(vb(x), 3) integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)) bra(vb(x)) hat(U) ket(vb(x)') hat(psi) (vb(x)')\
  &= integral dd(vb(x), 3)integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)) U(vb(x)) delta(vb(x) - vb(x)') hat(psi) (vb(x)')\
  &= integral dd(vb(x), 3) hat(psi)^dagger (vb(x)) U(vb(x)) hat(psi) (vb(x))\
$
*粒子之间的两体相互作用*$hat(V)$对应的两体算符，在坐标表象中的形式为
$
  hat(cal(V)) = 1/2 integral dd(vb(x)_1, 3) integral dd(vb(x)_2, 3) integral dd(vb(x)_3, 3) integral dd(vb(x)_4, 3) hat(psi)^dagger (vb(x)_1) hat(psi)^dagger (vb(x)_2) bra(vb(x)_1\, vb(x)_2) hat(V) ket(vb(x)_3\, vb(x)_4) hat(psi) (vb(x)_4) hat(psi) (vb(x)_3)\
$
对于通常的两体相互作用，有
$
  braket(vb(x)_1\, vb(x)_2, hat(V), vb(x)_3\, vb(x)_4) = V(vb(x)_3, vb(x)_4) delta(vb(x)_1 - vb(x)_3) delta(vb(x)_2 - vb(x)_4)
$
带入得到
$
  hat(cal(V)) = 1/2 integral dd(vb(x), 3) integral dd(vb(x)', 3) V(vb(x), vb(x)') hat(psi)^dagger (vb(x)) hat(psi)^dagger (vb(x)') hat(psi) (vb(x)') hat(psi) (vb(x))
$
多粒子系统的Hamilton量一般可以写为
$
  hat(H) = sum_i (hat(vb(p))_i^2/(2m) + hat(U)(hat(vb(x))_i)) + 1/2 sum_(i != j) hat(V)(hat(vb(x))_i, hat(vb(x))_j)
$
采用坐标表象，其二次量子化形式为
$
  hat(cal(H)) &= hat(cal(T)) + hat(cal(U)) + hat(cal(V))\
  &= integral dd(vb(x), 3) hat(psi)^dagger (vb(x)) (- (hbar^2)/(2m) nabla^2 + U(vb(x))) hat(psi) (vb(x)) + 1/2 integral dd(vb(x), 3) integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)) hat(psi)^dagger (vb(x)') V(vb(x), vb(x)') hat(psi) (vb(x)') hat(psi) (vb(x))\
$
总粒子数算符为
$
  hat(N) &= sum_i hat(a)_i^dagger hat(a)_i = integral dd(vb(x), 3) integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)) delta(vb(x) - vb(x)') hat(psi) (vb(x)')\
  &= integral dd(vb(x), 3) hat(psi)^dagger (vb(x)) hat(psi) (vb(x)) = integral dd(vb(x), 3) hat(n)(vb(x))
$
其中
$
  hat(n)(vb(x)) = hat(psi)^dagger (vb(x)) hat(psi) (vb(x))
$
称为粒子数密度算符。

利用对易关系
$
  [hat(psi) (vb(x)), hat(psi)^dagger (vb(x)')]_(plus.minus) = delta(vb(x) - vb(x)'), [hat(psi) (vb(x)), hat(psi) (vb(x)')]_(plus.minus) = 0, [hat(psi)^dagger (vb(x)), hat(psi)^dagger (vb(x)')]_(plus.minus) = 0
$
可以证明，不论对于Bose子还是Fermi子都有
$
  [hat(N), hat(psi)^dagger (vb(x))] = hat(psi)^dagger (vb(x)), [hat(N), hat(psi) (vb(x))] = - hat(psi) (vb(x))
$
所以，对于上式给出的Hamilton量，有
$
  [hat(H), hat(N)] = 0
$
即：*系统的总粒子数是一个守恒量*。这正是我们预期的结果：非相对论多粒子系统，总粒子数守恒。

=== 时间演化

在Schrödinger绘景中，多粒子系统的态矢$ket(Ψ(t))$的演化仍然满足Schrödinger方程
$
  i hbar pdv(, t) ket(Ψ(t)) = hat(cal(H)) ket(Ψ(t))
$
不过，$ket(Ψ(t))$是定义在 Fock 空间中的。当粒子数很大时，采用Schrödinger绘景描述时间演化变得很困难。所以，采用Heisenberg绘景更为方便。在Heisenberg绘景中，产生和湮灭算符变为
$
  hat(psi)^(dagger"H") (vb(x), t) = e^(i hat(cal(H)) t/hbar) hat(psi)^dagger (vb(x)) e^(- i hat(cal(H)) t/hbar)\
  hat(psi)^("H") (vb(x), t) = e^(i hat(cal(H)) t/hbar) hat(psi) (vb(x)) e^(- i hat(cal(H)) t/hbar)
$
由于时间宗量可以表示出Heisenberg绘景，下面略去上标H。由于Hamilton量不含时，所以不变。用$hat(psi) (vb(x), t)$和$hat(psi)^dagger (vb(x), t)$写出来的形式为
$
  hat(cal(H)) &= hat(cal(H))_0 + hat(cal(H))_"I"\
  hat(cal(H))_0 &= integral dd(vb(x), 3) hat(psi)^dagger (vb(x), t) (- (hbar^2)/(2m) nabla^2 + U(vb(x))) hat(psi) (vb(x), t)\
  hat(cal(H))_"I" &= 1/2 integral dd(vb(x), 3) integral dd(vb(x)', 3) hat(psi)^dagger (vb(x), t) hat(psi)^dagger (vb(x)', t) V(vb(x), vb(x)') hat(psi) (vb(x)', t) hat(psi) (vb(x), t)
$
由于绘景变换不改变对易或反对易关系，在Heisenberg绘景中，产生湮灭算符仍满足如下的等时对易或反对易关系
$
  [hat(psi) (vb(x), t), hat(psi)^dagger (vb(x)', t)]_(plus.minus) = delta(vb(x) - vb(x)'), [hat(psi) (vb(x), t), hat(psi) (vb(x)', t)]_(plus.minus) = 0, [hat(psi)^dagger (vb(x), t), hat(psi)^dagger (vb(x)', t)]_(plus.minus) = 0
$
在Heisenberg绘景中，算符$hat(psi) (vb(x), t)$满足的运动方程为
$
  pdv(, t) hat(psi) (vb(x), t) = 1/(i hbar) [hat(psi) (vb(x), t), hat(cal(H))]
$
计算对易关系，不论对于Bose子还是Fermi子都得到
$
  [hat(psi) (vb(x), t), hat(cal(H))_0] = (- (hbar^2)/(2m) nabla^2 + U(vb(x))) hat(psi) (vb(x), t)\
  [hat(psi) (vb(x), t), hat(cal(H))_"I"] = integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)', t) V(vb(x), vb(x)') hat(psi) (vb(x)', t) hat(psi) (vb(x), t)
$
所以，运动方程可以写为
$
  i hbar pdv(, t) hat(psi) (vb(x), t) = (- (hbar^2)/(2m) nabla^2 + U(vb(x))) hat(psi) (vb(x), t) \ + integral dd(vb(x)', 3) hat(psi)^dagger (vb(x)', t) V(vb(x), vb(x)') hat(psi) (vb(x)', t) hat(psi) (vb(x), t)
$
如果忽略粒子间的相互作用，则运动方程变为
$
  i hbar pdv(, t) hat(psi) (vb(x), t) = (- (hbar^2)/(2m) nabla^2 + U(vb(x))) hat(psi) (vb(x), t)
$
此方程在形式上与单粒子波函数$psi(vb(x), t)$满足的Schrödinger方程完全一致。

这并不是一种巧合，而是暗示我们：*可以把单粒子波函数$psi(vb(x), t)$看成是一个经典场，对此经典场进行量子化的手续等价于本章讨论的二次量子化方法*，而且给出了“二次量子化”的含义。

=== Schrödinger场

考虑单粒子波函数满足的Schrödinger方程
$
  i hbar pdv(, t) psi(vb(x), t) = (- (hbar^2)/(2m) nabla^2 + U(vb(x))) psi(vb(x), t)
$
若将$psi(vb(x), t)$看成是一个经典场，即“*Schrödinger场*”，则*Schrödinger方程即为此经典场的运动方程*。

与经典力学系统一样，经典场系统也存在一个Lagrange量，例如电磁场的Lagrange量密度就可以写为
$
  cal(L) = 1/2 (vb(E)^2 - vb(B)^2) = - 1/2 F_(mu nu) F^(mu nu)
$
利用最小作用量原理就可以得到Maxwell方程组，即电磁场的运动方程。

能不能构造出一个Lagrange量，利用最小作用量就可以得到Schrödinger方程呢？答案是肯定的。

(无相互作用)Schrödinger场的Lagrange量密度可以写为
$
  cal(L) = psi^* i hbar pdv(, t) psi - (hbar^2)/(2m) nabla psi^* dot nabla psi - psi^* U psi
$
其中，$psi = psi(vb(x), t)$。求作用量$S$还要对时空积分，
$
  S = integral dd(t) integral dd(vb(x), 3) cal(L)
$
因此$cal(L)$也可以写为
$
  cal(L) = psi^*(vb(x), t) (i hbar pdv(, t) + (hbar^2)/(2m) laplacian - U(vb(x))) psi(vb(x), t)
$
由于$psi$是复数场，须*将$psi$和$psi^*$看成独立的两个场*。利用变分原理，从$δ S = 0$就可以得到Schrödinger方程式。

按照从Lagrange形式到Hamilton形式的标准步骤，将$psi$和$psi^*$看成“*正则坐标*”，与之对应的“*正则共轭动量*”为
$
  pi(vb(x), t) = pdv(cal(L), dot(psi)(vb(x), t)) = i hbar psi^*(vb(x), t)\
  pi^*(vb(x), t) = pdv(cal(L), dot(psi)^*(vb(x), t)) = 0
$
可以理解为，每一个空间点$vb(x)$的场$psi$和$psi^*$都是体系的“正则坐标”。对于Schrödinger场，*与$psi^*$对应的正则共轭动量场恒为零*。所以，求得Hamilton量为
$
  H = integral dd(vb(x)) (pi dot(psi) - cal(L)) = integral dd(vb(x), 3) psi^*(vb(x), t) (- (hbar^2)/(2m) nabla^2 + U(vb(x))) psi(vb(x), t)
$
为了便于理解，可将空间点离散化（$integral dd(x) -> sum_i$）。

计算等时的“经典”Poisson括号可得到
$
  { psi(vb(x), t), psi(vb(x)', t) } = 0,\
  { psi^*(vb(x), t), psi^*(vb(x)', t) } = 0,\
  { psi(vb(x), t), psi^*(vb(x)', t) } = - i/hbar delta(vb(x) - vb(x)')
$
接下来就是标准的*正则量子化*步骤。将经典场替换为*场算符*
$
  psi(vb(x), t) -> hat(psi)(vb(x), t), pi(vb(x), t) -> hat(pi)(vb(x), t)
$
将Poisson括号替换为*对易关系*，根据以下对应规则
$
  {A, B}_"PB" = 1/(i hbar) [hat(A), hat(B)]
$
得到场算符满足的等时对易关系
$
  [hat(psi) (vb(x), t), hat(psi) (vb(x)', t)] = 0,\
  [hat(psi)^dagger (vb(x), t), hat(psi)^dagger (vb(x)', t)] = 0,\
  [hat(psi) (vb(x), t), hat(psi)^dagger (vb(x)', t)] = delta(vb(x) - vb(x)')
$
#newpara()
将场算符用一组正交完备的基函数${phi_i (vb(x))}$进行展开
$
  hat(psi) (vb(x), t) = sum_i hat(a)_i (t) phi_i (vb(x))\
  hat(psi)^dagger (vb(x), t) = sum_i hat(a)^dagger_i (t) phi_i^* (vb(x))
$
可以证明如下等时对易关系
$
  [hat(a)_i (t), hat(a)_j (t)] = 0,\
  [hat(a)^dagger_i (t), hat(a)^dagger_j (t)] = 0,\
  [hat(a)_i (t), hat(a)^dagger_j (t)] = delta_(i, j)
$
若${phi_i (vb(x))}$就是单粒子Hamilton量的本征态，本征值为${E_i}$，令
$
  hat(a)_i (t) = hat(a)_i e^(- i E_i t/hbar), hat(a)^dagger_i (t) = hat(a)^dagger_i e^(i E_i t/hbar)
$
则得到
$
  [hat(a)_i, hat(a)_j] = 0,
  [hat(a)^dagger_i, hat(a)^dagger_j] = 0,
  [hat(a)_i, hat(a)^dagger_j] = delta_(i, j)
$
以上讨论仅针对Bose子。对于Fermi子，要采用反对易关系。

正则量子化手续完成后，经典场的Hamilton量就成了
$
  hat(cal(H)) = integral dd(vb(x), 3) hat(psi)^dagger (vb(x), t) (- (hbar^2)/(2m) nabla^2 + U(vb(x))) hat(psi) (vb(x), t)
$
即无相互作用多粒子系统的哈密顿量。再进行表象变换，即可得到一般单粒子表象下的形式。

对于有相互作用的情形，拉格朗日量密度则取为
$
  cal(L) = psi^* i hbar pdv(, t) psi - (hbar^2)/(2m) nabla psi^* dot nabla psi - psi^* U psi \ - 1/2 integral dd(vb(x)', 3) psi^* (vb(x), t) psi^* (vb(x)', t) V(vb(x), vb(x)') psi (vb(x)', t) psi (vb(x), t)
$
采用同样的手续，可得到量子化的哈密顿量
$
  hat(cal(H))_"I" = 1/2 integral dd(vb(x), 3) integral dd(vb(x)', 3) hat(psi)^dagger (vb(x), t) hat(psi)^dagger (vb(x)', t) V(vb(x), vb(x)') hat(psi) (vb(x)', t) hat(psi) (vb(x), t)
$
#newpara()

综上，二次量子化有两种理解：
- 二次量子化是多粒子系统的量子化，把单粒子波函数看成多粒子系统的态矢，在Fock空间中引入产生湮灭算符，从而得到多粒子系统的量子力学描述。
- 二次量子化是场的量子化，把Schrödinger方程当作经典场，考虑其Lagrange量，以及Hamilton量，通过正则量子化手续将经典场量子化，

