#import "@preview/scripst:1.1.1": *

= 量子测量与量子力学的诠释问题

#note[
  本章节的内容颇具争议性，仅作为了解。
]

量子力学取得了巨大的成功，尤其是从Schrödinger方程和Born几率诠释演绎出的结果，经受了无数的实验检验。但是量子力学本身还存在一些尖锐的问题。

- 量子力学对微观世界的描述为什么是*概率性*的？是微观世界的固有属性，还是因为存在某些我们不知道的*隐变量*？对量子系统的物理量进行测量之前，这些物理量是否以确定的数值客观存在？
- 对量子系统的力学量进行测量之后，为什么系统的状态立刻*塌缩*到了力学量的某个本征态？这个塌缩的过程是怎样的？Copenhagen诠释将塌缩当成一个基本假设，认为测量仪器必须是经典的，外部的经典世界是诠释量子力学所必需的。
- 大量的量子单元构成的*宏观*系统，其量子效应是如何消失的？为什么我们没有看到过既死又活的Schrödinger的猫？从微观到宏观，从量子到经典，是怎么过渡的？

== 量子态和算符的直积

我们经常遇到这样的问题：由$N$个$(N ≤ 2)$量子系统构成的*量子复合系统*，例如$N$个电子构成的系统，也可能是微观粒子具有*内部自由度*，例如电子的自旋。这就需要由$N$个已知的Hilbert空间来构造复合系统的Hilbert空间，构造的方法就是*直积*。下面以$N = 2$为例进行讨论，推广到$N > 2$的情形是直接的。

=== 直积空间中的态矢量

假设Hilbert空间$cal(R)_1$中的态矢量为$ket(alpha)_1, ket(beta)_1, ...$，Hilbert空间$cal(R)_2$中的态矢量为$ket(psi)_2, ket(phi)_2, ...$。由$cal(R)_1$中的态矢量$ket(alpha)_1$和$cal(R)_2$中的态矢量$ket(psi)_2$构成的*直积态*写成
$
  ket(alpha)_1 times.o ket(psi)_2
$
在不引起混淆的情况下，直积符号$⊗$和代表空间的下标都可以省去。例如，两个自旋$1/2$粒子，自旋都向上的状态写为
$
  ket(arrow.t)_1 times.o ket(arrow.t)_2 = ket(arrow.t)_1 ket(arrow.t)_2 = ket(arrow.t) ket(arrow.t) = ket(arrow.t arrow.t)
$
- *加法*：复合系统的Hilbert空间就是直积空间$cal(R)_1 ⊗ cal(R)_2$。两个直积态相加
  $
    ket(alpha)_1 times.o ket(psi)_2 + ket(beta)_1 times.o ket(phi)_2
  $
  仍然是直积空间$cal(R)_1 ⊗ cal(R)_2$中的态矢量，但是一般不能写成*直积态*或者*可分离态*的形式，通常是一种*纠缠态*或者*不可分离态*。实际上，直积空间中的绝大部分量子态都是纠缠态，或者说纠缠态在直积空间中是稠密的。直积满足分配律：
  $
    ket(alpha)_1 times.o (ket(psi)_2 + ket(phi)_2) = ket(alpha)_1 times.o ket(psi)_2 + ket(alpha)_1 times.o ket(phi)_2
  $
- *数乘*：若$c$是一个复数，则
  $
    c (ket(alpha)_1 times.o ket(psi)_2) = (c ket(alpha)_1) times.o ket(psi)_2 = ket(alpha)_1 times.o (c ket(psi)_2)
  $
- *内积*：Hilbert空间$cal(R)_1$中的内积$braket(α, β)$和Hilbert空间$cal(R)_2$中的内积$braket(ϕ, ψ)$都已经定义好了，那么直积空间$cal(R)_1 ⊗ cal(R)_2$中的态矢量$ket(α)_1 times.o ket(ϕ)_2$和$ket(β)_1 times.o ket(ψ)_2$的内积定义为
  $
    (bra(alpha) times.o bra(phi)) (ket(beta) times.o ket(psi)) = braket(alpha, beta) braket(phi, psi)
  $
  如果两个态矢量不是简单的直积态，那么
  $
    &(bra(alpha_1) times.o bra(phi_1) + bra(alpha_2) times.o bra(phi_2)) (ket(beta_1) times.o ket(psi_1) + ket(beta_2) times.o ket(psi_2))\
    = & (bra(alpha_1) times.o bra(phi_1)) (ket(beta_1) times.o ket(psi_1)) + (bra(alpha_1) times.o bra(phi_1)) (ket(beta_2) times.o ket(psi_2)) \ &+ (bra(alpha_2) times.o bra(phi_2)) (ket(beta_1) times.o ket(psi_1)) + (bra(alpha_2) times.o bra(phi_2)) (ket(beta_2) times.o ket(psi_2)) \
    = & braket(alpha_1, beta_1) braket(phi_1, psi_1) + braket(alpha_1, beta_2) braket(phi_1, psi_2) + braket(alpha_2, beta_1) braket(phi_2, psi_1) + braket(alpha_2, beta_2) braket(phi_2, psi_2)
  $

=== 直积空间中的算符

Hilbert空间$cal(R)_1$中的算符为$hat(A), hat(B), ...$，Hilbert空间$cal(R)_2$中的算符为$hat(M), hat(N), ...$，这些算符对各自Hilbert空间中的态矢量的作用是明确的。定义直积空间中的算符$hat(A) times.o hat(M)$，它对直积空间中的态矢量的作用为
$
    & (hat(A) times.o hat(M)) (ket(alpha) times.o ket(psi) + ket(beta) times.o ket(phi) + ...) \
  = & hat(A) ket(alpha) times.o hat(M) ket(psi) + hat(A) ket(beta) times.o hat(M) ket(phi) + ...
$
算符运算有如下关系：
$
  (hat(A) + hat(B)) times.o hat(M) = hat(A) times.o hat(M) + hat(B) times.o hat(M)\
  (hat(A) times.o hat(M))(hat(B) times.o hat(N)) = (hat(A) hat(B)) times.o (hat(M) hat(N))\
$
上述定义实际上不需要去记忆，只需要记住： *态矢量只和自己空间中的态矢量做内积，算符只作用于自己空间中的态矢量*。在不引起混淆的情况下，直积符号$⊗$通常可以略去。

通常为了方便，我们在直积空间中也说算符$hat(A)$或者算符$hat(M)$，但这并不是指$cal(R)_1$或$cal(R)_2$中的算符，实际上，$hat(A)$是$hat(A)⊗hat(I)_2$的简写，$hat(M)$是$hat(I)_1⊗hat(M)$的简写，其中$hat(I)_1$和$hat(I)_2$分别是$cal(R)_1$和$cal(R)_2$中的单位算符。如果将两个空间的算符直接相加，例如$hat(A) + hat(M)$，实际上是
$
  hat(A) times.o hat(I)_2 + hat(I)_1 times.o hat(M)
$
例如，我们将两个电子的总自旋角动量写为
$
  hat(S) = hat(S)_1 + hat(S)_2 = (hat(S)_1 times.o hat(I)_2 + hat(I)_1 times.o hat(S)_2)
$

=== 直积空间的基矢

直积空间$cal(R)_1 ⊗ cal(R)_2$中存在一组完备的基矢，直积空间中任意的态矢量都可以用这组基矢展开。假设在$cal(R)_1$中取力学量$hat(K)$对应的$K$表象，正交归一的基矢是${ket(k_i)}$，在$cal(R)_2$中取力学量$hat(F)$对应的$F$表象，正交归一的基矢是${ket(f_j)}$，即
$
  hat(K) ket(k_i) = k_i ket(k_i), hat(F) ket(f_j) = f_j ket(f_j)
$
则$cal(R)_1$中的任意态矢量$ket(alpha)$和$cal(R)_2$中的任意态矢量$ket(psi)$可以展开为
$
  ket(alpha) = sum_i alpha_i ket(k_i), ket(psi) = sum_j psi_j ket(f_j)
$
直积运算得到
$
  ket(alpha) times.o ket(psi) = sum_i sum_j alpha_i psi_j (ket(k_i) times.o ket(f_j))
$
可见，构造直积空间$cal(R)_1 ⊗ cal(R)_2$的一组完备基矢，最朴素的做法就是将$cal(R)_1$和$cal(R)_2$的基矢两两直积起来，形成一组基矢，就可以叠加出直积空间中所有的态矢量。如果两个空间的维数分别为$n_1$和$n_2$，那么直积空间的维数为$n_1n_2$。

构造了直积空间$cal(R)_1 ⊗ cal(R)_2$的一组完备基矢后，就可以计算在这个表象下，量子态和力学量的矩阵形式 (如果都是离散表象) 。例如，在上述的 KF 表象中，将基矢$ket(k_i) ⊗ ket(f_j)$的排列顺序约定好以后，就可以计算量子态的列向量形式和力学量的矩阵元。例如
$
    & braket(k_i f_j, hat(A) times.o hat(M), k_m f_n) \
  = & (bra(k_i) times.o bra(f_j)) (hat(A) times.o hat(M)) (ket(k_m) times.o ket(f_n)) \
  = & braket(k_i, hat(A), k_m) braket(f_j, hat(M), f_n) \
  = & A_(i j) M_(m n)
$

#example[
  若$n_1 = 2， n_2 = 3$，那么在 KF 表象中量子态$ket(alpha) ⊗ ket(psi)$的矩阵形式为
  形式为
  $
    ket(alpha) times.o ket(psi) = mat(alpha_1; alpha_2) times.o mat(psi_1; psi_2; psi_3) = mat(alpha_1 psi_1; alpha_1 psi_2; alpha_1 psi_3; alpha_2 psi_1; alpha_2 psi_2; alpha_2 psi_3)
  $
  注意：第二个等号后面的形式并不唯一，依赖于基矢排列的顺序。

  力学量$hat(A) times.o hat(M)$在 KF 表象中的矩阵形式为
  $
    hat(A) times.o hat(M) = mat(A_(11), A_(12); A_(21), A_(22)) times.o mat(M_(11), M_(12), M_(13); M_(21), M_(22), M_(23); M_(31), M_(32), M_(33)) \
  $
  这实际上就是矩阵的直积。上面的两个式子可以简写为
  $
    ket(alpha) times.o ket(psi) = mat(alpha_1 psi; alpha_2 psi)\
    hat(A) times.o hat(M) = mat(A_(11) hat(M), A_(12) hat(M); A_(21) hat(M), A_(22) hat(M))
  $
]

=== 两粒子自旋态

考虑两个自旋$1/2$的粒子，每个粒子都有一个二维的自旋空间，两个自旋空间是独立的。那么两粒子的总自旋空间如何构造呢？这个操作就是前面提到的直积。

两个二维线性空间的直积是*四维线性空间*。首先构造这个四维空间的基矢，最简单直接的做法是*把原来两个二维空间的基矢两两直积起来*。如果粒子 1 的自旋空间基矢为$ket(arrow.t)_1$和$ket(arrow.b)_1$，粒子 2 的自旋空间基矢为$ket(arrow.t)_2$和$ket(arrow.b)_2$，那么四维直积空间的基矢的最简单构造方法就是
$
  ket(arrow.t arrow.t) = ket(arrow.t)_1 times.o ket(arrow.t)_2\
  ket(arrow.t arrow.b) = ket(arrow.t)_1 times.o ket(arrow.b)_2\
  ket(arrow.b arrow.t) = ket(arrow.b)_1 times.o ket(arrow.t)_2\
  ket(arrow.b arrow.b) = ket(arrow.b)_1 times.o ket(arrow.b)_2
$

#newpara()

*自旋单态和三重态*：基矢的其它取法就是将上述四个基矢做一个*幺正变换*。比如所谓的*耦合表象*，基矢取为
$
  ket(1\,1) = ket(arrow.t arrow.t)\
  ket(1\,0) = 1/sqrt(2) (ket(arrow.t arrow.b) + ket(arrow.b arrow.t))\
  ket(1\,-1) = ket(arrow.b arrow.b)\
  ket(0\,0) = 1/sqrt(2) (ket(arrow.t arrow.b) - ket(arrow.b arrow.t))
$
可以统一表示为$ket(S\, M_S)$，它们是两粒子总自旋平方$hat(S)^2$和总自旋$z$分量$hat(S)_z$的*共同本征态*。此即所谓的*两个自旋角动量的合成*。其中 $ket(1\,1)$ 和 $ket(1\,−1)$ 是两个粒子自旋态的简单直积，是直积态或者可分离态； $ket(1\,0)$ 和 $ket(0\,0)$ 是两个直积态的线性组合，不能写成两个自旋态的简单直积，这是*纠缠态*或者不可分离态。

*Bell基矢*：进一步地，我们可以将$ket(1\,1)$ 和 $ket(1\,-1)$再线性组合一下，构造出所谓的*Bell基矢*：
$
  ket(psi_+) = 1/sqrt(2) (ket(arrow.t arrow.b) + ket(arrow.b arrow.t))\
  ket(psi_-) = 1/sqrt(2) (ket(arrow.t arrow.b) - ket(arrow.b arrow.t))\
  ket(phi_+) = 1/sqrt(2) (ket(arrow.t arrow.t) + ket(arrow.b arrow.b))\
  ket(phi_-) = 1/sqrt(2) (ket(arrow.t arrow.t) - ket(arrow.b arrow.b))
$
这样一套基矢是${hat(S)_(1x) hat(S)_(2x), hat(S)_(1y) hat(S)_(2y), hat(S)_(1z) hat(S)_(2z)}$中任意两个算符的共同本征态。可以看到，每个Bell基矢都是*纠缠态*。

#example(subname: [Pauli矩阵对Bell基矢的作用])[
  计算$σ_(1x) σ_(2x)$对Bell基矢$ket(psi_plus.minus)$的作用。
]
#solution[
  $
    σ_(1x) σ_(2x) ket(psi_plus.minus) &= σ_(1x) σ_(2x) 1/sqrt(2) (ket(arrow.t)_1 times.o ket(arrow.b)_2 plus.minus ket(arrow.b)_1 times.o ket(arrow.t)_2) \
    &= 1/sqrt(2) (σ_(1x) ket(arrow.t)_1 times.o σ_(2x) ket(arrow.b)_2 plus.minus σ_(1x) ket(arrow.b)_1 times.o σ_(2x) ket(arrow.t)_2) \
    &= 1/sqrt(2) (ket(arrow.b)_1 times.o ket(arrow.t)_2 plus.minus ket(arrow.t)_1 times.o ket(arrow.b)_2) \
    &= plus.minus 1/sqrt(2) (ket(arrow.t arrow.b) plus.minus ket(arrow.b arrow.t)) \
    &= plus.minus ket(psi_plus.minus)
  $
  所以，$ket(psi_plus.minus)$ 是 $σ_(1x) σ_(2x)$ 的本征态。
]

#example(subname: [关联函数])[
  计算在Bell不等式中遇到的关联函数$P(vb(a), vb(b))$的量子力学结果。考虑量子态 $ket(psi_-)$ 即自旋单态，假设 $vb(a)$ 和 $vb(b)$ 是空间任意两个方向上的单位矢量。关联函数可以表示为
  $
    P(vb(a), vb(b)) = braket(psi_-, (vb(sigma)_1 dot vb(a)) times.o (vb(sigma)_2 dot vb(b)), psi_-)
  $
  即算符$(vb(sigma)_1 dot vb(a)) times.o (vb(sigma)_2 dot vb(b))$在自旋单态中的平均值。这样计算出来的结果是关联函数的量子力学结果。

  由于计算中公式太长，仅写出计算的概要。将算符$(vb(sigma)_1 dot vb(a)) times.o (vb(sigma)_2 dot vb(b))$写为明显形式
  $
    (vb(sigma)_1 dot vb(a)) times.o (vb(sigma)_2 dot vb(b)) = (sigma_(1 x) a_x plus sigma_(1 y) a_y plus sigma_(1 z) a_z) times.o (sigma_(2 x) b_x plus sigma_(2 y) b_y plus sigma_(2 z) b_z)
  $
  展开后一共 9 项，每项对量子态$ket(psi_-)$作用后的结果都可以按前面演示的方式耐心计算出来。其中的第一项是
  $
    (sigma_(1 x) times.o sigma_(2 x)) ket(psi_-) = - ket(psi_-)\
    (sigma_(1 y) times.o sigma_(2 y)) ket(psi_-) = - ket(psi_-)\
    (sigma_(1 z) times.o sigma_(2 z)) ket(psi_-) = - ket(psi_-)
  $
  其它 6 项，只需要按照定义耐心计算。由于正交性，结果都为 0。最终得到结果为
  $
    P(vb(a), vb(b)) = - (a_x b_x plus a_y b_y plus a_z b_z) = - vb(a) dot vb(b)
  $
]

#example(subname: [四粒子自旋态])[
  定义如下的四粒子自旋态
  $
    ket(Psi)_(1234) = ket(phi_+)_(12) times.o ket(phi_+)_(34)
  $
  证
  $
    ket(Psi)_(1234) = & 1/2 (ket(psi_+)_(14) times.o ket(psi_+)_(23) + ket(psi_-)_(14) times.o ket(psi_-)_(23) \
                      & + ket(phi_+)_(14) times.o ket(phi_+)_(23) + ket(phi_-)_(14) times.o ket(phi_-)_(23))
  $
  其中 (两粒子Bell基矢)
  $
    ket(phi_plus.minus)_(i j) = 1/sqrt(2) (ket(arrow.t)_i times.o ket(arrow.t)_j plus.minus ket(arrow.b)_i times.o ket(arrow.b)_j)\
    ket(psi_plus.minus)_(i j) = 1/sqrt(2) (ket(arrow.t)_i times.o ket(arrow.b)_j plus.minus ket(arrow.b)_i times.o ket(arrow.t)_j)
  $
]

#proof[
  将Bell基矢的明显形式
  $
    ket(phi_+)_(12) = 1/sqrt(2) (ket(arrow.t)_1 times.o ket(arrow.t)_2 + ket(arrow.b)_1 times.o ket(arrow.b)_2)\
    ket(phi_+)_(34) = 1/sqrt(2) (ket(arrow.t)_3 times.o ket(arrow.t)_4 + ket(arrow.b)_3 times.o ket(arrow.b)_4)
  $
  带入，我们得到 (为简洁起见略去直积符号$⊗$)
  $
    ket(Psi)_(1234) &= 1/2 (ket(arrow.t)_1 ket(arrow.t)_2 + ket(arrow.b)_1 ket(arrow.b)_2) (ket(arrow.t)_3 ket(arrow.t)_4 + ket(arrow.b)_3 ket(arrow.b)_4) \
    &= 1/2 (ket(arrow.t)_1 ket(arrow.t)_2 ket(arrow.t)_3 ket(arrow.t)_4 + ket(arrow.t)_1 ket(arrow.t)_2 ket(arrow.b)_3 ket(arrow.b)_4 + ket(arrow.b)_1 ket(arrow.b)_2 ket(arrow.t)_3 ket(arrow.t)_4 + ket(arrow.b)_1 ket(arrow.b)_2 ket(arrow.b)_3 ket(arrow.b)_4)
    &
  $
  每一项都是四个粒子自旋态的直积。由于量子态直积与顺序无关，调换顺序后
  $
    ket(Psi)_(1234) &= 1/2 (ket(arrow.t)_1 ket(arrow.t)_4 ket(arrow.t)_2 ket(arrow.t)_3 + ket(arrow.t)_1 ket(arrow.b)_4 ket(arrow.t)_2 ket(arrow.b)_3 \
      & + ket(arrow.b)_1 ket(arrow.t)_4 ket(arrow.b)_2 ket(arrow.t)_3 + ket(arrow.b)_1 ket(arrow.b)_4 ket(arrow.b)_2 ket(arrow.b)_3)
  $
  进一步采用简略写法$ket(arrow.t arrow.t)_(1 4) = ket(arrow.t)_1 ket(arrow.t)_4$等，可以写成
  $
    ket(Psi)_(1234) &= 1/2 (ket(arrow.t arrow.t)_(14) times.o ket(arrow.t arrow.t)_(23) + ket(arrow.t arrow.b)_(14) times.o ket(arrow.t arrow.b)_(23) \
      & + ket(arrow.b arrow.t)_(14) times.o ket(arrow.b arrow.t)_(23) + ket(arrow.b arrow.b)_(14) times.o ket(arrow.b arrow.b)_(23))
  $
  利用Bell基矢的定义，可反解出
  $
    ket(arrow.t arrow.t) & = 1/sqrt(2) (ket(phi_+) + ket(phi_-)) ,
                           ket(arrow.b arrow.b) & = 1/sqrt(2) (ket(phi_+) - ket(phi_-)) \
    ket(arrow.t arrow.b) & = 1/sqrt(2) (ket(psi_+) + ket(psi_-)) ,
                           ket(arrow.b arrow.t) & = 1/sqrt(2) (ket(psi_+) - ket(psi_-))
  $
  配上适当的粒子编号，得到
  $
    ket(Psi)_(1234) = & 1/2 (ket(psi_+)_(14) times.o ket(psi_+)_(23) + ket(psi_-)_(14) times.o ket(psi_-)_(23) \
                      & + ket(phi_+)_(14) times.o ket(phi_+)_(23) + ket(phi_-)_(14) times.o ket(phi_-)_(23))
  $
]

== Copenhagen诠释

在第一章中，我们将量子力学的测量公设表述为：

#definition(subname: [量子力学的测量公设])[
  设量子系统的力学量$hat(A)$的正交归一本征态为${ket(phi_n)}$，对应本征值为${a_n}$。若系统处于叠加态 $ket(j) = sum_n c_n ket(phi_n)$，则测量力学量$hat(A)$只能随机得到某个本征值$a_n$。若 $ket(j)$ 是归一的，那么得到结果$a_n$ 的概率为$P_n = abs(c_n)^2$。
]

这实际上就是*Born概率诠释*。

原始版本的Born概率诠释说，若粒子的波函数为$psi(r, t)$，则$t$时刻在$vb(r)$附近体积元 $dd(v)$内发现粒子的概率为 $abs((r; t))^2 dd(v)$。上面的表述将对位置的测量推广到了任意
力学量。对这条公设的表述在不同的书上可能是不同的。例如，有的书上干脆不提测量，说力学量*取值*为 $a_n$ 的概率为 $P_n = abs(c_n)^2$。也有的书上将这条公设称为*平均值公设*，并将其表述为：力学量$hat(A)$的(量子系综)*期望值*为
$
  expval(A) = sum_n a_n abs(c_n)^2 = braket(j, hat(A), j)
$

#newpara()

只提力学量的取值概率甚至只提力学量的期望值，就是所谓的*Shut up and calculate*。

第二章的量子动力学建立在波函数 (量子态) 公设、力学量公设和Schrödinger方程公设的基础上，量子系统在Schrödinger方程支配下的动力学时间演化是一个*幺正演化*。

如果抛开波函数对微观世界的描述是不是完备的 (上帝掷不掷骰子) 这个问题，那么量子动力学本身是自洽的体系，而且与经典物理中的动力学演化图像一脉相承。当然，波函数并不是可观测物理量。要将波函数的动力学演化与实验观测联系起来，我们通常只需要计算力学量的期望值 (实验系统通常是*量子系综*)。

因此，量子力学割裂成了两块：
+ 舒适区：在Schrödinger方程 + Born概率诠释 (量子系综) 的基础上发展起来的那部分量子理论，体系自洽，取得了巨大的成功
+ 争议区：波函数能否完备地描述微观世界？量子测量怎么定义？它到底是怎样的一个过程？这部分仍然还不清楚！

von Neumann进一步追问：测量以后量子系统的波函数是什么？并要求对紧接着的重复测量给出相同的结果。他先把测量理解为*相互作用产生了仪器 D 和系统 S 的关联或纠缠*，相互作用使总系统 D+S 的波函数从最初的可分离态演化为*量子纠缠态*：
$
  ket(psi(0)) = ket(S) times.o ket(D_0) -> ket(psi(t)) = sum_n c_n ket(phi_n) times.o ket(d_n)
$
根据这个假设，一旦发现仪器处于状态$ket(d_n)$，则整个波函数塌缩到 $ket(phi_n) ⊗ ket(d_n)$，从而由仪器状态$ket(d_n)$读出系统的状态$ket(phi_n)$，也就得到了力学量$hat(A)$的单次测量结果 $a_n$。这就是量子测量的*正交投影模型*或von Neumann模型。

进一步地，人们把这种波包塌缩现象简化为一个不能由Schrödinger方程描述的非幺正过程：在 $ket(S)$ 上测量$hat(A)$，一旦得到结果$a_n$，则测量以后的波函数变为$ket(S)$的一个分支$ket(phi_n)$。这个假设的确保证了紧接着的重复测量给出相同的结果。

Bohr不满足于物理层面上的直接描述和数学表达，他提出：*只有外部经典世界的存在，才能引起波包塌缩这种非幺正变化，即测量仪器必须是经典的。外部经典世界(包括测量仪器和观察者)的存在是诠释量子力学必不可少的。*

如果说前面的正交投影模型只是对量子测量的一个数学模型，那么Bohr则将波包塌缩的“神秘”现象进行了哲学层面的提升，成为*Copenhagen诠释*的核心内容。

Copenhagen诠释是由Bohr和Heisenberg等人在 1925-1927 年间提出的量子力学的一种诠释，至今也是量子力学的正统诠释。Copenhagen诠释的内容有不同的版本，不过都强调了量子力学对微观世界的描述本质上是概率性的，可能还包括不确定性原理、互补原理和对应原理。不过，其中最独特的部分就是对波包塌缩的诠释。

Copenhagen诠释中关于测量的说法导致了量子力学的诠释的二元论结构：微观系统的动力学演化服从Schrödinger方程，是一种幺正演化(U 过程)，但是对其测量过程的描述则*必须借助于外部经典世界*，量子测量是一个非幺正的突变(R 过程)。Copenhagen诠释对测量仪器和测量过程的经典处理显然是令人疑惑的，受到了很多物理学家的质疑。Weinberg在《Einstein的错误》一文中写道：

#quote(block: true)[
  Copenhagen诠释试图描述观测量子系统所发生的状况，经典地处理观察者与测量过程。这种处理方法肯定是*不对*的：观察者与他们的仪器也得遵守同样的量子力学规则，正如宇宙中的每一个量子系统都必须遵守量子力学规则......Copenhagen诠释明显地可以解释量子系统的量子行为，但它并没有达成解释的任务，那就是应用波函数演化确定性方程Schrödinger方程于观察者和他们的仪器。
]

由于仪器本身是由微观组元构成，每个组元服从量子力学，因而经典与量子之间的边界是模糊的。Copenhagen诠释为了保证自洽性，就要根据需要调整边界的位置，可以在仪器-系统之间，也可以在仪器-人类观察者之间，甚至可以在视觉神经和人脑之间。

Copenhagen诠释的波包塌缩假设与狭义相对论有表观上的冲突。例如，一个粒子在 $t = 0$ 时刻局域在某个空间点附近 (称为时空点$A$)，其波函数是一个波包，可以用平面波展开为
$
  psi(vb(r)) = integral dd(vb(p)) phi(vb(p)) e^(i/ hbar vb(p) dot vb(r))
$
在$t = T$时刻测量其动量得到确定的动量$vb(p)$，则波包塌缩为动量本征态$∼ e^(i/ hbar vb(p) dot vb(r))$，其空间分布在$t > T$时刻不再定域。因此，测量引起的波包塌缩导致了*定域性的破坏*：假设时空点 $A$ 和 $B$ 之间是类空的，那么在$t > T$时刻仍有可能在$B$点发现粒子。

== 测量过程的量子理论

测量是科学定量化描述客体的方法。物理学不仅以测量为关键方法，而且物理学理论也必须能够描述测量过程本身。测量的物理本质是通过*相互作用*，用一种*已经标定了的物理系统(测量仪器)对被测量系统进行定量化标定*。然而，测量对宏观经典系统和微观量子系统的影响却有很大区别。
- 对于*宏观经典系统*，测量仪器对被测量系统状态的影响可以努力做到忽略不计，或者通过系统误差分析在确定的精度上扣除测量的影响。
- 对于*微观量子系统*，测量仪器对被测量系统状态的影响不可忽略，甚至可以说是*破坏性的*，即波包塌缩。
如果我们相信，微观系统存在独立于主观意识之外的客观属性，那么也应该希望量子力学能够统一描述量子系统的测量过程，包括被测量系统和测量仪器。这个理论必须回答如下几个问题：什么是量子测量？如何描述量子测量？测量后的效果是什么？

=== 正交投影假设

关于测量后系统的状态是什么，波包塌缩假设已经给出了答案。对于处于叠加态
$
  ket(psi) = sum_n c_n ket(phi_n)
$
上的量子系统，测量力学量$hat(A)$，其中$ket(phi_n)$ 是对应于本征值$a_n$的本征态。一旦在一次测量中得到确定性结果$a_n$，则系统状态变为$ket(phi_n)$，即
$
  ket(psi) = sum_n c_n ket(phi_n) -> ket(phi_n)
$
这个假设称为*正交投影(projection)假设*或者*波函数约化(reduction)假设*。它是一个非幺正过程，不能由Schrödinger方程描述。它的出发点是要保证第一次测量后的后续测量给出相同的结果。

=== Heisenberg退相干解释

Heisenberg退相干解释：对于测量如何引起系统状态的改变，Heisenberg提出了一种模型，认为测量后，仪器会在系统的叠加态的每一个分支中引入随机相位，即测量引起的系统状态变化为
$
  ket(psi) = sum_n c_n ket(phi_n) -> sum_n c_n e^(i theta_n) ket(phi_n)
$
其中$θ_n$是随机相位，其相对相位因子
$
  e^(i (theta_m - theta_n)) = e^(i Delta_(m n))
$
也是随机的。完全随机性要求
$
  expval(e^(i theta_n)) = 0, expval(e^(i Delta_(m n))) = 0
$
注意：上述平均不是量子力学平均，而是随机变量求平均。

状态$ket(psi_f)$对应的密度矩阵为
$
  hat(rho)_f = ketbra(psi_f) = sum_m sum_n abs(c_n)^2 ketbra(phi_n) + sum_(m != n) c_m c_n^* e^(i Delta_(m n)) ketbra(phi_m, phi_n)
$
对随机变量取平均后，非对角元项消失，密度矩阵变为
$
  hat(rho)'_f = sum_n abs(c_n)^2 ketbra(phi_n)
$
Heisenberg的上述唯象描述意味着，测量是一个从相干的量子状态到经典随机态 (混合态) 的转变，即
$
  hat(rho) = ketbra(psi) = sum_n abs(c_n)^2 ketbra(phi_n) + sum_(m != n) c_m c_n^* ketbra(phi_m, phi_n) \ -> hat(rho)'_f = sum_n abs(c_n)^2 ketbra(phi_n)
$
这种密度矩阵的非对角项消失的过程叫做*量子退相干*(quantumdecoherence)。退相干后波函数描述的量子概率变成了*经典概率*，波包塌缩就转化为经典概率的随机实现问题。

前面提到，测量公设或者Born概率诠释有另外一种表述：观测结果是力学量的期望值。可以计算，力学量$hat(A)$在$ket(j)$态下的期望值为
$
  braket(psi, hat(A), psi) = Tr(hat(rho) hat(A)) = sum_n abs(c_n)^2 braket(phi_n, hat(A), phi_n) = sum_n a_n abs(c_n)^2
$
在$hat(rho)'_f$下描述的状态 (混合态) 下的期望值为
$
  expval(hat(A))_f = Tr(hat(rho)'_f hat(A)) = sum_n abs(c_n)^2 braket(phi_n, hat(A), phi_n) = sum_n a_n abs(c_n)^2
$
两个结果相同。所以，Heisenberg退相干假设与波包塌缩假设都能满足重复性的要求，即可重复性要求并没有排除波包塌缩以外的其他假设。

如果把测量过程描述为Heisenberg提出的退相干过程，则可以不借助任何经典理论而只用量子力学来描述测量过程。

=== von Neumann测量模型

von Neumann首先提出，测量仪器也应遵守量子力学规律，他把测量考虑为系统和仪器相互作用的结果。将被测量系统S和测量仪器D看成一个封闭系统，服从量子力学幺正演化。总系统的Hamilton量写为
$
  hat(H) = hat(H)_"S" + hat(H)_"D" + hat(V)
$
其中$hat(H)_"S"$和$hat(H)_"D"$分别是系统S和仪器D的Hamilton量，$hat(V)$是系统S和仪器D之间的相互作用。需要注意的是，上式应理解为
$
  hat(H) = hat(H)_"S" times.o hat(I)_"D" + hat(I)_"S" times.o hat(H)_"D" + hat(V)_"SD"
$
假设相互作用部分不影响系统的演化，即
$
  [hat(H)_"S", hat(V)] = 0
$
测量系统的力学量$hat(A)$，通常要求$hat(A)$是系统的守恒量，即
$
  [hat(H)_"S", hat(A)] = 0
$
设$hat(A)$和$hat(H)_"S"$的共同本征态为$|phi_n⟩$，即
$
  hat(A) ket(phi_n) = a_n ket(phi_n), hat(H)_"S" ket(phi_n) = E_n ket(phi_n)
$
由于$[hat(H)_"S", hat(V)] = 0$
$
  hat(V) ket(phi_n) = hat(V)_n (D) ket(phi_n)
$
这里需要注意的是$hat(V)$已经是定义在总Hilbert空间的算符，作用在$ket(phi_n)$上得到的本征值$hat(V)_n (D)$是仪器Hilbert空间的算符。设$t = 0$时刻系统处于叠加态，仪器处于状态$ket(D)$，因而总系统处于可分离态
$
  ket(psi(0)) = (sum_n c_n ket(phi_n)) times.o ket(D)
$
总系统的演化服从Schrödinger方程，我们得到
$
  ket(psi(t)) & = e^(- i/hbar (hat(H)_"S" + hat(H)_"D" + hat(V)) t) ket(psi(0)) \
              & = e^(- i/hbar (hat(H)_"S" + hat(H)_"D" + hat(V)) t) (sum_n c_n ket(phi_n)) times.o ket(D)
$
由于$[hat(H)_"S", hat(V)] = 0$ 且 $[hat(H)_"S", hat(H)_"D"] = 0$(两个独立的空间)，所以得到
$
  ket(psi(t)) & = e^(- i/hbar hat(H)_"S" t) e^(- i/hbar (hat(H)_"D" + hat(V)) t) (sum_n c_n ket(phi_n)) times.o ket(D) \
              & = sum_n c_n e^(- i/hbar E_n t) ket(phi_n) times.o (e^(- i/hbar (hat(H)_"D" + hat(V)_n (D)) t) ket(D)) \
              & = sum_n c_n ket(phi_n (t)) times.o ket(d_n (t))
$
其中
$
  ket(phi_n (t)) & = e^(- i/hbar E_n t) ket(phi_n) \
    ket(d_n (t)) & = e^(- i/hbar (hat(H)_"D" + hat(V)_n (D)) t) ket(D)
$
可见，在系统-仪器相互作用支配下，*初始的可分离态演化成了一个纠缠态*。如果在$t$时刻仪器的状态是正交的，即
$
  braket(d_m (t), d_n (t)) = delta_(m n)
$
则此时$ket(psi(t))$是一个关于基矢$ket(phi_n)$的Schmidt分解，$c_n$为Schmidt系数，我们称这个特殊的相互作用实现了一个理想的*预测量(premeasurement)*。显然，*预测量过程是一个幺正的演化过程，无须引入随机相位*。

如果观察者只关注系统的状态，不想了解系统以外的部分的状态，就可以对系统以外的自由度求迹，得到系统的*约化密度矩阵*
$
  hat(rho)_"S" (t) = Tr_"D"hat(rho)_"total" = Tr_"D" ketbra(psi(t))
$
计算上式中对仪器Hilbert空间的求迹。总系统密度矩阵为
$
  hat(rho)_"total" & = ketbra(psi(t)) \
                   & = sum_m sum_n c_n ket(phi_n (t)) times.o ket(d_n (t)) c_m^* bra(phi_m (t)) times.o bra(d_m (t)) \
                   & = sum_m sum_n c_n c_m^* ketbra(phi_n (t), phi_m (t)) times.o ketbra(d_n (t), d_m (t))
$
设$ket(α)$是仪器Hilbert空间的一组完备正交基矢，则可以计算
$
  Tr_"D" ketbra(d_n (t), d_m (t)) = sum_α braket(α, d_n (t)) braket(d_m (t), α) = braket(d_m (t), d_n (t))
$
最终得到系统的约化密度矩阵为
$
  hat(rho)_"S" (t) = sum_n abs(c_n)^2 ketbra(phi_n (t)) + sum_(m != n) c_n c_m^* F_(m n) (t) ketbra(phi_n (t), phi_m (t))
$
其中系统外部 (仪器) 波函数的交叠
$
  F_(m n) (t) = braket(d_m (t), d_n (t))
$
称为*退相干因子*。理想的量子测量要求
$
  F_(m n) (t) -> 0
$
即：测量导致量子退相干
$
  hat(rho)_"S" (t) -> sum_n abs(c_n)^2 ketbra(phi_n (t))
$

#newpara()

#example(subname: [Stern-Gerlach实验])[
  上面的理论描述过于抽象，可以用Stern-Gerlach实验来进行具体计算。

  在Stern-Gerlach实验中，“系统”的状态是银原子的自旋态，“仪器”的状态则是银原子的空间波函数。考虑自旋$1/2$、磁矩为$\mu$、质量为$m$的粒子沿着$x$方向通过沿$z$方向的非均匀磁场$B_z = lambda z$，Hamilton量为
  $
    hat(H) = hat(p)_z^2 / (2 m) - mu hat(S)_z B_z sigma_z = hat(p)_z^2 / (2 m) - mu lambda z sigma_z
  $
  粒子在$x, y$方向上的运动是平庸的，为了方便这里不考虑。粒子的初始状态为
  $
    ket(psi(0)) = (alpha ket(arrow.t) + beta ket(arrow.b)) times.o ket(phi(z, 0))
  $
  其中$ket(phi(0))$为粒子的空间量子态。
]
采用坐标表象，初态则写为
$
  braket(z, psi(0)) = (alpha ket(arrow.t) + beta ket(arrow.b)) times.o braket(z, phi(0)) =
  (alpha ket(arrow.t) + beta ket(arrow.b)) phi(z, 0)
$
假设粒子的空间波函数$phi(z, 0)$是一个Gauss波包，即
$
  phi(z, 0) = (1/(2 pi a^2))^(1/4) exp(- z^2 / (4a^2))
$
下面计算粒子状态的时间演化
$
  ket(psi(t)) & = e^(- i/hbar hat(H) t) ket(psi(0)) \
              & = e^(- i/hbar hat(H) t) (alpha ket(arrow.t) + beta ket(arrow.b)) times.o ket(phi(0)) \
$
将Hamilton量带入，计算得到
$
  e^(-i / hbar (hat(p)_z^2 / (2 m) - mu lambda z sigma_z) t) ket(arrow.t) = ket(arrow.t) times.o e^(-i / hbar (hat(p)_z^2 / (2 m) - mu lambda z) t) \
  e^(-i / hbar (hat(p)_z^2 / (2 m) - mu lambda z sigma_z) t) ket(arrow.b) = ket(arrow.b) times.o e^(-i / hbar (hat(p)_z^2 / (2 m) + mu lambda z) t)
$
从而得到
$
  ket(psi(t)) & = e^(- i/hbar (hat(p)_z^2 / (2 m) - mu lambda z sigma_z) t) (alpha ket(arrow.t) + beta ket(arrow.b)) times.o ket(phi(0)) \
  & = alpha ket(arrow.t) times.o e^(-i / hbar (hat(p)_z^2 / (2 m) - mu lambda z) t) ket(phi(0)) + beta ket(arrow.b) times.o e^(-i / hbar (hat(p)_z^2 / (2 m) + mu lambda z) t) ket(phi(0)) \
  & = alpha ket(arrow.t) times.o ket(phi_+(t)) + beta ket(arrow.b) times.o ket(phi_-(t))
$
其中
$
  ket(phi_+(t)) = e^(-i / hbar (hat(p)_z^2 / (2 m) - mu lambda z) t) ket(phi(0)) \
  ket(phi_-(t)) = e^(-i / hbar (hat(p)_z^2 / (2 m) + mu lambda z) t) ket(phi(0))
$
可见，*在非均匀磁场提供的相互作用下，粒子的状态从可分离态演化为量子纠缠态*。进入坐标表象，写为
$
  braket(z, psi(t)) & = alpha ket(arrow.t) times.o braket(z, phi_+(t)) + beta ket(arrow.b) times.o braket(z, phi_-(t)) \
                    & = alpha ket(arrow.t) phi_+(z, t) + beta ket(arrow.b) phi_-(z, t)
$
接下来还需要计算两个矩阵元
$
  phi_plus.minus (z, t) &= braket(z, e^(-i / hbar (hat(p)_z^2 / (2 m) minus.plus mu lambda z) t), phi(0)) \
  &= integral dd(z') braket(z, e^(-i / hbar (hat(p)_z^2 / (2 m) minus.plus mu lambda z) t), z') braket(z', phi(0))\
  &= integral dd(z') K_minus.plus (z, t; z', 0) phi(z', 0)
$
其中$phi(z', 0)$为Gauss波包，形式前面已经给出，而$K_minus.plus (z, t; z', 0)$为粒子在线性势$V(z) = minus.plus mu lambda z$下的传播子，其结果为 (见第四章)
$
  K_minus.plus (z, t; z', 0) = sqrt(m / (2 pi hbar i t)) exp(i / hbar ((m (z - z')^2) / (2 t) minus.plus (mu lambda t)/2 (z + z') - (mu^2 lambda^2 t^2) / (24)))
$
带入，这是一个关于$z'$的Gauss积分，只需耐心计算即可。最终结果为
$
  phi_plus.minus (z, t) = (a^2/(2 pi))^(1/4)/(a^2 + (i tilde(t)) / (2 m))^(1/2) exp(- (i lambda^2 mu^2 tilde(t)^3) / (6 m) - ((z minus.plus lambda mu tilde(t)^2 / (2 m))^2) / (4 sqrt(a^2 + (i tilde(t)) / (2 m))) plus.minus i lambda mu z tilde(t))
$
其中$tilde(t) = t/hbar$。由此可以计算出交叠积分
$
  F(t) = braket(phi_+(t), phi_-(t)) = integral dd(z) phi_-^* (z, t) phi_+ (z, t)
$
最终结果为
$
  abs(F(t)) = exp(- 2 a^2 lambda^2 mu^2 tilde(t)^2 - (lambda^2 mu^2 tilde(t)^4) / (8 m^2 a^2))
$
*总结*：在Stern-Gerlach实验中，$z$方向的非均匀磁场提供了粒子的自旋自由度和空间自由度之间的相互作用，这使得我们可以将粒子的自旋部分视作被测量系统，而空间部分视为测量仪器。在相互作用下，初始的可分离态演化为量子纠缠态
$
  ket(psi(0)) = (alpha ket(arrow.t) + beta ket(arrow.b)) times.o ket(phi(0)) \ -> ket(psi(t)) = alpha ket(arrow.t) times.o ket(phi_+(t)) + beta ket(arrow.b) times.o ket(phi_-(t))
$
经过时间足够长的时间，“仪器”波函数的交叠
$
  braket(phi_+(t), phi_-(t)) -> 0
$
则Stern-Gerlach实验实现了一个*理想的预测量*。最终，从胶片上斑点的位置可以读出粒子自旋状态的测量结果。

=== von Neumann测量模型的问题

仪器-系统的相互作用导致总系统的状态演化为量子纠缠态，从而导致了预测量。其解释是基于von Neumann本人提出的*波包塌缩假设*：在纠缠态
$
  ket(psi(t)) = sum_n c_n ket(phi_n (t)) times.o ket(d_n (t))
$
中，如果观察到仪器处于$ket(d_n)$态，则总体便塌缩到$ket(phi_n) times.o ket(d_n)$，并由此推断出系统处于$ket(phi_n)$态。

但是问题来了：怎么判断仪器处于$ket(d_n)$态？必须有第二个仪器D(2)去测量第一个仪器D(1)=D(例如，Stern-Gerlach实验中，最终还需要由胶片来测量粒子的空间状态)。以此类推，就存在一个仪器链，即*von Neumann链*。那么von Neumann链的末端是什么呢？Wigner等人认为必须是一个*有意识的观察者*，即通过测量得到结果的过程是一个有观察者介入的过程，因此通过包含主观要素的测量确定下来的微观系统的性质，不能独立于人的意识而存在。

上述量子测量理论还存在一个问题： ${ket(phi_n)}$ 只是系统Hilbert空间的一组基矢，也可以选取其他的基矢 ${ket(xi_m)}$，两者之间通过幺正变换联系。因此，总系统演化出的纠缠态
$
  ket(psi(t)) = sum_n c_n ket(phi_n (t)) times.o ket(d_n (t))
$
也可以用另一组基矢${ket(xi_m)}$表达为
$
  ket(psi(t)) = sum_m a_m ket(xi_m (t)) times.o ket(d'_m (t))
$
即：可以根据系统任何给定的基矢进行Schmidt分解，从而相互作用导致的抽象的纠缠态无法直接描述量子测量，因为我们无法预先知道为什么要选基矢$ket(phi_n)$而不是$ket(xi_m)$。这就是所谓的预选基问题或者基矢偏好问题。

为了解决这个问题， Zurek 等人在von Neumann量子测量模型中进一步考虑*环境*的作用。考虑环境E，它也是量子系统，与被测量系统S和测量仪器D之间存在相互作用。E可以是仪器的一部分，它具有大量的自由度。*它的作用是吸收密度矩阵非对角元所包含的信息从而导致波函数塌缩*，从而确定被测量的可观测量(*指针可观测量*)及其本征态(*指针基*)。假设*SDE复合系统*的初态为
$
  ket(psi(0)) = (sum_n c_n (0) ket(phi_n)) times.o ket(D(t=0)) times.o ket(E(t=0))
$
由于S和D的相互作用，演化到某时刻$t_1$，S和D之间产生了纠缠，但此时与环境E的相互作用尚未启动。此时SDE复合系统的状态为
$
  ket(psi(t_1)) = (sum_n c_n (t_1) ket(phi_n) times.o ket(d_n (t_1)) )times.o ket(E(t_1))
$
此时指针观测量尚未确定，SD系统的纠缠态也可以用其他基矢来展开。环境与测量仪器的相互作用从$t_1$时刻开启，到$t_2$时刻$(t_2 > t_1)$，SDE复合系统的状态演化为
$
  ket(psi(t_2)) = sum_n c_n (t_2) ket(phi_n) times.o ket(d_n (t_2)) times.o ket(e_n (t_2))
$
此时，环境与被测量系统和测量仪器建立了关联，从而确定了指针基${ket(phi_n))$。$t > t_2$时D+S系统的约化密度矩阵要对环境求迹：
$
  hat(rho) (t) & = Tr_"E" ketbra(psi(t)) \
  & = sum_(m n) c_n c^*_m ketbra(phi_n, phi_m) times.o ketbra(d_n (t), d_m (t)) times.o Tr_"E" braket(e_m (t), e_n (t))\
  &= sum_(m n) c_n c^*_m braket(e_n (t), e_m (t)) ketbra(phi_n, phi_m) times.o ketbra(d_n (t), d_m (t))
$
如果用${ket(xi_n)}$为基矢重新展开S的态，则
$
  hat(rho) (t) & = sum_(m n) a_n a^*_m braket(e'_n (t), e'_m (t)) ketbra(xi_n, xi_m) times.o ketbra(d'_n (t), d'_m (t))
$
不同基矢展开的区别就在于环境态的交叠积分不同，只有指针基${ket(phi_n)}$可以使得环境态的交叠积分随时间迅速衰减到零，从而实现退相干，解决了预选基问题。

对一些简单的可解模型的研究表明，当$t$足够大时
$
  braket(e_n (t), e_m (t)) -> delta_(m n)
$
通常是以指数形式快速趋近。这样，D+S系统的密度矩阵$hat(ρ)$的非对角元随时间指数衰减，即指针可观测量的本征态之间的量子相干性随时间指数衰减，使被测量系统进入指针可观测量的任何一个本征态而不是它们的叠加，从而完成退相干，
$
  hat(rho) (t) -> sum_n abs(c_n)^2 ketbra(phi_n) times.o ketbra(d_n (t))
$
这样，就可以利用系综的观念解释现有实验的一切问题，不必强调单次测量中的随机塌缩，也无须引入人的主观意识。环境导致的退相干也可以解释Schrödinger的猫这样的佯谬。

#footnote[
  具体模型计算可参考： W. H. Zurek, Physical Review D 24, 1516 (1980); 26, 1862
  (1982)。
]
