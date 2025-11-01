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
  ket(alpha)_1 times.circle ket(psi)_2
$
在不引起混淆的情况下，直积符号$⊗$和代表空间的下标都可以省去。例如，两个自旋$1/2$粒子，自旋都向上的状态写为
$
  ket(arrow.t)_1 times.circle ket(arrow.t)_2 = ket(arrow.t)_1 ket(arrow.t)_2 = ket(arrow.t) ket(arrow.t) = ket(arrow.t arrow.t)
$
- *加法*：复合系统的Hilbert空间就是直积空间$cal(R)_1 ⊗ cal(R)_2$。两个直积态相加
  $
    ket(alpha)_1 times.circle ket(psi)_2 + ket(beta)_1 times.circle ket(phi)_2
  $
  仍然是直积空间$cal(R)_1 ⊗ cal(R)_2$中的态矢量，但是一般不能写成*直积态*或者*可分离态*的形式，通常是一种*纠缠态*或者*不可分离态*。实际上，直积空间中的绝大部分量子态都是纠缠态，或者说纠缠态在直积空间中是稠密的。直积满足分配律：
  $
    ket(alpha)_1 times.circle (ket(psi)_2 + ket(phi)_2) = ket(alpha)_1 times.circle ket(psi)_2 + ket(alpha)_1 times.circle ket(phi)_2
  $
- *数乘*：若$c$是一个复数，则
  $
    c (ket(alpha)_1 times.circle ket(psi)_2) = (c ket(alpha)_1) times.circle ket(psi)_2 = ket(alpha)_1 times.circle (c ket(psi)_2)
  $
- *内积*：Hilbert空间$cal(R)_1$中的内积$braket(α, β)$和Hilbert空间$cal(R)_2$中的内积$braket(ϕ, ψ)$都已经定义好了，那么直积空间$cal(R)_1 ⊗ cal(R)_2$中的态矢量$ket(α)_1 times.circle ket(ϕ)_2$和$ket(β)_1 times.circle ket(ψ)_2$的内积定义为
  $
    (bra(alpha) times.circle bra(phi)) (ket(beta) times.circle ket(psi)) = braket(alpha, beta) braket(phi, psi)
  $
  如果两个态矢量不是简单的直积态，那么
  $
    &(bra(alpha_1) times.circle bra(phi_1) + bra(alpha_2) times.circle bra(phi_2)) (ket(beta_1) times.circle ket(psi_1) + ket(beta_2) times.circle ket(psi_2))\
    = & (bra(alpha_1) times.circle bra(phi_1)) (ket(beta_1) times.circle ket(psi_1)) + (bra(alpha_1) times.circle bra(phi_1)) (ket(beta_2) times.circle ket(psi_2)) \ &+ (bra(alpha_2) times.circle bra(phi_2)) (ket(beta_1) times.circle ket(psi_1)) + (bra(alpha_2) times.circle bra(phi_2)) (ket(beta_2) times.circle ket(psi_2)) \
    = & braket(alpha_1, beta_1) braket(phi_1, psi_1) + braket(alpha_1, beta_2) braket(phi_1, psi_2) + braket(alpha_2, beta_1) braket(phi_2, psi_1) + braket(alpha_2, beta_2) braket(phi_2, psi_2)
  $

=== 直积空间中的算符

Hilbert空间$cal(R)_1$中的算符为$hat(A), hat(B), ...$，Hilbert空间$cal(R)_2$中的算符为$hat(M), hat(N), ...$，这些算符对各自Hilbert空间中的态矢量的作用是明确的。定义直积空间中的算符$hat(A) times.circle hat(M)$，它对直积空间中的态矢量的作用为
$
    & (hat(A) times.circle hat(M)) (ket(alpha) times.circle ket(psi) + ket(beta) times.circle ket(phi) + ...) \
  = & hat(A) ket(alpha) times.circle hat(M) ket(psi) + hat(A) ket(beta) times.circle hat(M) ket(phi) + ...
$
算符运算有如下关系：
$
  (hat(A) + hat(B)) times.circle hat(M) = hat(A) times.circle hat(M) + hat(B) times.circle hat(M)\
  (hat(A) times.circle hat(M))(hat(B) times.circle hat(N)) = (hat(A) hat(B)) times.circle (hat(M) hat(N))\
$
上述定义实际上不需要去记忆，只需要记住： *态矢量只和自己空间中的态矢量做内积，算符只作用于自己空间中的态矢量*。在不引起混淆的情况下，直积符号$⊗$通常可以略去。

通常为了方便，我们在直积空间中也说算符$hat(A)$或者算符$hat(M)$，但这并不是指$cal(R)_1$或$cal(R)_2$中的算符，实际上，$hat(A)$是$hat(A)⊗hat(I)_2$的简写，$hat(M)$是$hat(I)_1⊗hat(M)$的简写，其中$hat(I)_1$和$hat(I)_2$分别是$cal(R)_1$和$cal(R)_2$中的单位算符。如果将两个空间的算符直接相加，例如$hat(A) + hat(M)$，实际上是
$
  hat(A) times.circle hat(I)_2 + hat(I)_1 times.circle hat(M)
$
例如，我们将两个电子的总自旋角动量写为
$
  hat(S) = hat(S)_1 + hat(S)_2 = (hat(S)_1 times.circle hat(I)_2 + hat(I)_1 times.circle hat(S)_2)
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
  ket(alpha) times.circle ket(psi) = sum_i sum_j alpha_i psi_j (ket(k_i) times.circle ket(f_j))
$
可见，构造直积空间$cal(R)_1 ⊗ cal(R)_2$的一组完备基矢，最朴素的做法就是将$cal(R)_1$和$cal(R)_2$的基矢两两直积起来，形成一组基矢，就可以叠加出直积空间中所有的态矢量。如果两个空间的维数分别为$n_1$和$n_2$，那么直积空间的维数为$n_1n_2$。

构造了直积空间$cal(R)_1 ⊗ cal(R)_2$的一组完备基矢后，就可以计算在这个表象下，量子态和力学量的矩阵形式 (如果都是离散表象) 。例如，在上述的 KF 表象中，将基矢$ket(k_i) ⊗ ket(f_j)$的排列顺序约定好以后，就可以计算量子态的列向量形式和力学量的矩阵元。例如
$
    & braket(k_i f_j, hat(A) times.circle hat(M), k_m f_n) \
  = & (bra(k_i) times.circle bra(f_j)) (hat(A) times.circle hat(M)) (ket(k_m) times.circle ket(f_n)) \
  = & braket(k_i, hat(A), k_m) braket(f_j, hat(M), f_n) \
  = & A_(i j) M_(m n)
$

#example[
  若$n_1 = 2， n_2 = 3$，那么在 KF 表象中量子态$ket(alpha) ⊗ ket(psi)$的矩阵形式为
  形式为
  $
    ket(alpha) times.circle ket(psi) = mat(alpha_1; alpha_2) times.circle mat(psi_1; psi_2; psi_3) = mat(alpha_1 psi_1; alpha_1 psi_2; alpha_1 psi_3; alpha_2 psi_1; alpha_2 psi_2; alpha_2 psi_3)
  $
  注意：第二个等号后面的形式并不唯一，依赖于基矢排列的顺序。

  力学量$hat(A) times.circle hat(M)$在 KF 表象中的矩阵形式为
  $
    hat(A) times.circle hat(M) = mat(A_(11), A_(12); A_(21), A_(22)) times.circle mat(M_(11), M_(12), M_(13); M_(21), M_(22), M_(23); M_(31), M_(32), M_(33)) \
  $
  这实际上就是矩阵的直积。上面的两个式子可以简写为
  $
    ket(alpha) times.circle ket(psi) = mat(alpha_1 psi; alpha_2 psi)\
    hat(A) times.circle hat(M) = mat(A_(11) hat(M), A_(12) hat(M); A_(21) hat(M), A_(22) hat(M))
  $
]

=== 两粒子自旋态

考虑两个自旋$1/2$的粒子，每个粒子都有一个二维的自旋空间，两个自旋空间是独立的。那么两粒子的总自旋空间如何构造呢？这个操作就是前面提到的直积。

两个二维线性空间的直积是*四维线性空间*。首先构造这个四维空间的基矢，最简单直接的做法是*把原来两个二维空间的基矢两两直积起来*。如果粒子 1 的自旋空间基矢为$ket(arrow.t)_1$和$ket(arrow.b)_1$，粒子 2 的自旋空间基矢为$ket(arrow.t)_2$和$ket(arrow.b)_2$，那么四维直积空间的基矢的最简单构造方法就是
$
  ket(arrow.t arrow.t) = ket(arrow.t)_1 times.circle ket(arrow.t)_2\
  ket(arrow.t arrow.b) = ket(arrow.t)_1 times.circle ket(arrow.b)_2\
  ket(arrow.b arrow.t) = ket(arrow.b)_1 times.circle ket(arrow.t)_2\
  ket(arrow.b arrow.b) = ket(arrow.b)_1 times.circle ket(arrow.b)_2
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

#example()[
  计算$σ_(1x) σ_(2x)$对Bell基矢$ket(psi_plus.minus)$的作用。
]
#solution[
  $
    σ_(1x) σ_(2x) ket(psi_plus.minus) &= σ_(1x) σ_(2x) 1/sqrt(2) (ket(arrow.t)_1 times.circle ket(arrow.b)_2 plus.minus ket(arrow.b)_1 times.circle ket(arrow.t)_2) \
    &= 1/sqrt(2) (σ_(1x) ket(arrow.t)_1 times.circle σ_(2x) ket(arrow.b)_2 plus.minus σ_(1x) ket(arrow.b)_1 times.circle σ_(2x) ket(arrow.t)_2) \
    &= 1/sqrt(2) (ket(arrow.b)_1 times.circle ket(arrow.t)_2 plus.minus ket(arrow.t)_1 times.circle ket(arrow.b)_2) \
    &= plus.minus 1/sqrt(2) (ket(arrow.t arrow.b) plus.minus ket(arrow.b arrow.t)) \
    &= plus.minus ket(psi_plus.minus)
  $
  所以，$ket(psi_plus.minus)$ 是 $σ_(1x) σ_(2x)$ 的本征态。
]

#example()[
  计算在Bell不等式中遇到的关联函数$P(vb(a), vb(b))$的量子力学结果。考虑量子态 $ket(psi_-)$ 即自旋单态，假设 $vb(a)$ 和 $vb(b)$ 是空间任意两个方向上的单位矢量。关联函数可以表示为
  $
    P(vb(a), vb(b)) = braket(psi_-, (vb(sigma)_1 dot vb(a)) times.circle (vb(sigma)_2 dot vb(b)), psi_-)
  $
  即算符$(vb(sigma)_1 dot vb(a)) times.circle (vb(sigma)_2 dot vb(b))$在自旋单态中的平均值。这样计算出来的结果是关联函数的量子力学结果。

  由于计算中公式太长，仅写出计算的概要。将算符$(vb(sigma)_1 dot vb(a)) times.circle (vb(sigma)_2 dot vb(b))$写为明显形式
  $
    (vb(sigma)_1 dot vb(a)) times.circle (vb(sigma)_2 dot vb(b)) = (sigma_(1 x) a_x plus sigma_(1 y) a_y plus sigma_(1 z) a_z) times.circle (sigma_(2 x) b_x plus sigma_(2 y) b_y plus sigma_(2 z) b_z)
  $
  展开后一共 9 项，每项对量子态$ket(psi_-)$作用后的结果都可以按前面演示的方式耐心计算出来。其中的第一项是
  $
    (sigma_(1 x) times.circle sigma_(2 x)) ket(psi_-) = - ket(psi_-)\
    (sigma_(1 y) times.circle sigma_(2 y)) ket(psi_-) = - ket(psi_-)\
    (sigma_(1 z) times.circle sigma_(2 z)) ket(psi_-) = - ket(psi_-)
  $
  其它 6 项，只需要按照定义耐心计算。由于正交性，结果都为 0。最终得到结果为
  $
    P(vb(a), vb(b)) = - (a_x b_x plus a_y b_y plus a_z b_z) = - vb(a) dot vb(b)
  $
]
