#import "@preview/scripst:1.1.1": *

= 对称性理论

*对称性*在物理学中具有极其重要的地位。

- 在经典力学中，对称性变换是对力学量的变换。*Noether定理*告诉我们，如果系统的作用量在某个连续整体对称性变换下不变，则存在对应的*守恒量*，利用守恒量往往可以简化运动方程的求解。在经典物理中，对称性还可能导致物理量之间存在某种关系，从而使问题简化(例如：介质的应力张量、介电张量等)。
- 在量子力学中，由于量子系统的状态用波函数或者Hilbert空间的态矢来描述，对称性有了更加丰富的含义和物理后果。例如，量子系统的能级简并一般就反映了系统的对称性。如何定义*量子系统的对称性变换*？如何定义*量子系统在对称性变换下的不变性*？*量子系统的对称性会带来哪些物理后果*？这是本章要探讨的问题。

== 对称性变换

在经典力学中，对称性变换是对力学量直接变换。在量子力学中，描述系统状态的是Hilbert空间中的态矢，力学量成了算符。我们预期，量子系统的对称性变换包含了对态矢和力学量的变换。

#definition(subname: [量子系统的对称性变换])[
  量子系统的对称性变换，就是在保持系统的*物理观测结果不变*的情况下，把它的力学量和态矢都变换为新的力学量和态矢。具体来说，对称性变换$cal(Q)$将系统的力学量$hat(Ω)$变换为系统的另一个力学量$hat(Omega)'$，
  $
    hat(Omega) ->^(cal(Q)) hat(Omega)'
  $
  同时将系统任意可能的状态$ket(u)$变换为另一个状态$ket(u')$(*整体相位因子未定*)
  $
    ket(u) ->^(cal(Q)) ket(u')
  $
]

什么叫做*物理观测结果不变*？根据量子力学的测量公设，量子系统中可以观测的是系统的力学量，观测结果是*力学量的本征值和取值几率*。所以，物理观测结果不变很自然地要求对称性变换要满足如下两个条件。
+ 任意的力学量$hat(Ω)$和变换后的力学量$hat(Ω)'$具有*相同的本征值谱*。设$hat(Ω)$和$hat(Ω)'$的本征方程为
  $
    hat(Omega) ket(omega_n) = omega_n ket(omega_n), hat(Omega)' ket(omega'_n) = omega'_n ket(omega'_n)
  $
  则要求本征值谱${omega_n}$和${omega'_n}$一致。不失一般性，设$omega'_n = omega_n$，即$ket(omega_n) ->^cal(Q) ket(omega'_n)$。
+ 在任意的状态$ket(psi)$下对任意的力学量$hat(Ω)$进行测量，根据*Born几率诠释*，变换前后*得到任意本征值的几率不变*，即
  $
    abs(braket(n, psi))^2 = abs(braket(n', psi'))^2
  $
  假设$ket(u)$和$ket(v)$是系统的任意两个状态，对称性变换$cal(Q)$将它们变换为$ket(u')$和$ket(v')$，即
  $
    ket(u) ->^(cal(Q)) ket(u'),
    ket(v) ->^(cal(Q)) ket(v')
  $
  那么Born几率诠释要求
  $
    abs(braket(u, v))^2 = abs(braket(u', v'))^2
  $
  为什么可以将上面任意力学量$hat(Ω)$的本征态推广到*任意状态*？实际上，系统的任意状态$ket(psi)$都可以是某个Hermite算符的本征态。当系统处于$ket(v)$态时，测量该Hermite算符对应的力学量，得到$ket(v)$对应的本征值的概率为$abs(braket(v, psi))^2$。

  最简单的例子是投影算符$hat(P)_u = ketbra(u)$，变换后为$hat(P)_u' = ketbra(u')$，变换前后$hat(P)_u$的平均值不变，得到
  $
    bra(v) hat(P)_u ket(v) = bra(v') hat(P)_u' ket(v') => abs(braket(u', v'))^2 = abs(braket(u, v))^2
  $

*满足Born几率诠释要求的变换具有怎样的性质？*Wigner证明了著名的Wigner定理：
#theorem(subname: [Wigner定理])[
  如果系统Hilbert空间中的态矢之间存在一一对应$cal(Q)$，它保持任意两个态矢之间的内积的模方不变，则总可以适当地选择相位因子，使得$cal(Q)$要么是*线性幺正变换*，要么是*反线性幺正变换*。
]
可以具体表述为：假设$ket(u)$和$ket(v)$是系统Hilbert空间中的任意两个态矢，可逆变换对应的算符$hat(Q)$将它们变换为
$
  ket(u') = hat(Q) ket(u), ket(v') = hat(Q) ket(v)
$
如果算符$hat(Q)$保持$ket(u)$和$ket(v)$之间的内积的模方不变，即
$
  abs(braket(u, v))^2 = abs(braket(u', v'))^2
$
那么总可以合理选择态矢的相位因子，使得$hat(Q)$要么是*线性幺正算符*，要么是*反线性幺正算符*，两种情况都有$hat(Q)^dagger hat(Q) = hat(Q) hat(Q)^dagger = 1$。

#note(subname: [线性和反线性幺正算符])[
  - 线性幺正算符$hat(Q)$满足线性条件
    $
      hat(Q)(a ket(u) + b ket(v)) = a hat(Q) ket(u) + b hat(Q) ket(v)
    $
    其中$a, b in CC$是任意复数。
  - 反线性幺正算符$hat(Q)$满足反线性条件
    $
      hat(Q)(a ket(u) + b ket(v)) = a^* hat(Q) ket(u) + b^* hat(Q) ket(v)
    $
    其中$a, b in CC$是任意复数，$a^*$表示复共轭。

  Dirac符号中，都认为$hat(Q)$是线性幺正算符
  $
    braket(psi_1, hat(A), psi_2) = bra(psi_1) (hat(A) ket(psi_2)) = (bra(psi_1) hat(A)^dagger) ket(psi_2)
  $
  而反线性算符不能用这种方式表示，其向左和向右作用是不一样的。
]
#newpara()

#proof[
  任意选取一组正交归一和完备的基矢$ket(phi_1), ket(phi_2), ket(phi_3), · · ·$，在所说的变换下，得到一组对应的态矢$ket(phi'_1), ket(phi'_2), ket(phi'_3), · · ·$，即$ket(phi_n) ->^cal(Q) ket(phi'_n)$。

  *第一步：*证明$ket(phi'_n)$也是正交归一的。根据$braket(phi_m, phi_n) = delta_(m n)$，以及变换保持内积模方不变得到
  $
    abs(braket(phi'_m, phi'_n))^2 = abs(braket(phi_m, phi_n))^2 = delta_(m n)
  $
  利用$braket(phi'_m, phi'_m) ≥ 0$，即可得到
  $
    braket(phi'_m, phi'_n) = delta_(m n)
  $
  #newpara()

  *第二步：*证明${ket(phi'_n)}$也是完备的。采用反证法，假设存在非零态矢 $ket(w')$ 与所有 $ket(phi'_n)$ 正交，利用 $abs(braket(phi_n, w))^2 = abs(braket(phi'_n, w'))^2$ 得到$braket(phi_n, w) = 0$。由于 $ket(phi_n)$ 是完备的，所以 $ket(w) = 0$，从而 $ket(w') = 0$，  与假设矛盾。

  *第三步：*选择基矢$ket(phi'_n)$的相位。对于任意$n != 1$，定义
  $
    ket(alpha_n) = ket(phi_1) + ket(phi_n)
  $
  变换后的态矢$ket(alpha'_n)$必然可以用基矢$ket(phi'_n)$展开，即
  $
    ket(alpha'_n) = sum_m ket(phi'_m) braket(phi'_m, alpha'_n)
  $
  由于
  $
    abs(braket(phi'_m, alpha'_n))^2 = abs(braket(phi_m, alpha_n))^2
  $
  所以$ket(alpha'_n)$的展开式也只包含$m = 1$和$m = n$两个分量，并且展开系数是模为$1$的相位因子，即
  $
    ket(alpha'_n) = e^(i theta_1) ket(phi'_1) + e^(i theta_n) ket(phi'_n)
  $
  即
  $
    e^(-i theta_1) ket(alpha'_n) = ket(phi'_1) + e^(i (theta_n - theta_1)) ket(phi'_n)
  $
  将相位因子吸收入态矢量，定义$(n!=1)$
  $
    ket(alpha''_n) = e^(-i theta_1) ket(alpha'_n), ket(phi''_n) = e^(i (theta_n - theta_1)) ket(phi'_n)
  $
  有
  $
    ket(alpha''_n) = ket(phi'_1) + ket(phi''_n)
  $
  重新选取基矢${ket(phi'_n)}$的相位，将新的基矢记为$ket(phi''_n)$，其中
  $
    ket(phi''_1) = ket(phi'_1), ket(phi''_n) = e^(i (theta_n - theta_1)) ket(phi'_n) (n != 1)
  $
  请注意*$ket(phi''_1)$的选取是任意的*，这里为了方便把它的下标设定为1(完全可以换一个下标，例如$l$)，作为相位调整的基准。

  *第四步：*确定新旧展开系数之间的关系。设 $ket(u)$ 是系统的任意态矢，$ket(u) = hat(Q) ket(u')$。变换前后的态矢可以展开为
  $
    ket(u) = sum_n c_m ket(phi_n) = c_1 ket(phi_1) + sum_(m != 1) c_m ket(phi_m)\
    ket(u') = sum_n c'_m ket(phi''_m) = c'_1 ket(phi''_1) + sum_(m != 1) c'_m ket(phi''_m)
  $
  研究$c_m$和$c'_m$ 之间的关系。利用任意两个态矢的内积模方不变，得
  $
    abs(braket(phi_m, u))^2 = abs(braket(phi'_m, u'))^2 = abs(braket(phi''_m, u'))^2
  $
  由此得到(对任意$m$)
  $
    abs(c_m)^2 = abs(c'_m)^2
  $
  态矢$ket(u)$的整体相位可以自由选择。设 $c_1 = abs(c_1) e^(- i gamma_1)$，定义新的态矢
  $
    e^(i gamma_1) ket(u) = abs(c_1) ket(phi_1) + sum_(m != 1) c_m e^(i gamma_1) ket(phi_m)
  $
  因此，我们可以适当选择相位因子，*使得$c_1$为实数*。这相当于重新定义$ket(phi_1)$和$c_m (m != 1)$。态矢$ket(u')$的整体相位也可以自由选择。由于$abs(c_1) = abs(c'_1)$，所以$c_1/c'_1$是一个模为 1 的复数。令$c_1/c'_1 = e^(i delta_1)$，可将$ket(u')$的展开式改写为
  $
    e^(i delta_1) ket(u') = c_1 ket(phi''_1) + sum_(m != 1) c'_m e^(i delta_1) ket(phi''_m)
  $
  因此，我们可以适当选择相位因子，*使得$c'_1 = c_1$*。这相当于重新定义$ket(u')$和$c'_m (m != 1)$，实际上也相当于重新定义变换算符。

  接下来，直接计算得到($n != 1$)
  $
    braket(alpha_n, u) = braket(phi_1, u) + braket(phi_n, u) = c_1 + c_n\
    braket(alpha''_n, u') = braket(phi''_1, u') + braket(phi''_n, u') = c'_1 + c'_n
  $
  利用内积模方不变，得到
  $
       & abs(c_1 + c_n)^2 = abs(c'_1 + c'_n)^2 \
    => & abs(c_1)^2 + abs(c_n)^2 + c_1 c_n^* + c_1 c_n = abs(c'_1)^2 + abs(c'_n)^2 + c'_1 c'^*_n + c'^*_1 c'_n \
  $
  由于$abs(c_m)^2 = abs(c'_m)^2$ (对任意 $m$)，所以
  $
    c_1 c_n^* + c_1 c_n = c'_1 c'^*_n + c'^*_1 c'_n
  $
  根据前面的相位因子选择，有
  $
    c_1 = c'_1 = c_1^* = c'^*_1
  $
  所以
  $
    c_n + c_n^* = c'_n + c'^*_n
  $
  设
  $
    c_n = a_n + i b_n, c'_n = a'_n + i b'_n
  $
  其中$a_n, b_n, a'_n, b'_n$均为实数。则上式给出
  $
    a_n = a'_n
  $
  再利用$abs(c_n)^2 = abs(c'_n)^2$，得到
  $
    b_n = plus.minus b'_n
  $
  + 对于$b'_n = b_n$的情形，我们得到
    $
      c'_m = c_m, forall m\
      ket(u') = sum_m c_m ket(phi''_m)
    $
    对应的变换为*线性幺正变换*
  + 对于$b'_n = - b_n$的情形，我们得到
    $
      c'_m = c_m^*, forall m\
      ket(u') = sum_m c^*_m ket(phi''_m)
    $
    对应的变换为*反线性幺正变换*

  接下来要解决的问题是：变换算符$hat(Q)$如何定义？对于一个算符，只要知道了它作用在任意一组基矢上的结果，就给出了这个算符的定义。我们将变换算符$hat(Q)$定义为
  $
    hat(Q) ket(phi_n) = ket(phi''_n)
  $
  同时，我们适当选择相位因子，使得对于任意态矢$ket(u)$，有
  $
    hat(Q) ket(u) = ket(u')
  $
  将$ket(u)$和$ket(u')$的展开式带入上式，对于线性幺正变换，有
  $
    hat(Q) (sum_m c_m ket(phi_m)) = sum_m c_m hat(Q) ket(phi_m) = sum_m c_m ket(phi''_m)
  $
  对于反线性幺正变换，有
  $
    hat(Q) (sum_m c_m ket(phi_m)) = sum_m c^*_m hat(Q) ket(phi_m) = sum_m c^*_m ket(phi''_m)
  $
]

事实上我们还得说明这个变换是线性变换。

*线性幺正变换*： 对应情形(1)，通常简称为幺正变换。

首先我们证明 $hat(Q)$ 是线性变换，即：对于任意态矢$ket(psi_1)$、$ket(psi_2)$和任意复数$c_1$、$c_2$有
$
  hat(Q) (c_1 ket(psi_1) + c_2 ket(psi_2)) = c_1 hat(Q) ket(psi_1) + c_2 hat(Q) ket(psi_2)
$
#proof[
  将$ket(psi_1)$和$ket(psi_2)$展开为
  $
    ket(psi_1) = sum_m x_m ket(phi_m), ket(psi_2) = sum_m y_m ket(phi_m)
  $
  带入计算，得到
  $
      & hat(Q) (c_1 ket(psi_1) + c_2 ket(psi_2)) = hat(Q) sum_m (c_1 x_m + c_2 y_m) ket(phi_m) \
    = & sum_m (c_1 x_m + c_2 y_m) hat(Q) ket(phi_m) = sum_m (c_1 x_m + c_2 y_m) ket(phi''_m) \
    = & c_1 sum_m x_m ket(phi''_m) + c_2 sum_m y_m ket(phi''_m) = c_1 ket(psi'_1) + c_2 ket(psi'_2) \
    = & c_1 hat(Q) ket(psi_1) + c_2 hat(Q) ket(psi_2)
  $
]
接下来我们证明：*幺正变换保证任意两个态矢之间的内积不变*(不仅仅是模方)，
$
  braket(psi_1, psi_2) = braket(psi'_1, psi'_2)
$
证明是直截了当的，利用展开式直接计算即可。
#proof[
  计算过程如下
  $
    ket(psi_1) = sum_m x_m ket(phi_m), ket(psi_2) = sum_m y_m ket(phi_m)\
    ket(psi'_1) = sum_m x_m ket(phi''_m), ket(psi'_2) = sum_m y_m ket(phi''_m)\
    => braket(psi_1, psi_2) = sum_m x^*_m y_m = braket(psi'_1, psi'_2)
  $
]
最后我们证明：$hat(Q)$ 是幺正算符，即 $Q^(-1) = Q^dagger$。对于任意的线性算符$hat(L)$，我们有
$
  bra(psi_1) (hat(L) ket(psi_2)) = (bra(psi_1) hat(L)) ket(psi_2) = braket(psi_1, hat(L), psi_2)
$
即：对于线性算符，不用区分向前还是向后作用。为了与后面的反线性算符对比，我们简要说明这一性质。利用展开式以及$hat(L)$是线性算符以及内积的线性性质，得到
$
  bra(psi_1) (hat(L) ket(psi_2)) & = (sum_m x^*_m bra(phi_m)) (hat(L) sum_n y_n ket(phi_n)) \
                                 & = (sum_m x^*_m bra(phi_m)) (sum_n y_n hat(L) ket(phi_n)) \
                                 & = sum_(m, n) x^*_m y_n bra(phi_m) (hat(L) ket(phi_n)) \
  (bra(psi_1) hat(L)) ket(psi_2) & = sum_(m, n) x^*_m y_n (bra(phi_m) hat(L)) ket(phi_n) \
$
两者相等实际上是基于Dirac的*结合律约定*：对于线性算符，先向前和先向后作用得到相同的结果
$
  bra(phi_m) (hat(L) ket(phi_n)) = (bra(phi_m) hat(L)) ket(phi_n)
$
#proof[
  由于$hat(Q)$是线性算符，其Hermite共轭算符$hat(Q)^dagger$亦为线性算符，所以
  $
    braket(psi'_1, psi'_2) &= (bra(psi_1) hat(Q)^dagger) (hat(Q) ket(psi_2)) = bra(psi_1) (hat(Q)^dagger hat(Q) ket(psi_2))\
    &= braket(psi_1, hat(Q)^dagger hat(Q), psi_2)\
    &= braket(psi_1, psi_2)
  $
  由于$ket(psi_1)$和$ket(psi_2)$是任意态矢，我们得到
  $
    hat(Q)^dagger hat(Q) = 1
  $
  由于$hat(Q)$变换是可逆的(定理的条件设定)，所以$hat(Q)^(-1) = hat(Q)^dagger$，即 $hat(Q)$是幺正算符(同时也得到 $hat(Q) hat(Q)^dagger = 1$)。
]

*反线性幺正变换*：对应情形(2)，通常简称为反幺正变换。

首先证明$hat(Q)$是反线性变换，即：对于任意态矢$ket(psi_1)$、$ket(psi_2)$和任意复数$c_1$、$c_2$，有
$
  hat(Q) (c_1 ket(psi_1) + c_2 ket(psi_2)) = c_1^* hat(Q) ket(psi_1) + c_2^* hat(Q) ket(psi_2)
$
#proof[
  将$ket(psi_1)$和$ket(psi_2)$展开为
  $
    ket(psi_1) = sum_m x_m ket(phi_m), ket(psi_2) = sum_m y_m ket(phi_m)
  $
  带入计算，得到
  $
      & hat(Q) (c_1 ket(psi_1) + c_2 ket(psi_2)) = hat(Q) sum_m (c_1 x_m + c_2 y_m) ket(phi_m) \
    = & sum_m (c_1 x_m + c_2 y_m)^* hat(Q) ket(phi_m) = sum_m (c_1^* x^*_m + c_2^* y^*_m) ket(phi''_m) \
    = & c_1^* sum_m x^*_m ket(phi''_m) + c_2^* sum_m y^*_m ket(phi''_m) = c_1^* ket(psi'_1) + c_2^* ket(psi'_2) \
    = & c_1^* hat(Q) ket(psi_1) + c_2^* hat(Q) ket(psi_2)
  $
]
接下来我们证明：*反幺正变换将任意两个态矢之间的内积变换为其复共轭*
$
  braket(psi_1, psi_2) = braket(psi'_1, psi'_2)^*
$
证明是直截了当的，利用展开式直接计算即可。
#proof[
  计算过程如下
  $
    ket(psi_1) = sum_m x_m ket(phi_m), ket(psi_2) = sum_m y_m ket(phi_m)\
    ket(psi'_1) = sum_m x^*_m ket(phi''_m), ket(psi'_2) = sum_m y^*_m ket(phi''_m)\
    => braket(psi_1, psi_2) = sum_m x^*_m y_m = braket(psi'_1, psi'_2)^*
  $
]
最后我们证明：$hat(Q)$是幺正算符，即 $hat(Q)^(-1) = hat(Q)^dagger$。与线性算符不同，对于任意的反线性算符$hat(A)$，我们约定
$
  bra(psi_1) (hat(A) ket(psi_2)) = ((bra(psi_1) hat(A)) ket(psi_2))^*
$
这意味着，*对于反线性算符，Dirac结合律约定通常是不成立的*，需要特别说明反线性算符是向前作用还是向后作用。利用展开式，以及$hat(L)$是反线性算符以及内积的线性性质，得到
$
    bra(psi_1) (hat(A) ket(psi_2)) & = (sum_m x^*_m bra(phi_m)) (hat(A) sum_n y_n ket(phi_n)) \
                                   & = (sum_m x^*_m bra(phi_m)) (sum_n y^*_n hat(A) ket(phi_n)) \
                                   & = sum_(m, n) x^*_m y^*_n bra(phi_m) (hat(A) ket(phi_n)) \
  ((bra(psi_1) hat(A)) ket(psi_2)) & = sum_(m, n) x_m y_n ((bra(phi_m) hat(A)) ket(phi_n)) \
$
约定$bra(phi_m)(hat(A) ket(phi_n)) = ((bra(phi_m) hat(A)) ket(phi_n))^*$，我们就自洽地得到上面两个式子相等。

#proof[
  由于$hat(Q)$是反线性算符，其Hermite共轭算符$hat(Q)^dagger$亦为反线性算符，但是$hat(Q)^dagger hat(Q)$则为线性算符。所以
  $
    braket(psi'_1, psi'_2) & = ((bra(psi_1) hat(Q)^dagger) (hat(Q) ket(psi_2))) = (bra(psi_1) (hat(Q)^dagger hat(Q) ket(psi_2)))^*\
    & = braket(psi_1, hat(Q)^dagger hat(Q), psi_2)^*\
    & = braket(psi_1, psi_2)^*
  $
  由于$ket(psi_1)$和$ket(psi_2)$是任意态矢，我们得到
  $
    hat(Q)^dagger hat(Q) = 1
  $
  由于$hat(Q)$变换是可逆的(定理的条件设定)，所以$hat(Q)^(-1) = hat(Q)^dagger$，即 $hat(Q)$是幺正算符(同时也得到 $hat(Q) hat(Q)^dagger = 1$)。
]

接下来要解决的问题是：*在对称性变换下，任意力学量$hat(Ω)$如何变换？*根据物理观测结果不变的要求，力学量在变换前后本征值和取值几率都不变，所以其平均值必然是不变的，即
$
  braket(psi, hat(Omega), psi) = braket(psi', hat(Omega)', psi')
$
根据态矢的变换关系得到
$
  (bra(psi) hat(Q)^dagger) hat(Omega)' (hat(Q) ket(psi)) = braket(psi, hat(Omega), psi)
$
+ $hat(Q)$是线性算符。此时可将上式写为
  $
    bra(psi) (hat(Q)^dagger hat(Omega)' hat(Q)) ket(psi) = braket(psi, hat(Omega), psi)
  $
  由于态矢$ket(psi)$是任意的，所以
  $
    hat(Q)^dagger hat(Omega)' hat(Q) = hat(Omega) => hat(Omega)' = hat(Q) hat(Omega) hat(Q)^dagger
  $
+ $hat(Q)$是反线性算符。此时上式写为
  $
    (bra(psi) hat(Q)^dagger) hat(Omega)' (hat(Q) ket(psi)) = (bra(u) (hat(Q)^dagger hat(Omega)' hat(Q) ket(psi)))^*
  $
  由于$hat(Q)^dagger hat(Omega)' hat(Q)$是线性算符且是Hermite算符，所以得到
  $
    (bra(u) hat(Q)^dagger) hat(Omega)' (hat(Q) ket(u)) = braket(u, (hat(Q)^dagger hat(Omega)' hat(Q)), u)^* = braket(u, (hat(Q)^dagger hat(Omega)' hat(Q)), u)
  $
  最后一个等号成立是由于前后的态矢都是$ket(u)$，即Hermite算符的平均值必为实数。由于态矢$ket(u)$是任意的，所以仍然得到
  $
    hat(Q)^dagger hat(Omega)' hat(Q) = hat(Omega) => hat(Omega)' = hat(Q) hat(Omega) hat(Q)^dagger
  $

#newpara()

*总结：量子系统的对称性变换包含了对态矢和力学量的变换。*

#proposition(subname: [量子系统的对称性变换])[
  对态矢$ket(psi)$的变换为
  $
    ket(psi) ->^(cal(Q)) ket(psi') = hat(Q) ket(psi)
  $
  对力学量$hat(Omega)$的变换为
  $
    hat(Omega) ->^(cal(Q)) hat(Omega)' = hat(Q) hat(Omega) hat(Q)^dagger = hat(Q) hat(Omega) hat(Q)^(-1)
  $
  有
  $
    hat(Omega) ket(psi) ->^(cal(Q)) hat(Omega)' ket(psi') = hat(Q) hat(Omega) hat(Q)^dagger (hat(Q) ket(psi)) = hat(Q) (hat(Omega) ket(psi))
  $
]

#newpara()
注意：*目前仅从物理观测结果不变出发，得到了量子系统可能的对称性变换的性质。这并不意味着量子系统一定具有该变换下的对称性。*

例如，任意系统都可以进行空间平移变换，这并不表明该系统具有空间平移不变性。

== 对称性群

现在来解决第二个大问题：*什么叫做量子系统具有某种对称性？*或者具体一点，如何定义量子系统在某种对称性变换下的不变性。

=== 量子系统的不变性

量子系统的运动方程是Schrödinger方程。量子系统具有某种对称性，可以定义为*Schrödinger方程在此对称性变换下是不变的*，即变换前后的态矢遵循同样的演化方程。变换前的态矢$ket(psi)$遵守Schrödinger方程
$
  i hbar pdv(, t) ket(psi(t)) = hat(H) ket(psi(t))
$
考虑对称性变换$hat(Q)$，在最一般的情形下，它可能是含时的。很容易推导出变换后的态矢$ket(psi') = hat(Q) ket(psi)$满足的运动方程为
$
  i hbar pdv(, t) ket(psi'(t)) & = i hbar pdv(, t) (hat(Q) ket(psi(t))) \
                               & = i hbar (pdv(hat(Q), t) ket(psi(t)) + hat(Q) pdv(, t) ket(psi(t))) \
                               & = i hbar pdv(hat(Q), t) ket(psi(t)) + hat(Q) hat(H) ket(psi(t)) \
                               & = (hat(Q) hat(H) hat(Q)^dagger + i hbar pdv(hat(Q), t) hat(Q)^dagger) ket(psi'(t))
$
因此，如果
$
  hat(H)' = hat(Q) hat(H) hat(Q)^dagger + i hbar pdv(hat(Q), t) hat(Q)^dagger = hat(H)\
  <=> i hbar pdv(hat(Q), t) = [hat(H), hat(Q)]
$
则变换后的态$ket(psi') = hat(Q) ket(psi)$也满足同样的Schrödinger方程
$
  i hbar pdv(, t) ket(psi'(t)) = hat(H) ket(psi'(t))
$
此时我们称系统具有$hat(Q)$变换下的不变性。可见，*是否具有这种对称性取决于系统的Hamilton量*。

在绝大部分情况下，我们所研究的对称性变换是不含时的，此时退化为
$
  [hat(H), hat(Q)] = 0
$
即：系统的Hamilton量和$hat(Q)$对易，则系统具有该变换下的*不变性*。

=== 对称性变换群

在数学上，某一类对称性变换的集合${hat(Q)_i}$可以构成一个*变换群*，群元的乘法很自然地定义为算符的乘法。设$ket(j)$为任意态矢，则
$
  (hat(Q)_i hat(Q)_j) ket(j) ≡ hat(Q)_i (hat(Q)_j ket(j))
$
根据群论，这类变换需要满足如下条件：
- *封闭性*：对于任意两个变换$hat(Q)_i$和$hat(Q)_j$，它们的乘积$hat(Q)_i hat(Q)_j$仍然是该类变换中的一个变换。
- *结合律*：对于任意三个变换$hat(Q)_i$、$hat(Q)_j$和$hat(Q)_k$，有
  $
    (hat(Q)_i hat(Q)_j) hat(Q)_k = hat(Q)_i (hat(Q)_j hat(Q)_k)
  $
- *单位元*：存在单位变换$hat(I)$，使得对于任意变换$hat(Q)_i$，有
  $
    hat(I) hat(Q)_i = hat(Q)_i hat(I) = hat(Q)_i
  $
  这个单位元就是单位算符 (恒等变换)。
- *逆元*：对于任意变换$hat(Q)_i$，存在唯一逆变换$hat(Q)_i^(-1)$，使得
  $
    hat(Q)_i hat(Q)_i^(-1) = hat(Q)_i^(-1) hat(Q)_i = hat(I)
  $


*连续对称性变换*：连续对称性变换$hat(U)$总可以从恒等变换出发，连续作用无穷小变换得到，所以必然是*幺正变换*。考虑无穷小变换
$
  hat(U) = 1 - (i epsilon)/hbar hat(G) + O(epsilon^2)
$
有限变换可以通过无穷小变换连续作用得到$(alpha in RR)$
$
  hat(U)(alpha) = lim_(N -> oo) [1 - (i alpha)/(N hbar) hat(G)]^N = e^(- i/hbar alpha hat(G))
$
根据幺正性要求$hat(U)^dagger hat(U) = 1$得到
$
  hat(G)^dagger = hat(G)
$
所以$hat(G)$是*Hermite算符*，对应系统的某个力学量。我们将$hat(G)$称为对称性变换$hat(U)$的*生成元*。

可以验证，当$alpha$取遍所有实数时，所有的$hat(U)(alpha)$变换构成的集合${hat(U){(alpha)}}$构成一个群 (连续群)。群元之间的乘法是可交换的，因此是一个Abel群。如果$hat(U)$是不含时的，系统具有$hat(U)$变换下的不变性意味着
$
  [hat(H), hat(U)] = 0 <=> hat(U) hat(H) hat(U)^dagger = hat(H) <=> [hat(H), hat(G)] = 0
$
这意味着生成元$hat(G)$是守恒量。

在$hat(U)$变换下，任意力学量$hat(Omega)$变为$hat(Omega)' = hat(U) hat(Omega) hat(U)^dagger$，利用无穷小变换或者Baker-Campbell-Hausdorff公式，我们得到
$
  hat(Omega)' - hat(Omega) = i/ hbar alpha [hat(G), hat(Omega)] + O(alpha^2)
$
如果对基本力学量的变换做出要求，则可以得到基本对易关系(见空间平移变换)。

上述讨论的是单参数连续对称性变换，可以很自然地推广到多个参数的情形，此时有限变换可以用一组参数${alpha_a in RR}$一组生成元${hat(G)_a}$表达为(默认对$a = 1, 2,..., n$求和，即Einstein求和约定)
$
  hat(U)(alpha_a) = e^(- i/hbar alpha_a hat(G)_a)
$
无穷小变换为
$
  hat(U)(epsilon_a) = 1 - (i epsilon_a)/hbar hat(G)_a + O(epsilon_a^2)
$
根据幺正性要求得到
$
  hat(U)^dagger hat(U) = 1 + 1/hbar epsilon_a (hat(G)_a^dagger - hat(G)_a) + O(epsilon_a^2) <=> hat(G)_a^dagger = hat(G)_a
$
所以每个生成元$hat(G)_a$都是Hermite算符
$
  hat(G)_a^dagger = hat(G)_a
$
当$alpha_a$取遍所有实数时，所有的$hat(U)(alpha_a)$变换构成的集合${hat(U)(alpha_a)}$也可能构成一个群。如果生成元$hat(G)_a$之间是两两对易的，则为$n$个单参数群的直积，群元之间的乘法是可交换的。

如果生成元$hat(G)_a$之间并不两两对易，则群元之间的乘法一般不可交换，即非Abel群。

要构成非Abel群，生成元$hat(G)_a$之间的对易关系要满足一定的条件。考虑任意群元$hat(U)(alpha_a)$和$hat(U)(beta_a)$，其一阶无穷小形式为
$
  hat(U)(alpha_a) = 1 - (i alpha_a)/hbar hat(G)_a + O(alpha_a^2)\
  hat(U)(beta_a) = 1 - (i beta_a)/hbar hat(G)_a + O(beta_a^2)
$
很显然，精确到一阶无穷小，$hat(U)(alpha_a)$和$hat(U)(beta_a)$是可交换的。因此，要研究非Abel效应，必须考虑二阶无穷小项。将群元$hat(U)(alpha_a)$展开到二阶无穷小
$
  hat(U)(alpha_a) = 1 - (i alpha_a)/hbar hat(G)_a - 1/(2 hbar^2) alpha_a alpha_b hat(G)_a hat(G)_b + ...\
  hat(U)^(-1)(alpha_a) = 1 + (i alpha_a)/hbar hat(G)_a - 1/(2 hbar^2) alpha_a alpha_b hat(G)_a hat(G)_b + ...
$
则可以计算 (到二阶项)
$
  hat(U)^(-1)(alpha_a) hat(U)^(-1)(beta_a) hat(U)(alpha_a) hat(U)(beta_a) & = 1 - 1/(hbar^2) alpha_a beta_b [hat(G)_a, hat(G)_b] + ...
$
这是四个群元相乘，仍然是某个群元，所以(注意要对$c$求和)
$
  [hat(G)_a, hat(G)_b] = i hbar c_(a b c) hat(G)_c
$
即生成元之间的对易关系(*代数*)必须是封闭的。

#example(subname: [空间平移变换])[
  先考虑最简单的一维情形，沿$x$方向的平移。此时，生成元$hat(G)$即为$hat(p)_x$。沿$x$方向平移$a$，对应的变换算符为
  $
    hat(U)(a) = e^(- i/hbar a hat(p)_x)
  $
  在此变换下，坐标算符$hat(x)$变换为
  $
    hat(x)' = hat(U)^dagger(a) hat(x) hat(U)(a) = e^(-i/hbar a hat(p)_x) hat(x) e^(i/hbar a hat(p)_x) = hat(x) - a
  $
  计算中利用了Baker-Campbell-Hausdorff公式和基本对易关系$[hat(x), hat(p)_x] = i hbar$。

  反过来，若要求$hat(U)(a)$为平移算符，使得$hat(x)' = hat(x) - a$，则得到基本对易关系$[hat(x), hat(p)_x] = i hbar$。

  三维空间的变换则为多参数连续变换，对应的变换算符为
  $
    hat(U)(vb(a)) = e^(- i/hbar sum_(i = x, y, z) a_i p_i) = e^(- i/hbar vb(a) dot vb(hat(p)) )
  $
  由于$[hat(p)_i, hat(p)_j] = 0$，这实际上是三个方向的平移群的直积。可以计算出 $hat(U)(vb(a))$ 对坐标表象基矢的作用为
  $
    hat(U)(vb(a)) ket(vb(x)) & = integral dd(vb(p)) hat(U)(vb(a)) ket(vb(p)) braket(vb(p), vb(x)) \
                             & = integral dd(vb(p)) e^(- i/hbar vb(a) dot vb(p)) ket(vb(p)) e^(- i/hbar vb(p) dot vb(x)) \
                             & = ket(vb(x) + vb(a))
  $
  系统状态的变换为$ket(psi) -> hat(U)(vb(a)) ket(psi)$，所以波函数的变换为
  $
    psi(vb(x)) = braket(vb(x), psi) -> psi'(vb(x)) = braket(vb(x), hat(U)(vb(a)), psi) = braket(vb(x) + vb(a), psi) = psi(vb(x) + vb(a))
  $
]

#example(subname: [空间转动变换])[
  空间转动变换的生成元为系统的角动量。对于绕$z$轴转动角度$phi$的操作，转动变换算符为
  $
    hat(U)_3 (phi) = e^(- i/hbar phi hat(J)_3)
  $
  对于一般的转动变换，变换算符可写为
  $
    hat(U)(phi_1, phi_2, phi_3) = e^(- i/hbar sum_(i = 1)^3 phi_i hat(J)_i)
  $
  这些变换构成$"SO"(3)$群，生成元的对易关系即为群的代数
  $
    [hat(J)_i, hat(J)_j] = i hbar epsilon_(i j k) hat(J)_k
  $
]

#example(subname: [时间平移变换])[
  根据前面的讨论，我们可能会认为时间平移变换的生成元是系统的Hamilton量，时间平移$tau$对应的变换算符为
  $
    hat(U)(tau) = e^(i/hbar tau hat(H))\
    hat(U)(tau) ket(psi(t)) = ket(psi(t - tau))
  $
  取$tau = t$，则得到时间演化$ket(psi(t)) = e^(-i/hbar hat(H) t) ket(psi(0))$。

  这说明上面写出的算符只有在体系真实的时间演化中才成立，并不是一般意义上的时间平移。实际上，时间平移算符可以定义为
  $
    hat(U)(tau) = e^(- tau pdv(, t))
  $
  显然，$hat(H) = i hbar pdv(, t)$只有对真实的时间演化才成立。
]

#definition(subname: [量子系统的对称性群])[
  我们将“量子系统具有某种变换下的不变性” 定义为“Schrödinger方程在变换下形式不变”，由此我们推导出变换$hat(Q)$满足方程
  $
    i hbar pdv(hat(Q), t) = [hat(H), hat(Q)] <=> hat(Q) hat(H) hat(Q)^dagger= hat(H) - i hbar pdv(hat(Q), t) hat(Q)^dagger
  $
  对于给定Hamilton量$hat(H)$的量子系统，所有保持系统的运动方程 (即方程)形式不变的对称性变换$hat(Q)$构成一个群，称为系统的*一般对称性群*。
]

#proof[
  + 封闭性：设$hat(Q)_1$和$hat(Q)_2$均满足方程，则
    $
      [hat(H), hat(Q)_1 hat(Q)_2] & = hat(Q)_1 [hat(H), hat(Q)_2] + [hat(H), hat(Q)_1] hat(Q)_2 \
                                  & = i hbar hat(Q)_1 pdv(hat(Q)_2, t) + i hbar pdv(hat(Q)_1, t) hat(Q)_2 \
                                  & = i hbar pdv(hat(Q)_1 hat(Q)_2, t)
    $
    所以$hat(Q)_1 hat(Q)_2$也满足方程。
  + 结合律：算符的乘法满足结合律。
  + 单位元：单位算符$hat(I)$满足方程，因为$[hat(H), hat(I)] = 0$且$pdv(hat(I), t) = 0$。
  + 逆元：设$hat(Q)$满足方程，则
    $
      0 = [hat(Q) hat(Q)^(-1), hat(H)] & = hat(Q) [hat(Q)^(-1), hat(H)] + [hat(Q), hat(H)] hat(Q)^(-1) \
                                       & = hat(Q) [hat(Q)^(-1), hat(H)] - i hbar pdv(hat(Q), t) hat(Q)^(-1) \
                                       & = hat(Q) [hat(Q)^(-1), hat(H)] + hat(Q) i hbar pdv(hat(Q)^(-1), t)
                                         => [hat(Q)^(-1), hat(H)] = i hbar pdv(hat(Q)^(-1), t)
    $
    所以$hat(Q)^(-1)$也满足方程。
]

含时的一般对称性群只对Hamilton量含时的系统有意义。通常我们只关心*Hamilton量不含时*的量子系统，此时只需要考虑不含时的对称性变换$hat(Q)$。系统在变换下具有不变性意味着
$
  [hat(H), hat(Q)] = 0 <=> hat(Q) hat(H) hat(Q)^dagger = hat(H)
$
满足上述条件的所有*不含时*对称性变换$hat(Q)$构成系统的对称性群。这个结论显然是含时情况的特例，无需再次证明。

对于Hamilton量不含时的系统，求解$hat(H)$的本征方程 (定态方程) 得到系统的能量本征值 (能级) 和正交归一能量本征态：
$
  hat(H) ket(phi_(n i)) = E_n ket(phi_(n i))\
  braket(phi_(n i), phi_(m j)) = delta_(n m) delta_(i j)
$
设能量本征值$E_n$的简并度为$f_n$，即$i = 1, 2, ... , f_n$。如果Hamilton量在某种变换下保持不变，那么群表示论就派上用场了。

== 对称性和简并

假设系统具有$hat(Q)$变换的不变性，那么$hat(H) hat(Q) = hat(Q) hat(H)$，所以
$
  hat(H) hat(Q) ket(phi_(n i)) = hat(Q) hat(H) ket(phi_(n i)) = E_n (hat(Q) ket(phi_(n i)))
$
因此，$hat(Q) ket(psi_(n i))$亦是本征值为$E_n$的能量本征态。

如果 $E_n$ 没有简并，那么必然有 $hat(Q) ket(phi_(n i)) = e^(i theta_n) ket(phi_(n i))$。不失一般性，可取 $theta_n = 0$。

如果$E_n$有简并，则$hat(Q) ket(phi_(n i))$ 必然可以写为简并子空间中的 $f_n$ 个基矢的线性叠加，即
$
  hat(Q) ket(phi_(n i)) = sum_(j = 1)^(f_n) M_(j i)^((n))(hat(Q)) ket(phi_(n j))
$
对应的数学语言是：$f_n$维的简并子空间是$hat(Q)$变换的*不变子空间*，构成系统对称性群的一个*表示空间*，这些$f_n times f_n$的矩阵${M^((n))}$就是群的一个$f_n$维矩阵*表示*。对于无简并的能级，表示即为 1。

将展开系数形成的矩阵 ${M^((n))}$ 称为对称性群的表示，其物理含义其实是非常形象的。考虑任意两个变换$hat(Q)^1$和$hat(Q)^2$，则可以计算 (为方便，略去上下标 $n$)
$
  hat(Q)^1 hat(Q)^2 ket(phi_i) & = hat(Q)^1 (sum_j M_(j i)(hat(Q)^2) ket(phi_j)) \
                               & = sum_j M_(j i)(hat(Q)^2) hat(Q)^1 ket(phi_j) \
                               & = sum_(j, k) M_(j i)(hat(Q)^2) M_(k j)(hat(Q)^1) ket(phi_k) \
                               & = sum_k (sum_j M_(k j)(hat(Q)^1) M_(j i)(hat(Q)^2)) ket(phi_k)
$
另一方面，$hat(Q)^1 hat(Q)^2$必然也是对称性群中的元素，所以也可以这么计算
$
  hat(Q)^1 hat(Q)^2 ket(phi_i) = sum_k M_(k i)(hat(Q)^1 hat(Q)^2) ket(phi_k)
$
利用矩阵乘法的定义，我们得到
$
  M(hat(Q)^1 hat(Q)^2) = M(hat(Q)^1) M(hat(Q)^2)
$
进一步可以验证，矩阵集合${M(hat(Q))}$按照矩阵乘法也构成一个群(*矩阵群*)，与系统的对称性群*同态*。因此，将矩阵集合${M(hat(Q))}$称为对称性群的一个*矩阵表示*，*表示空间即为简并能量本征态构成的简并子空间，表示的维数即为简并子空间的维数*。

#definition(subname: [群的表示])[
  群的表示在数学上是一个抽象概念，这里采用一种形象的表述。假设群$cal(G)$是由某些*变换操作*${hat(R)}$构成的，这些变换操作*作用的对象是一个线性空间*$V$。假设$V$的某个$m$维子空间$V_s$对于操作${hat(R)}$是封闭的(即不变子空间)，在这个$m$维子空间中取一组基矢${w_1k, w_2, ... , w_m}$，则对于群中的任意变换操作$hat(R)$必然有
  $
    hat(R) w_i = sum_(j = 1)^m M_(j i)(hat(R)) w_j
  $
]
我们称$V_s$是群的一个$m$维*表示空间*，这些$m times m$的矩阵${M(hat(R))}$就是群的一个$m$维*矩阵表示*。由于$V_s$中基矢的取法并不唯一，另一种不同的取法将给出不同的表示矩阵$M'$，两者之间必然可以通过相似变换联系起来，$M' = S^(-1) M S$，称两个表示是*等价*的。

通俗地说：不变子空间$V_s$*荷载*了群$G$的一个$m$维表示。

要将群表示与量子系统的对称性联系起来，则需要引入不可约表示的概念。
#definition(subname: [可约表示])[
  假设$V_a$和$V_b$是$V$的两个不相交的不变子空间，维数分别为$m_a$和$m_b$，它们都可以荷载群的表示，即$m_a$维表示${M_a}$和$m_b$维表示${M_b}$。取两个子空间的*直和*$V_s = V_a plus.o V_b$，则它也是群$cal(G)$的$m_a + m_b$维表示空间，表示矩阵亦为矩阵的直和 (*直和表示*)
  $
    M_s = M_a plus.o M_b = mat(M_a, 0; 0, M_b)
  $
  如果群的表示${M}$存在与之等价的直和表示 $M tilde.eq M_a plus.o M_b plus.o dots$，
  则称${M}$是(完全)*可约表示*，表示空间也可以分解为若干个表示空间的直和$V_a plus.o_b plus.o dots$。
]
反之，如果群的表示${M}$不存在与之等价的直和表示，则称为*不可约表示*，其表示空间是*不可约表示空间*，不能进一步分解为若干个表示空间的直和。

*所以，研究群的表示，只需要研究它的不等价的不可约表示。*

接下来将探讨第三个大问题：*对称性群带来怎样的物理结果？*

首先，虽然群的表示是数学上的抽象概念，但是我们在量子系统中已经找到了*荷载它的表示空间*。根据
$
  hat(Q) ket(phi_(n i)) = sum_(j = 1)^(f_n) M_(j i)^((n))(hat(Q)) ket(phi_(n j))
$
简并的能量本征态构成的简并子空间是系统对称性群的一个*表示空间*，矩阵集合${M^{(n)}}$就是群的一个$f_n$维*表示*。

然后，利用群表示论可以证明几个定理：

#theorem[
  设$cal(G)$是Hamilton量$hat(H)$的对称性群，构成$cal(G)$的不可约表示空间的能量本征态必然属于同一个能级(拥有相同的能量本征值)。
]

#proof[
  用反证法。

  设$hat(H)$的$l$个本征态${ket(phi_(i))}$构成群$cal(G)$的第$n$个不可约表示的表示空间，而它们的能量本征值并不相同。取最简单的例子，前$l - 1$个属于能级$E_1$，最后一个属于能级$E_2$，即
  $
    E_i = & E_1, i = 1, 2, ... , l - 1 \
    E_i = & E_2, i = l
  $
  将群$cal(G)$中的任意变换算符$hat(Q)$作用到本征方程$hat(H) ket(phi_(l)) = E_2 ket(phi_(l))$的两边，左边和右边可以分别计算为
  $
    hat(Q) hat(H) ket(phi_(l)) & = hat(H) (hat(Q) ket(phi_(l))) = hat(H) (sum_(i = 1)^l M_(i l) ket(phi_(i))) = sum_(i = 1)^(l) M_(i l) E_i ket(phi_(i)) \
  $
  $
    E_2 (hat(Q) ket(phi_(l))) & = sum_(i = 1)^(l) M_(i l) E_2 ket(phi_(i))
  $
  两边作用的结果相等，得到
  $
    sum_(i = 1)^(l) M_(i l) (E_i - E_2) ket(phi_(i)) = 0
  $
  由于$l$个本征态${ket(phi_(i))}$是线性无关的，所以对任意$i = 1, 2, ... , l$都
  $
    M_(i l) (E_i - E_2) = 0
  $
  但是，当$j = 1, 2, ... , l - 1$时，$E_j = E_1 != E_2$，所以 $M_(j l) = 0$。同理，对本征方程$H ket(phi_(j)) = E_1 ket(phi_(j)) (j = 1, 2, ... , l - 1)$ 两边作用 $hat(Q)$，可以证明当 $j = 1, 2, ... , l - 1$ 时， $M_(l j) = 0$。因此 $M$ 是分块对角的形式，是可约表示。这与假设 ${ket(phi_(i))}$ 构成群的不可约表示空间矛盾。所以只能有 $E_1 = E_2$。
]

#theorem[
  设Hamilton量$hat(H)$的对称性群为$cal(G)$，则$hat(H)$的正交归一本征态荷载的表示必为*幺正表示*，即表示矩阵是幺正矩阵。
]

#proof[
  设$m$个正交归一的能量本征态${ket(phi_i)}$荷载了$G$的一个$m$维表示，即 $braket(psi_i, psi_j)=delta_(i j), i, j = 1, 2, dots, m$。由于群变换$hat(Q)$ 保持内积的模方不变，所以$braket(psi'_i, psi'_j)=delta_(i j)$。利用
  $
    ket(psi'_i) = hat(Q) ket(psi_i) = sum_(k = 1)^m M_(k i)(hat(Q)) ket(psi_k)
  $
  即可证明$M(hat(Q))$是幺正矩阵。
]
因此对于对称性群$cal(G)$，只需要研究它的*不等价的不可约幺正表示*就可以了。事实上，物理学中遇到的绝大部分群，其任何一个表示都存在与之等价的幺正表示。

#theorem[
  设Hamilton量$hat(H)$的对称性群为$cal(G)$，正交归一能量本征态${ket(psi_(m mu))}$荷载$cal(G)$的不可约表示$M^((m))(hat(Q))$，正交归一能量本征态${ket(psi_(n nu))}$荷载$cal(G)$的不可约表示$M^((n))(hat(Q))$，即
  $
    hat(Q) ket(psi_(m mu)) = sum_(mu' = 1)^(f_m) M_(mu' mu)^((m))(hat(Q)) ket(psi_(m mu'))\
    hat(Q) ket(psi_(n nu)) = sum_(nu' = 1)^(f_n) M_(nu' nu)^((n))(hat(Q)) ket(psi_(n nu'))
  $
  若两个不可约 (幺正) 表示不等价，则必然有
  $
    braket(psi_(m mu), psi_(n nu)) = delta_(m n) delta_(mu nu) C
  $
]
证明从略。这个定理说明，*两个不等价的不可约幺正表示必然分别由两个正交的子空间所荷载*。

上述定理为我们勾勒出一幅大致的*物理图像：量子系统的能量本征态按照对称性群的不可约表示进行分类*，即：
```
对称性群 -> 不可约表示 -> 能级结构
```
不可约表示的维数是群的固有数学性质。

因此我们可以先不用暴力求解Hamilton量$hat(H)$的本征值和本征态，利用其对称性群的不可约表示就可以对系统的能级简并度做出理论预测，从而与实验结果或者暴力数值结果进行对比。

不过，在粒子物理等领域的研究中，我们更感兴趣的是，能不能根据实验上对能级结构的观测结果反推Hamilton量的对称性群？例如，在粒子物理发展的早期，Gell-Mann就根据实验中观测到的强子质量谱反推出强相互作用具有$"SU"(3)$味对称性。也就是说，我们想把逻辑反过来：
```
能级结构 -> 不可约表示 -> 对称性群
```
想要这个逻辑成立，我们需要证明如下的结论：
#definition[
  系统Hamilton量属于任一本征值的简并子空间，都荷载着其对称性群的一个不可约表示。
]
可惜的是，这个结论无法给出一般性的证明，因此不能成为一个定理。不过，经过对大量物理系统的研究，人们相信它是对的。你可以把这个结论当成一种*信仰*。

简并子空间必然是对称性群的表示空间，但是对某些系统，我们发现表示是可约的。一种解释是*偶然简并*，即例外情况。不过，如果我们坚持信仰，就有了另一种解释：对称性群没有找全，只找到了它的一个*子群*。如果采用完整的对称性群，简并子空间荷载的表示就成为不可约表示了。后面将以氢原子为例进行说明。

#example(subname: [氢原子的动力学对称性])[
  对于在一般中心势场$V(r)$中运动的粒子，Hamilton量的对称性群是三维正当转动群SO(3)。Hamilton量的本征方程的解为
  $
    & hat(H) ket(n l m) = E_n ket(n l m) \
    & psi_(n l m)(vb(r)) = braket(vb(r), n l m) = R_(n l)(r) Y_(l m)(theta, phi) \
    & m = -l, -l + 1, ... , l - 1, l
  $
  能量本征值$E_(n l)$依赖量子数$n$和$l$，*能级是$(2l + 1)$重简并的*，而$"SO"(3)$群的不可约表示的维数亦为$(2l + 1)$。对于给定的$n$和$l$，这个$(2l + 1)$维的简并子空间荷载了$"SO"(3)$群的$(2l + 1)$维不可约表示$M^((l))$，即$hat(D)$为量子态转动算符)
  $
    hat(D) ket(n l m) = sum_(m' = -l)^l M_(m' m)^((l))(hat(D)) ket(n l m')
  $
  #newpara()
  对于氢原子，$V(r)$为Coulomb势，我们已经求出其解为
  $
    & hat(H) ket(n l m) = E_n ket(n l m) \
    & psi_(n l m)(vb(r)) = braket(vb(r), n l m) = R_(n l)(r) Y_(l m)(theta, phi) \
    & l = 0, 1, ... , n - 1; \
    & m = -l, -l + 1, ... , l - 1, l
  $
  能量本征值$E_n$只依赖量子数$n$，因此能级简并度为
  $
    g_n = sum_(l = 0)^(n - 1) (2l + 1) = n^2
  $
  对于给定的$n$，这个$n^2$维的本征子空间也是$"SO"(3)$群的一个$n^2$维表示的表示空间。这个表示显然是可约的，可以直和分解为
  $
    M = sum_(l = 0)^(n - 1) plus.o M^((l))
  $
]

- 问题：$n^2$维的表示$M$不是不可约表示，氢原子能级具有更高的简并度是偶然的吗？还是因为氢原子系统具有比空间转动不变性($"SO"(3)$)更高的对称性？
- 群论的语言：Hamilton量的简并子空间是荷载着其对称性群的不可约表示呢，还是也可以荷载可约表示？
- 回答：经过大量的研究，人们普遍相信，系统的能级之所以有更高的简并度，是因为系统具有更高的对称性。
- 群论的语言：Hamilton量的属于任一本征值的简并子空间，都荷载着其对称性群的一个不可约表示。
- 对于具体的研究：找到系统的对称性群，以此预言能级的简并度。如果真实的简并度比预言的要高，那么很可能意味着对称性群没有找全。

下面我们寻找比$"SO"(3)$更高的对称性群。类氢原子的Hamilton量为
$
  hat(H) = vb(hat(p))^2/(2 mu) - alpha/r, alpha = (Z e^2)/(4 pi epsilon_0)
$
我们预期它应该具有比$"SO"(3)$更大的对称性群，而更大的对称性群通常意味着需要新的守恒量作为生成元。在经典力学中，对于平方反比力场中的Kepler问题，除了能量和角动量，还有一个守恒量，即*Runge-Lenz矢量*。考虑质点在势场$- alpha/r$中的运动，运动方程为
$
  vb(dot(p)) = - alpha/r^3 vb(r)
$
两边与角动量$vb(L)$叉乘，得到
$
  vb(dot(p)) times vb(L) = - alpha/r^3 vb(r) times (vb(r) times m vb(dot(r)))
$
由角动量守恒
$
  dv(vb(L), t) = 0 => vb(dot(p)) times vb(L) = dv(, t) (vb(p) times vb(L))
$
右边通过计算得到
$
  - alpha/r^3 vb(r) times (vb(r) times m vb(dot(r))) = m alpha (vb(dot(r))/r - (dot(r) vb(r))/r^2) = m alpha dv(, t) (vb(r)/r)
$
所以
$
  dv(, t) (vb(p) times vb(L) - m alpha vb(r)/r) = 0
$
这样我们就得到了一个新的守恒量。

定义*Runge-Lenz*矢量
$
  vb(M) = vb(p) times vb(L) - m alpha vb(r)/r
$
它是一个守恒量。一个显然成立的关系是
$
  vb(M) dot vb(L) = 0
$
即$vb(M)$位于轨道平面内。计算$vb(M)$与$vb(r)$的点乘得到
$
  vb(M) dot vb(r) = 1/m (vb(p) times vb(L)) dot vb(r) - alpha r = vb(L)^2/m - alpha r
$
令$vb(M)$与$vb(r)$的夹角为$theta$，上式即为轨道方程
$
  1/r = (m alpha)/(vb(L)^2) (1 + abs(vb(M))/alpha cos theta)
$
其中$vb(M)$和$vb(L)$是守恒量，由初始条件决定。

在量子力学中，上面定义的定义的Runge-Lenz矢量并不是Hermite算符。采用*对称化排序*，将Runge-Lenz矢量定义为
$
  hat(vb(M)) = 1/(2m) (vb(hat(p)) times vb(hat(L)) - vb(hat(L)) times vb(hat(p))) - alpha hat(vb(r))/r
$
可以证明：
- 它仍是一个守恒量，即
  $
    [hat(H), hat(vb(M))] = 0
  $
- 它与角动量$vb(L)$之间仍然满足
  $
    hat(vb(M)) dot vb(hat(L)) = vb(hat(L)) dot hat(vb(M)) = 0
  $
- 直接计算可以得到
  $
    hat(vb(M))^2 = 2/m hat(H) (vb(hat(L))^2 + hbar^2) + alpha^2
  $
系统态矢的转动变换算符可以写为(沿$vu(e)_n$转动$phi$角)
$
  hat(D)_n (phi) = e^(- i/hbar phi vu(e)_n dot vb(hat(L))) = e^(- i/hbar sum_(i=1)^3 phi_i L_i)
$
这些变换构成三维空间正当转动群$"SO"(3)$，等价于$"SU"(2)$。角动量$hat(L)$的三个分量为群的*生成元*，满足*对易关系(代数)*
$
  [hat(L)_i, hat(L)_j] = i hbar epsilon_(i j k) hat(L)_k
$
引入Runge-Lenz矢量$hat(vb(M))$。直接计算或利用任意矢量算符与角动量之间的对易关系，得到
$
  [hat(M)_i, hat(L)_j] = i hbar epsilon_(i j k) hat(M)_k
$
此外，可以计算出
$
  [hat(M)_i, hat(M)_j] = - (2 i hbar)/(m) epsilon_(i j k) hat(H) hat(L)_k
$
由于出现Hamilton量算符$hat(H)$
$
  [hat(L)_i, hat(L)_j], [hat(M)_i, hat(L)_j], [hat(M)_i, hat(M)_j]
$
并不构成封闭的代数，所以无法将$hat(L)$和$hat(M)$与更大的连续对称变换的生成元对应起来。因此，我们退而求其次，只考虑能量为$E<0$ 的束缚态，即*在能量本征值为$E < 0$的不变子空间中讨论*。这样，可将Hamilton量$hat(H)$用能量本征值$E$代替。定义
$
  hat(vb(N)) = (- m/(2 E))^(1/2) hat(vb(M))
$
这样定义的$hat(vb(N))$矢量与角动量有相同的量纲。

利用$hat(vb(N))$矢量，可将对易关系写为
$
  [hat(L)_i, hat(L)_j] & = i hbar epsilon_(i j k) hat(L)_k \
  [hat(N)_i, hat(L)_j] & = i hbar epsilon_(i j k) hat(N)_k \
  [hat(N)_i, hat(N)_j] & = i hbar epsilon_(i j k) hat(L)_k
$
以上对易关系形成一个封闭的代数，对应的是四维空间正当转动群$"SO"(4)$。

由于$"SO"(n)$的生成元数目为$n(n − 1)/2$，因此$"SO"(4)$生成元为6个，即$hat(vb(L))$和$hat(vb(N))$的共6个分量。引入第四维坐标$x_4$和动量$p_4$，定义
$
  tilde(L)_(i j) = hat(x)_i hat(p)_j - hat(x)_j hat(p)_i
$
则6个生成元为
$
  tilde(L)_(2 3) = hat(L)_1, tilde(L)_(3 1) = hat(L)_2, tilde(L)_(1 2) = hat(L)_3\
  tilde(L)_(1 4) = hat(N)_1, tilde(L)_(2 4) = hat(N)_2, tilde(L)_(3 4) = hat(N)_3
$
更为简便的做法是将$hat(vb(L))$和$hat(vb(N))$线性组合，拆解为两组独立的$"SU"(2)$生成元。定义
$
  hat(vb(J))_1 = 1/2 (hat(vb(L)) + hat(vb(N)))\
  hat(vb(J))_2 = 1/2 (hat(vb(L)) - hat(vb(N)))
$
则
$
  [hat(J)_(1 i), hat(J)_(1 j)] & = i hbar epsilon_(i j k) hat(J)_(1 k), [hat(vb(J))^2_1, hat(J)_(1 i)] = 0 \
  [hat(J)_(2 i), hat(J)_(2 j)] & = i hbar epsilon_(i j k) hat(J)_(2 k), [hat(vb(J))^2_2, hat(J)_(2 i)] = 0 \
  [hat(J)_(1 i), hat(J)_(2 j)] & = 0
$
这说明$hat(vb(J))_1$ 和 $hat(vb(J))_2$ 分别独立地生成两个$"SU"(2)$群，它们都是守恒
量，即
$
  [hat(H), hat(vb(J))_1] = 0, [hat(H), hat(vb(J))_2] = 0
$
所以，体系的对称性群是两个$"SU"(2)$群的直积群$"SU"(2) times.o "SU"(2)$，与$"SO"(4)$同构。利用$hat(vb(M)) dot hat(vb(L)) = 0$
$
  hat(J)_1^2 - hat(J)_2^2 = hat(vb(L)) dot hat(vb(M)) = 0
$
同时可以计算出
$
  hat(J)_1^2 + hat(J)_2^2 = 1/2 (hat(vb(L))^2 + hat(vb(N))^2) = 1/2 (hat(vb(L))^2 - m/(2 E) hat(vb(M))^2)
$
从而
$
  hat(J)_1^2 = hat(J)_2^2 = 1/4 (hat(vb(L))^2 - m/(2 E) hat(vb(M))^2)
$
所以$hat(vb(J))_1^2$和$hat(vb(J))_2^2$必然有相同的本征值，记为$j(j + 1)hbar^2$，其中$j = 0, 1/2, 1, 3/2, ...$。

求解$hat(H), hat(vb(J))_1^2, hat(J)_(1 z), hat(vb(J))_2^2, hat(J)_(2 z)$的共同完备本征态。在能量为$E < 0$的简并子空间中，得到
$
          hat(H) ket(psi) & = E ket(psi) \
  hat(vb(J))_1^2 ket(psi) & = hat(vb(J))_2^2 ket(psi) = j(j + 1) hbar^2 ket(psi) \
    hat(J)_(1 z) ket(psi) & = m_1 hbar ket(psi) \
    hat(J)_(2 z) ket(psi) & = m_2 hbar ket(psi) \
                 m_1, m_2 & = -j, -j + 1, ... , j - 1, j
$
考虑到
$
  hat(vb(M))^2 = 2/m hat(H) (vb(hat(L))^2 + hbar^2) + alpha^2
$
带入
$
  hat(J)_1^2 = hat(J)_2^2 = 1/4 (hat(vb(L))^2 - hat(H)/E (vb(hat(L))^2 + hbar^2) - alpha^2 m/(2 E))
$
带入本征方程
$
  hat(vb(J))_1^2 ket(psi) = j(j + 1) hbar^2 ket(psi)
$
得到
$
  1/4 (hat(vb(L))^2 - hat(H)/E (vb(hat(L))^2 + hbar^2) - (alpha^2 m)/(2 E)) ket(psi) = j(j + 1) hbar^2 ket(psi)
$
注意到$hat(H) ket(psi) = E ket(psi)$，最终得到
$
  1/4 (- hbar^2 - (alpha^2 m)/(2 E)) ket(psi) = j(j + 1) hbar^2 ket(psi)
$
由此解出
$
  E = - (m alpha^2)/(2 hbar^2 (2 j + 1)^2)
$
令$n = 2j + 1$，则有
$
  E_n = - (m alpha^2)/(2 hbar^2 n^2), n = 1, 2, 3, ...
$
由于$j = 0, 1/2, 1, 3/2, ...$，所以$n$为正整数，所以结果与类氢原子能级公式一致。

*能级简并度*：若$hat(vb(J))_1^2$和 $hat(vb(J))_2^2$的本征值为$j_1(j_1 + 1)hbar^2$和$j_2(j_2 + 1)hbar^2$，则能级简并度为
$
  (2 j_1 + 1)(2 j_2 + 1)
$
另一方面，系统的对称性群$"SU"(2) times.o "SU"(2)$的不可约表示的维数亦为$(2j_1 + 1)(2j_2 + 1)$。这实际上就是两个角动量的合成。最终，由于$hat(vb(J))_1^2 = hat(vb(J))_2^2$导致$j_1 = j_2 = j$，因此能级简并度为$(2j+ 1)^2 = n^2$。

== 空间反演

=== 空间反演算符的定义与性质

在经典力学中，单粒子系统的空间反演变换定义为
$
  vb(r) -> vb(r)' = - vb(r) \
  vb(p) -> vb(p)' = - vb(p)
$
在量子力学中，我们需要定义态矢$ket(psi)$的空间反演变换
$
  ket(psi) -> ket(psi)' = hat(P) ket(psi)
$
很自然地，我们要求这个变换保持量子态的内积的模方不变。根据 Wigner 定理，$hat(P)^(-1) = hat(P)^dagger$。同时，我们要求力学量的平均值不变，即
$
  braket(psi', vb(hat(r))', psi') & = braket(psi, vb(hat(r)), psi) \
  braket(psi', vb(hat(p))', psi') & = braket(psi, vb(hat(p)), psi)
$
得到力学量的变换为
$
  hat(vb(r)) -> hat(vb(r))' = hat(P) hat(vb(r)) hat(P)^dagger\
  hat(vb(p)) -> hat(vb(p))' = hat(P) hat(vb(p)) hat(P)^dagger
$
但是我们实际上还在原地踏步，还不知道$hat(P)$的定义，更不知道它是线性的还是反线性的。
#note[
  连续变换都可以从恒等变换连续地变化而来，因此它们必然是线性的算符。而离散变换则不一定是线性的，例如空间反演和时间反演在量子力学中通常被定义为反线性算符。
]
很自然地，我们要求经典力学中的空间反演变换在平均值的意义上仍然成立，即
$
  braket(psi', hat(vb(r)), psi') & = - braket(psi, hat(vb(r)), psi) \
  braket(psi', hat(vb(p)), psi') & = - braket(psi, hat(vb(p)), psi)
$
由此得到
$
  hat(P) hat(vb(r)) hat(P)^dagger & = - hat(vb(r)) \
  hat(P) hat(vb(p)) hat(P)^dagger & = - hat(vb(p))
$
所以可以计算
$
  hat(P) [hat(x)_i, hat(p)_j] hat(P)^dagger & = [hat(P) hat(x)_i hat(P)^dagger, hat(P) hat(p)_j hat(P)^dagger] \
                                            & = [- hat(x)_i, - hat(p)_j] \
                                            & = [hat(x)_i, hat(p)_j]
$
利用$[hat(x)_i, hat(p)_j] = i hbar delta_(i j)$，得到
$
  hat(P) i hat(P)^dagger = i => hat(P) i = i hat(P)
$
$hat(P)$ 和虚数单位$i$可交换，因此$hat(P)$必为线性算符。对$hat(P) hat(vb(r)) hat(P)^dagger & = - hat(vb(r)), hat(P) hat(vb(p)) hat(P)^dagger & = - hat(vb(p))$式再做一次$hat(P)$变换，得到$hat(P)^2$与$hat(vb(r)), hat(vb(p))$都对易。总可以吸收相位因子适当定义$hat(P)$，使得 $hat(P)^2 = 1$。这样就得到
$
  hat(P)^dagger = hat(P) = hat(P)^(-1)
$
由此得到如下结论：
- $hat(P)$是Hermite算符，可以作为一个力学量， *无经典对应*。
- $hat(P)$的本征值为$±1$，对应本征态分别称为偶宇称态和奇宇称态。任意态矢必然可以展开为偶宇称态和奇宇称态的叠加
  $
    ket(psi) = 1/2 (ket(psi) + hat(P) ket(psi)) + 1/2 (ket(psi) - hat(P) ket(psi))
  $

=== 空间反演算符对基矢的作用

*空间反演算符$hat(P)$对基矢的作用如何？*
