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
