#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

= 角动量理论

本章系统研究量子力学中的角动量理论。
- 角动量与空间转动有着密切的联系。在经典力学中，角动量是空间转动变换的生成元，在量子力学中亦是如此。
- 在量子力学中，角动量是算符，也是三维欧式空间中的矢量。本章将从矢量的变换性质出发导出角动量的对易关系，同时可以按照转动变换下的性质对算符进行分类。
- 求解角动量的本征值问题和角动量的合成，结果与本科量子力学并无不同。不过，可以从群论角度得到新的理解。
- 量子系统的状态由态矢来描写，量子态的转动与空间矢量的转动是不一样的。本章将研究量子态的转动，尤其是自旋态的转动。

== 空间转动与角动量

在三维Euclidean空间中，转动可以用一个$3 times 3$的*正交矩阵*来刻画。在直角坐标系中，绕着过原点的任意轴转动后，任意点的坐标$(x, y, z) = (x_1, x_2, x_3)$变为(重复指标求和，下同)
$
  x_i' = R_(i j) x_j, i = 1, 2, 3
$
其中转动矩阵$R$是正交矩阵
$
  R^TT R = 1 <=> R^(-1) = R^TT
$
三维欧氏空间中的*矢量*满足与坐标一样的变换关系。对于任意矢量$vb(A) = mat(A_1; A_2; A_3)$，其在转动后的坐标系中的表示为
$
  A_i' = R_(i j) A_j
$
显然，*转动变换保持两个矢量之间的内积不变*，尤其是保持矢量的长度不变。

需要注意的是，当我们说转动时有两种等价的观点：
- *主动观点*：矢量在转动，坐标轴 (基矢) 不变。规定逆时针转动的转角为正；
- *被动观点*：矢量固定，坐标轴 (基矢) 在转动。对转角正负的规定与主动观点相反。
本讲义中采用主动观点。

具体来说，绕$z$轴转动角度$phi$，转动矩阵为
$
  R_3(phi) = mat(
    cos phi, -sin phi, 0;
    sin phi, cos phi, 0;
    0, 0, 1
  )
$
同样，可以写出绕$x$和$y$轴转动的矩阵
$
  R_1(phi) = mat(
    1, 0, 0;
    0, cos phi, -sin phi;
    0, sin phi, cos phi
  ),
  R_2(phi) = mat(
    cos phi, 0, sin phi;
    0, 1, 0;
    -sin phi, 0, cos phi
  )
$
考虑*无穷小转动*，即$phi -> 0$。将上述三个矩阵展开到$phi$的一阶无穷小，得到
$
  R_i = 1 - i phi G_i + O(phi^2), i = 1, 2, 3
$
矩阵$G_i$即为转动矩阵的生成元，不难求出其结果为
$
  G_1 = mat(
    0, 0, 0;
    0, 0, -i;
    0, i, 0
  ),
  G_2 = mat(
    0, 0, i;
    0, 0, 0;
    -i, 0, 0
  ),
  G_3 = mat(
    0, -i, 0;
    i, 0, 0;
    0, 0, 0
  )
$
可以验证，这三个Hermite矩阵满足对易关系
$
  [G_i, G_j] = i epsilon_(i j k) G_k
$
考虑任意方向的单位矢量
$
  vu(n) = n_1 vu(e)_1 + n_2 vu(e)_2 + n_3 vu(e)_3
$
利用生成元，可将沿$vu(n)$方向的转动矩阵写为
$
  R(vu(n), phi) = e^(- i phi vu(n) dot vb(G)) = e^(- i phi (n_1 G_1 + n_2 G_2 + n_3 G_3))
$
#newpara()

上述转动是三维Euclidean空间的固有属性，所有转动矩阵${R_n(ϕ)}$构成三维正当转动群$"SO"(3)$，满足$det R = 1$。

在量子力学中，量子系统的状态由态矢$ket(psi)$来描述。当系统发生转动时，系统的态矢也将发生变化，这种变化是一种*幺正变换*。假设对系统进行转动操作$R$，态矢$ket(psi)$在转动后变为$ket(psi')$，定义转动算符$hat(D)(R)$为
$
  ket(psi') = hat(D)(R) ket(psi)
$
注意：转动操作$R$是$3 times 3$的正交矩阵，它只作用在空间矢量上。转动算符$hat(D)(R)$作用在态矢上，它的矩阵表示的维数取决于Hilbert空间的维数。转动变换是连续变换，所以必为幺正变换，
$
  hat(D)^(-1)(R) = hat(D)^(dagger)(R)
$
转动变换保持态矢的内积不变
$
  braket(psi'_1, psi'_2) = braket(psi_1, hat(D)^dagger (R) hat(D)(R), psi_2) = braket(psi_1, psi_2)
$
由于$hat(D)(R)$是连续的幺正变换，它必然可以通过某个厄米算符来生成。考虑绕着任意方向$vu(n)$的转动，转角为$phi$。对于无穷小转动，转动算符必然可以写为
$
  hat(D)(R) = 1 - i /hbar hat(J)_n phi + O(phi^2)
$
由幺正性，生成元$hat(J)_n$是Hermite算符，对应体系的某个力学量。我们将这个力学量定义为系统角动量$hat(vb(J))$沿$vu(n)$方向的分量
$
  hat(J)_n = vu(n) dot hat(vb(J)) = n_1 hat(J)_1 + n_2 hat(J)_2 + n_3 hat(J)_3
$
*有限转动变换*可以由无穷小转动连续生成，结果为
$
  hat(D)(R) = lim_(N->oo) (1 - i /hbar phi / N hat(J)_n)^N = e^(- i /hbar phi hat(J)_n )
$
对于无自旋的系统，可以验证其正确性。在转动操作下，位矢$vb(r)$变为$vb(r) = R vb(r)$，则位置本征态的变化为
$
  ket(vb(r)) -> ket(vb(r)') = ket(R vb(r)) = hat(D)(R) ket(vb(r))
$
进入坐标表象，波函数的变化为
$
  psi(vb(r)) = braket(vb(r), psi) -> psi'(vb(r)) = braket(vb(r), hat(D)(R) psi) = braket(vb(r), psi')
$
对于无穷小转动，得到
$
  R vb(r) = vb(r) + phi vu(n) times vb(r) + O(phi^2)
$
转动操作$R$的逆$R^(-1)$就是保持$vu(n)$不变，$phi -> -phi$
$
  R^(-1) vb(r) = vb(r) - phi vu(n) times vb(r) + O(phi^2)
$
再利用
$
  hat(D)^dagger (R) = hat(D)(R^(-1))
$
得到
$
  psi'(vb(r)) & = braket(vb(r), hat(D)(R), psi) = braket(R^(-1) vb(r), psi) = psi(R^(-1) vb(r)) \
              & = psi(vb(r) - phi vu(n) times vb(r) + O(phi^2)) \
              & = psi(vb(r)) - phi (vu(n) times vb(r)) dot grad psi(vb(r)) + O(phi^2) \
              & = psi(vb(r)) - phi vu(n) dot (vb(r) times (vb(r) times grad)) psi(vb(r)) + O(phi^2) \
              & = (1 - i/hbar phi vu(n) dot hat(vb(L))) psi(vb(r)) + O(phi^2)
$
其中
$
  hat(vb(L)) = - i hbar vb(r) times grad
$
就是轨道角动量算符在坐标表象的形式。

现在，我们认为$hat(vb(J))$可以是任意的角动量。不论是轨道角动量还是自旋角动量，亦或是总角动量，都是三维Euclidean空间中的矢量。但由于$hat(vb(J))$是Hilbert空间中算符，所以我们要求$hat(vb(J))$*在平均值的意义上*满足与空间矢量同样的变换关系，即：对于任意态矢$ket(psi)$，对系统进行转动操作$R$后，有
$
  braket(psi', hat(J)_i, psi') = R_(i j) braket(psi, hat(J)_j, psi)
$
#note[
  这里和上一章要求内积不变看似是矛盾的。但上章我们采用的是被动变换的观点，而这里我们采用的是主动变换的观点。
]
写成矩阵形式就是
$
  braket(psi', hat(vb(J)), psi') = R braket(psi, hat(vb(J)), psi)
$
从这样的假设出发，我们可以推导出角动量的对易关系
$
  [hat(J)_i, hat(J)_j] = i hbar epsilon_(i j k) hat(J)_k
$
考虑沿$vu(n)$方向的转动，态矢的变化为
$
  ket(psi) -> ket(psi') = hat(D)(R) ket(psi) = e^(- i/hbar phi hat(J)_n ) ket(psi)
$
转动后$hat(J)_i$的平均值为
$
  braket(psi', hat(J)_i, psi') = braket(psi, e^(i/hbar phi hat(J)_n ) hat(J)_i e^(- i/hbar phi hat(J)_n ), psi)
$
由于态矢$ket(psi)$是任意的，有
$
  e^(i/hbar phi hat(J)_n ) hat(J)_i e^(- i/hbar phi hat(J)_n ) = R_(i j) hat(J)_j
$
显然，讨论无穷小变换即可。左边采用Baker-Hausdorff公式展开，得到
$
  e^(i/hbar phi hat(J)_n ) hat(J)_i e^(- i/hbar phi hat(J)_n ) = hat(J)_i + i/hbar phi [hat(J)_n, hat(J)_i] + O(phi^2)
$
右边的一阶无穷小形式为
$
  R_(i j) hat(J)_j = (1 - i phi vb(G) dot vu(n))_(i j) hat(J)_j + O(phi^2)
$
一阶展开项系数相等，得到
$
  1/hbar sum_(k=1)^3 n_k [hat(J)_k, hat(J)_i] = - i sum_(j=1)^3 n_k (G_k)_(i j) hat(J)_j
$
三个$n_k$是独立的，所以
$
  [hat(J)_k, hat(J)_i] = - i hbar sum_(j=1)^3 (G_k)_(i j) hat(J)_j = i hbar epsilon_(k i j) hat(J)_j
$
将下标改写一下即得到熟悉的对易关系
$
  [hat(J)_i, hat(J)_j] = i hbar epsilon_(i j k) hat(J)_k
$
可见，角动量的空间矢量属性决定了它的对易关系。我们可以反向检验一下。考虑沿$z$轴的转动，利用Baker-Hausdorff公式展开
$
  e^(i/hbar phi hat(J)_z) hat(J)_x e^(- i/hbar phi hat(J)_z ) = hat(J)_x + (i phi)/hbar [hat(J)_z, hat(J)_x] + 1/2! ((i phi) / hbar)^2 [hat(J)_z, [hat(J)_z, hat(J)_x]] + \ 1/3! ((i phi) / hbar)^3 [hat(J)_z, [hat(J)_z, [hat(J)_z, hat(J)_x]]] + ... \
$
利用角动量的对易关系，进一步得到
$
  e^(i/hbar phi hat(J)_z) hat(J)_x e^(- i/hbar phi hat(J)_z ) & = hat(J)_x (1 - (phi^2)/2! + (phi^4)/4! - ...) - hat(J)_y (phi - (phi^3)/3! + (phi^5)/5! - ...) \
  & = hat(J)_x cos phi - hat(J)_y sin phi
$
所以
$
  braket(psi', hat(J)_x, psi') = braket(psi, e^(i/hbar phi hat(J)_z) hat(J)_x e^(- i/hbar phi hat(J)_z ), psi) = braket(psi, hat(J)_x, psi) cos phi - braket(psi, hat(J)_y, psi) sin phi
$
类似地，可以计算
$
  braket(psi', hat(J)_y, psi') = braket(psi, hat(J)_x, psi) sin phi + braket(psi, hat(J)_y, psi) cos phi\
  braket(psi', hat(J)_z, psi') = braket(psi, hat(J)_z, psi)
$
写成矩阵形式，即
$
  braket(psi', hat(vb(J)), psi') = R_3(phi) braket(psi, hat(vb(J)), psi)
$

#newpara()

顺着这个思路，我们可以很自然地将量子力学中的算符按照其在转动变换下的性质进行分类。
- 标量算符：算符$hat(S)$的平均值在转动变换下与标量一致，即
  $
    braket(psi', hat(S), psi') = braket(psi, hat(S), psi)\
    => e^(i/hbar phi hat(J)_n) hat(S) e^(- i/hbar phi hat(J)_n) = hat(S) => [hat(J)_n, hat(S)] = 0
  $
- 矢量算符：算符$hat(vb(V)) = (hat(V)_1, hat(V)_2, hat(V)_3)$的平均值在转动变换下与矢量一致，即
  $
    braket(psi', hat(V)_i, psi') = R_(i j) braket(psi, hat(V)_j, psi)\
    braket(psi', hat(vb(V)), psi') = R braket(psi, hat(vb(V)), psi)\
    => e^(i/hbar phi hat(J)_n) hat(V)_i e^(- i/hbar phi hat(J)_n) = R_(i j) hat(V)_j\
    e^(i/hbar phi hat(J)_n) hat(vb(V)) e^(- i/hbar phi hat(J)_n) = R hat(vb(V))
  $
  对应的对易关系为
  $
    [hat(J)_i, hat(V)_j] = i hbar epsilon_(i j k) hat(V)_k
  $
- 张量算符：算符$hat(T)_(i j)$的平均值在转动变换下与张量一致，即
  $
    braket(psi', hat(T)_(i j), psi') = R_(i k) R_(j l) braket(psi, hat(T)_(k l), psi)\
    => e^(i/hbar phi hat(J)_n) hat(T)_(i j) e^(- i/hbar phi hat(J)_n) = R_(i k) R_(j l) hat(T)_(k l)
  $
  对应的对易关系为
  $
    [hat(J)_i, hat(T)_(j k)] = i hbar (epsilon_(i j l) hat(T)_(l k) + epsilon_(i k l) hat(T)_(j l))
  $

另一种引入角动量对易关系的方式：一些教材中利用*群的同态*关系引入角动量的对易关系。

可以证明，转动算符满足如下条件：
- 单位元：$R dot 1 = R => hat(D)(R) hat(D)(1) = hat(D)(R) => hat(D)(1) = 1$
- 封闭性：$R_a R_b = R_c => hat(D)(R_a) hat(D)(R_b) = hat(D)(R_c)$
- 逆元素：$R R^(-1) = 1 => hat(D)(R) hat(D)(R^(-1)) = hat(D)(1) = 1 => hat(D)(R^(-1)) = hat(D)^(-1)(R)$
- 结合律：显然满足。$(R_a R_b) R_c = R_a (R_b R_c) => (hat(D)(R_a) hat(D)(R_b)) hat(D)(R_c) = hat(D)(R_a) (hat(D)(R_b) hat(D)(R_c))$
所有转动算符$hat(D)(R)$也构成群，与转动矩阵${R}$构成的群同态。因此，转动算符$hat(D)(R)$保持转动矩阵$R$之间的乘法关系。

考虑无穷小转动，转动矩阵可以展开为
$
  R_i (phi) = 1 - i phi G_i - 1/2 phi^2 G_i^2 + O(phi^3)\
  R_i^(-1) (phi) = 1 + i phi G_i - 1/2 phi^2 G_i^2 + O(phi^3)
$
考虑转动矩阵群中四个元素的乘积$R^(-1)_x (phi) R^(-1)_y (theta) R_x (phi) R_y (theta)$，利用群的封闭性，有
$
  &R^(-1)_x (phi) R^(-1)_y (theta) R_x (phi) R_y (theta)\
  =& (1 + i phi G_x - 1/2 phi^2 G_x^2)(1 + i theta G_y - 1/2 theta^2 G_y^2)(1 - i phi G_x - 1/2 phi^2 G_x^2)(1 - i theta G_y - 1/2 theta^2 G_y^2)\
  =& 1 - i phi theta [G_x, G_y] + ...\
  =& 1 - i phi theta G_z + ...\
  =& R_z (phi theta)
$
群同态保持乘法关系，得到
$
  hat(D)_x^(-1)(phi) hat(D)_y^(-1)(theta) hat(D)_x (phi) hat(D)_y (theta) = hat(D)_z (phi theta)
$
利用转动算符的展开式
$
  hat(D)_i (phi) = 1 - (i phi)/(hbar) hat(J)_i - (phi^2)/(2hbar^2) hat(J)_i^2 + O(phi^3)\
  hat(D)_i^(-1) (phi) = 1 + (i phi)/(hbar) hat(J)_i - (phi^2)/(2hbar^2) hat(J)_i^2 + O(phi^3)
$
带入计算，保留到二阶无穷小项，得到
$
  1 - (phi theta)/(hbar^2) [hat(J)_x, hat(J)_y] + ... = 1 - (i phi theta)/(hbar) hat(J)_z + ...\
  => [hat(J)_x, hat(J)_y] = i hbar hat(J)_z
$
同理可得其他两个对易关系
$
  [hat(J)_y, hat(J)_z] = i hbar hat(J)_x,\
  [hat(J)_z, hat(J)_x] = i hbar hat(J)_y
$
Sakurai书上的做法与此处有些不同。Sakurai书上计算的是群元的对易关系。对于转动操作$R$，计算得到
$
  R_x (epsilon) R_y (epsilon) - R_y (epsilon) R_x (epsilon) = R_z (epsilon^2)
$
然后认为群同态保持群元间的对易关系，从而有
$
  hat(D)_x (epsilon) hat(D)_y (epsilon) - hat(D)_y (epsilon) hat(D)_x (epsilon) = hat(D)_z (epsilon^2)
$
带入转动算符展开式得到$[hat(J)_x, hat(J)_y] = i hbar hat(J)_z$等关系。群同态保持乘法关系是没问题的，但是保持对易关系，则需要Lie群Lie代数的知识去证明。

== 角动量的本征值

=== 角动量的平方、升降算符以及共同本征态

由于角动量的不同分量之间不对易，因此不可能有*共同本征态*。定义*角动量的平方*
$
  hat(vb(J))^2 = hat(J)_x^2 + hat(J)_y^2 + hat(J)_z^2 = sum_(i=1)^3 hat(J)_i^2
$
可以证明它与角动量的任意分量对易
$
  [hat(vb(J))^2, hat(J)_i] & = sum_i (hat(J)_j [hat(J)_j, hat(J)_i] + [hat(J)_j, hat(J)_i] hat(J)_j) \
                           & = i hbar sum_(j, k) epsilon_(j i k) (hat(J)_j hat(J)_k + hat(J)_k hat(J)_j) \
                           & = i hbar sum_(j, k) (epsilon_(j i k) hat(J)_j hat(J)_k + epsilon_(k j i) hat(J)_k hat(J)_j) \
                           & = 0
$
所以，$hat(vb(J))^2$与任意分量$hat(J)_i$有*共同本征态*。不失一般性，求解$hat(vb(J))^2$ 与$hat(J)_z$的共同本征态和对应的本征值。考虑到$hat(vb(J))^2$与$hat(J)_z$的量纲分别为$[hbar^2]$ 和 $[hbar]$，将本征方程写为
$
  hat(vb(J))^2 ket(lambda\, m) = lambda hbar^2 ket(lambda\, m)\
  hat(J)_z ket(lambda\, m) = m hbar ket(lambda\, m)
$
其中$lambda$和$m$是标记共同本征态的两个无量纲量子数。结果为
$
  lambda = j(j+1), m = -j, -j+1, ..., j-1, j\
  m = -j, -j+1, ..., j-1, j
$
对于给定的$j$，$2j + 1$个简并的本征态$ket(j\, m)$就构成三维正当转动群$"SO"(3)$的$2j + 1$维表示空间。下面采用类似求解谐振子本征值问题的代数方法(*升降算符方法*)来求解。

引入两个新的算符
$
  hat(J)_+ = hat(J)_x + i hat(J)_y,\
  hat(J)_- = hat(J)_x - i hat(J)_y = hat(J)_+^dagger
$
则
$
  hat(vb(J))^2 & = hat(J)^2_z + 1/2 (hat(J)_+ hat(J)_- + hat(J)_- hat(J)_+) \
               & = hat(J)^2_z + 1/2 (hat(J)_+ hat(J)_+^dagger + hat(J)_+^dagger hat(J)_+)
$
易证如下对易关系
$
  [hat(J)_+, hat(J)_-] = 2 hbar hat(J)_z, [hat(J)_z, hat(J)_plus.minus] = plus.minus hbar hat(J)_plus.minus, [hat(vb(J))^2, hat(J)_plus.minus] = [hat(vb(J))^2, hat(J)_z] = 0
$
考虑$hat(vb(J))^2 - hat(J)_z^2$在任意本征态中的平均值，得到
$
  braket(j\, m, hat(vb(J))^2 - hat(J)_z^2, j\, m) & = braket(j\, m, 1/2 (hat(J)_+ hat(J)_- + hat(J)_- hat(J)_+), j\, m) \
  & = 1/2 (braket(j\, m, hat(J)_-^dagger hat(J)_-, j\, m) + braket(j\, m, hat(J)_+^dagger hat(J)_+, j\, m)) \
  & = 1/2 (norm(hat(J)_- ket(j\, m))^2 + norm(hat(J)_+ ket(j\, m))^2) >= 0
$
而
$
  braket(j\, m, hat(vb(J))^2 - hat(J)_z^2, j\, m) = (lambda - m^2) hbar^2 braket(j\, m, j\, m)
$
所以
$
  lambda - m^2 >= 0 => lambda >= m^2
$
利用对易关系得到
$
  hat(vb(J))^2 hat(J)_+ ket(lambda\, m) & = hat(J)_+ hat(vb(J))^2 ket(lambda\, m) = lambda hbar^2 hat(J)_+ ket(lambda\, m)\
  hat(J)_z hat(J)_+ ket(lambda\, m) & = (hat(J)_+ hat(J)_z + hbar hat(J)_+) ket(lambda\, m) = (m + 1) hbar hat(J)_+ ket(lambda\, m)
$
以及
$
  hat(vb(J))^2 hat(J)_- ket(lambda\, m) & = hat(J)_- hat(vb(J))^2 ket(lambda\, m) = lambda hbar^2 hat(J)_- ket(lambda\, m)\
  hat(J)_z hat(J)_- ket(lambda\, m) & = (hat(J)_- hat(J)_z - hbar hat(J)_-) ket(lambda\, m) = (m - 1) hbar hat(J)_- ket(lambda\, m)
$
所以，若$ket(lambda\, m)$是$hat(vb(J))^2$和$hat(J)_z$的共同本征态，则$hat(J)_+ ket(lambda\, m)$和$hat(J)_- ket(lambda\, m)$也是它们的共同本征态，对应的本征值分别为$(lambda, m + 1)$和$(lambda, m - 1)$。这些本征值和本征态的关系为：
#align(center)[#three-line-table[
  | 共同本征态 | $hat(J)_z$的本征值 | $hat(vb(J))^2$的本征值 |
  | ------------ | ------------------ | --------------------- |
  | ... | ... | ... |
  | $hat(J)_+^2 ket(lambda\, m)$ | $(m + 2) hbar$ | $lambda hbar^2$ |
  | $hat(J)_+ ket(lambda\, m)$ | $(m + 1) hbar$ | $lambda hbar^2$ |
  | $ket(lambda\, m)$ | $m hbar$ | $lambda hbar^2$ |
  | $hat(J)_- ket(lambda\, m)$ | $(m - 1) hbar$ | $lambda hbar^2$ |
  | $hat(J)_-^2 ket(lambda\, m)$ | $(m - 2) hbar$ | $lambda hbar^2$ |
  | ... | ... | ... |
]]
与谐振子代数解法类比，称$hat(J)_+$为*上升算符*，$hat(J)_-$为*下降算符*。它们的作用是：作用在共同本征态$ket(lambda\, m)$上，使得$hat(J)_z$的本征值增加或减少$hbar$，而$hat(vb(J))^2$的本征值不变。因为$m^2 <= lambda$，所以存在最大值$b$和最小值$a$，因此不妨设
$
  m = a, a + 1, ..., b - 1, b
$
对于最大值$b$对应的本征态$ket(lambda\, b)$，有
$
  hat(J)_z hat(J)_+ ket(lambda\, b) = (b + 1) hbar hat(J)_+ ket(lambda\, b)
$
与$b$是最大值矛盾，所以必然有
$
  hat(J)_+ ket(lambda\, b) = 0
$
利用
$
  0 = hat(J)_- hat(J)_+ ket(lambda\, b) = (hat(vb(J))^2 - hat(J)_z^2 - hbar hat(J)_z) ket(lambda\, b) = (lambda - b^2 - b) hbar^2 ket(lambda\, b)
$
得到
$
  lambda = b(b + 1)
$
对于最小值$a$对应的本征态$ket(lambda\, a)$，有
$
  hat(J)_z hat(J)_- ket(lambda\, a) = (a - 1) hbar hat(J)_- ket(lambda\, a)
$
与$a$是最小值矛盾，所以必然有
$
  hat(J)_- ket(lambda\, a) = 0
$
利用
$
  0 = hat(J)_+ hat(J)_- ket(lambda\, a) = (hat(vb(J))^2 - hat(J)_z^2 + hbar hat(J)_z) ket(lambda\, a) = (lambda - a^2 + a) hbar^2 ket(lambda\, a)
$
得到
$
  lambda = a(a - 1)
$
由于$lambda$的结果是唯一的，所以
$
  a(a - 1) = b(b + 1) => (a + b)(a - b - 1) = 0
$
由于$a$为最小值，所以只能有
$
  a + b = 0 => b = -a
$
故$m$的取值为
$
  m = -b, -b + 1, ..., b - 1, b
$
最后的问题是求解$b$的取值。由于$lambda = b(b + 1)$，因此可以将本征态$ket(lambda\, m)$记为$ket(b\, m)$。从最小$m$取值的本征态$ket(b\, -b)$出发，依次作用上升算符，最后必然得到最大$m$取值的本征态$ket(b\, b)$，即
$
  (hat(J)_+)^(2b) ket(b\, -b) = ket(b\, b)
$
作用的次数$2b$的取值必然为$2b = 0, 1, 2, ...$，按照习惯，令$b = j$，则
$
  j & = 0, 1/2, 1, 3/2, 2, ... \
  m & = -j, -j + 1, ..., j - 1, j
$
$m$的可取值为$2j + 1$个。可以看到，$j$可以取$1/2, 3/2, ...$这样的半奇数，必然没有轨道角动量这样的经典对应。

=== 角动量的矩阵形式

$hat(vb(J))^2$与$hat(J)_z$ 的共同本征态${ket(j\, m)}$ 构成一个表象(角动量表象)。在这个表象中，角动量的各个分量呈现为矩阵形式。首先，$hat(vb(J))^2$与$hat(J)_z$显然是对角的，即
$
  braket(j'\, m', hat(vb(J))^2, j\, m) & = j(j + 1) hbar^2 delta_(j' j) delta_(m' m) \
      braket(j'\, m', hat(J)_z, j\, m) & = m hbar delta_(j' j) delta_(m' m)
$
接下来计算$hat(J)_x$和$hat(J)_y$ 的矩阵形式。根据上升算符$hat(J)_+$和下降算符
$hat(J)_-$的性质，可知
$
  hat(J)_+ ket(j\, m) = a_(j m) ket(j\, m + 1) => bra(j\, m) a^*_(j m) = bra(j\, m + 1) hat(J)_- \
  hat(J)_- ket(j\, m) = b_(j m) ket(j\, m - 1) => bra(j\, m) b^*_(j m) = bra(j\, m - 1) hat(J)_+
$
因此，可以计算内积
$
  braket(j\, m, hat(J)_- hat(J)_+, j\, m) & = braket(j\, m, hat(vb(J))^2 - hat(J)_z^2 - hbar hat(J)_z, j\, m) \
                                          & = (j(j + 1) - m(m + 1)) hbar^2 \
$
另一方面
$
  braket(j\, m, hat(J)_- hat(J)_+, j\, m) & = braket(j\, m+1, a^*_(j m) a_(j m) j\, m + 1) \
                                          & = abs(a_(j m))^2
$
选取适当的相位，可以约定
$
  a_(j m) = sqrt(j(j + 1) - m(m + 1)) hbar = sqrt((j - m)(j + m + 1)) hbar
$
同样可以得到
$
  b_(j m) = sqrt(j(j + 1) - m(m - 1)) hbar = sqrt((j + m)(j - m + 1)) hbar
$
上述结果总结为
$
  hat(J)_+ ket(j\, m) = sqrt((j - m)(j + m + 1)) hbar ket(j\, m + 1)\
  hat(J)_- ket(j\, m) = sqrt((j + m)(j - m + 1)) hbar ket(j\, m - 1)
$
所以，$hat(J)_+$和$hat(J)_-$的矩阵元为
$
  braket(j'\, m', hat(J)_+, j\, m) & = sqrt((j - m)(j + m + 1)) hbar delta_(j' j) delta_(m' m + 1) \
  braket(j'\, m', hat(J)_-, j\, m) & = sqrt((j + m)(j - m + 1)) hbar delta_(j' j) delta_(m' m - 1)
$
利用$hat(J)_plus.minus$与$hat(J)_x, hat(J)_y$的关系，可以得到
$
  braket(j'\, m', hat(J)_x, j\, m) & = 1/2 (braket(j'\, m', hat(J)_+ + hat(J)_-, j\, m)) \
  & = 1/2 (sqrt((j - m)(j + m + 1)) hbar delta_(j' j) delta_(m' m + 1) + sqrt((j + m)(j - m + 1)) hbar delta_(j' j) delta_(m' m - 1)) \
  braket(j'\, m', hat(J)_y, j\, m) & = 1/(2i) (braket(j'\, m', hat(J)_+ - hat(J)_-, j\, m)) \
  & = 1/(2i) (sqrt((j - m)(j + m + 1)) hbar delta_(j' j) delta_(m' m + 1) - sqrt((j + m)(j - m + 1)) hbar delta_(j' j) delta_(m' m - 1))
$
可以看到，这些矩阵对于$j$量子数都是对角的，也就是说这些矩阵是*分块对角*矩阵，每个分块对应一个$2j + 1$维的*不变子空间*。

#example(subname: [$j = 1/2$ (例如电子自旋)])[
  由于$j = 1/2$，则$m$的取值为$1/2, -1/2$，这是一个二维Hilbert空间，基矢为$ket(1/2\, 1/2)$和$ket(1/2\, -1/2)$，即$ket(arrow.t)$和$ket(arrow.b)$。按此基矢排列顺序，利用前面的矩阵元公式，得到
  $
    hat(vb(J))^2 = 3/4 hbar^2 mat(1, 0; 0, 1) = 3/4 hbar^2 I_2,\
    hat(J)_z = (1/2) hbar mat(1, 0; 0, -1) = 1/2 hbar sigma_z,\
    hat(J)_x = (1/2) hbar mat(0, 1; 1, 0) = 1/2 hbar sigma_x,\
    hat(J)_y = (1/2) hbar mat(0, -i; i, 0) = 1/2 hbar sigma_y,\
    hat(J)_+ = hbar mat(0, 1; 0, 0), = hbar sigma_+,\
    hat(J)_- = hbar mat(0, 0; 1, 0) = hbar sigma_-
  $
]

#example(subname: [$j = 1$ (例如光子自旋)])[
  由于$j = 1$，则$m$的取值为$1, 0, -1$，这是一个Hilbert维希尔伯特空间，基矢为$ket(1\, 1)$，$ket(1\, 0)$和$ket(1\, -1)$。按此基矢排列顺序，利用前面的矩阵元公式，得到
  $
    hat(vb(J))^2 = 2 hbar^2 mat(1, 0, 0; 0, 1, 0; 0, 0, 1),
    hat(J)_z = hbar mat(1, 0, 0; 0, 0, 0; 0, 0, -1),\
    hat(J)_x = (hbar)/(sqrt(2)) mat(0, 1, 0; 1, 0, 1; 0, 1, 0),
    hat(J)_y = (hbar)/(sqrt(2)) mat(0, -i, 0; i, 0, -i; 0, i, 0),\
  $
]
#newpara()

问题：为什么粒子的自旋量子数$j$是固定的？

对于轨道角动量，根据量子力学的态叠加原理，粒子的状态可以是$l$取不同值的本征态$ket(l\, m)$的叠加态，即存在这样的态
$
  ket(psi) = sum_(l=0)^oo sum_(m=-l)^l c_(l m) ket(l\, m)
$
但是，对于自旋角动量，每种粒子的自旋量子数$j$的取值是固定的。例如，电子的$j = 1/2$，光子的$j = 1$。你找不到一种粒子，它处于不同$j$的叠加态！

这个问题，量子力学本身回答不了。若想解开这个谜团，需要进一步学习量子场论等课程。

== 量子态的转动

对系统进行转动操作：将系统沿任意方向$vu(n)$转动角度$phi$，对应的转动矩阵为
$
  R_vu(n) ( phi) = e^(- i phi vb(G) dot vu(n)) = e^(- i phi (G_1 n_1 + G_2 n_2 + G_3 n_3))
$
它只作用在空间矢量上。对应的量子态的转动算符为
$
  hat(D)(R) = hat(D)(R_vu(n) ( phi)) = e^(- i/hbar phi hat(J) dot vu(n)) = e^(- i/hbar phi hat(J)_n)
$
由于$hat(vb(J))^2$和$hat(J)_n$对易，因此$hat(vb(J))^2$和$hat(D)(R)$也对易，所以有
$
  hat(vb(J))^2 hat(D)(R) ket(j\, m) = hat(D)(R) hat(vb(J))^2 ket(j\, m) = j(j + 1) hbar^2 hat(D)(R) ket(j\, m)
$
这说明$hat(D)(R) ket(j\, m)$仍然是$hat(vb(J))^2$的本征态，对应的本征值不变。对于确定的$j, hat(D) (R)ket(j\, m)$ 必然也在$2j+1$个正交归一基矢$ket(j\, m)$张成的线性空间中，即：这个$2j + 1$维的子空间是转动变换的*不变子空间*。所以必然有
$
  hat(D)(R) ket(j\, m) = sum_(m'=-j)^j cal(D)_(m' m)^((j)) (R) ket(j\, m')
$
叠加系数$cal(D)_(m' m)^((j)) (R)$是什么呢？将上式两边与左矢$bra(j\, m'')$做内积，得到
$
  cal(D)_(m_1 m_2)^((j)) (R) = braket(j\, m_1, hat(D)(R), j\, m_2)\
  m_1, m_2 = -j, -j + 1, ..., j - 1, j
$
从而
$
  hat(D)(R) ket(j\, m) = sum_(m'=-j)^j ket(j\, m') braket(j\, m', hat(D)(R), j\, m)
$
这说明，在$(2j + 1)$维的不变子空间中，有*完备性关系*
$
  hat(D)(R) = sum_(m'-j)^j ketbra(j\, m')
$
叠加系数$cal(D)_(m' m)^((j)) (R)$ (又称为 Wigner 函数) 构成的矩阵$cal(D)^((j))(R)$ 称为量子态的转动矩阵。以上讨论说明，研究量子态的转动，归结为计算转动矩阵$cal(D)^((j))(R)$。
