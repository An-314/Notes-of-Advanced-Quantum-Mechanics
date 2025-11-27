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
