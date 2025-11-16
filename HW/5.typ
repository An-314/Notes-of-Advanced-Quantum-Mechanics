#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第4次作业],
  author: "AnZrew",
  time: "2025年11月",
)

#exercise(subname: [Steck 7.1])[
  If we define the vector $vb(J) := vu(x) J_x + vu(y) J_y + vu(z) J_z$, where $vu(x), vu(y), vu(z)$ are the usual unit vectors, show that
  $
    vb(J) times vb(J) = i hbar vb(J)
  $
  This is really just a component-free restatement of the defining commutation relations for angular momentum.
]

#proof[
  角动量分量满足
  $
    [J_i, J_j] = i hbar epsilon_(i j k) J_k \
  $
  按分量定义算符“叉积”：
  $
    (vb(J) times vb(J))_k = epsilon_(k i j) J_i J_j \
  $
  由于$epsilon_(k i j)$反对称，乘积中的对称部分不贡献，于是
  $
    epsilon_(k i j) J_i J_j = 1/2 epsilon_(k i j) [J_i, J_j]
  $
  代入角动量对易式：
  $
    1/2 epsilon_(k i j) [J_i, J_j] = 1/2 epsilon_(k i j) i hbar epsilon_(i j l) J_l = i hbar J_k
  $
  从而
  $
    vb(J) times vb(J) = i hbar vb(J) \
  $
]

#exercise(subname: [Steck 7.2])[
  Consider an operator $hat(A)$, which satisfies $[hat(A), hat(J)_x] = [hat(A), hat(J)_y] = 0$. Show that $[hat(A), hat(J)_z] = 0$.
]

#proof[
  #lemma(subname: [Jacobi恒等式])[
    Jacobi 恒等式
    $
      [hat(A), [hat(B), hat(C)]] + [hat(B), [hat(C), hat(A)]] + [hat(C), [hat(A), hat(B)]] = 0 \
    $
  ]
  #proof[
    Jacobi 恒等式证明
    $
      [hat(A), [hat(B), hat(C)]] = hat(A) hat(B) hat(C) - hat(A) hat(C) hat(B) - hat(B) hat(C) hat(A) + hat(C) hat(B) hat(A) \
      [hat(B), [hat(C), hat(A)]] = hat(B) hat(C) hat(A) - hat(B) hat(A) hat(C) - hat(C) hat(A) hat(B) + hat(A) hat(C) hat(B) \
      [hat(C), [hat(A), hat(B)]] = hat(C) hat(A) hat(B) - hat(C) hat(B) hat(A) - hat(A) hat(B) hat(C) + hat(B) hat(A) hat(C) \
    $
    三式相加，发现各项均抵消，得证。
  ]
  回到原题，利用Jacobi恒等式
  $
    [hat(A), [hat(J)_x, hat(J)_y]] + [hat(J)_x, [hat(J)_y, hat(A)]] + [hat(J)_y, [hat(A), hat(J)_x]] = 0 \
  $
  由于$[hat(A), hat(J)_x] = [hat(A), hat(J)_y] = 0$，第二项和第三项为零，得
  $
    [hat(A), [hat(J)_x, hat(J)_y]] = 0 \
  $
  代入角动量对易式$[hat(J)_x, hat(J)_y] = i hbar hat(J)_z$，得
  $
    [hat(A), hat(J)_z] = 0 \
  $
]

#exercise(subname: [Steck 7.3])[
  Let $vb(a)$ be a vector that commutes with $vb(J)$. Then show that
  $
    [vb(a) dot vb(J),vb(J)] = - i hbar vb(a) times vb(J)
  $
]

#proof[
  对易式按分量写成
  $
    [vb(a) dot vb(J),J_k] = sum_i a_i [J_i, J_k] \
  $
  其中用到了
  $
    [a_i, J_k] = 0 \
  $
  利用角动量对易关系$[J_i, J_k] = i hbar epsilon_(i k j) J_j$，得
  $
    [vb(a) dot vb(J),J_k] = sum_i a_i i hbar epsilon_(i k j) J_j = - i hbar epsilon_(k i j) a_i J_j = - i hbar (vb(a) times vb(J))_k \
  $
  从而
  $
    [vb(a) dot vb(J),vb(J)] = - i hbar vb(a) times vb(J) \
  $
]
