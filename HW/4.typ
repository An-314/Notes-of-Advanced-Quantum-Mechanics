#import "@preview/scripst:1.1.1": *

#show: scripst.with(
  title: [高等量子力学第4次作业],
  author: "AnZrew",
  time: "2025年10月",
)

#exercise(subname: [Steck 10.7])[
  Show that under Hamiltonian evolution,
  $
    pdv(, t) rho = (- i)/hbar [hat(H), rho]
  $
  the purity $Tr(rho^2)$ is a constant of the motion.
]

#solution[
  纯度
  $
    dv(, t) Tr(rho^2) = Tr(dv(, t) (rho^2)) = 2 Tr(rho dv(rho, t)) \
  $
  代入 von Neumann 方程
  $
    dv(rho, t) = (- i)/hbar [hat(H), rho] \
  $
  得
  $
    dv(, t) Tr(rho^2) = (- 2 i)/hbar Tr(rho [hat(H), rho]) = (- 2 i)/hbar (Tr(rho hat(H) rho) - Tr(rho rho hat(H))) = 0 \
  $
  由于迹的循环不变性，故纯度是守恒量。
]

#exercise(subname: [Sakurai 3.10])[
  + Consider a pure ensemble of identically prepared spin $1/2$ systems. Suppose the expectation values $expval(S_x)$ and $expval(S_z)$ and the sign of $expval(S_y)$ are known. Show how we may determine the state vector. Why is it unnecessary to know the magnitude of $expval(S_y)$?
  + Consider a mixed ensemble of spin $1/2$ systems. Suppose the ensemble averages $expval(expval(S_x))$, $expval(expval(S_y))$, and $expval(expval(S_z))$ are all known. Show how we may construct the $2 times 2$ density matrix that characterizes the ensemble.
]

#solution[
  + 用 Bloch 矢量表示纯态
    $
      vb(P) = mat(P_x; P_y; P_z), abs(vb(P)) = 1 \
    $
    $S_i$是自旋算符，满足
    $
      S_i = (hbar/2) sigma_i \
    $
    其中$sigma_i (i = x; y; z)$是Pauli矩阵。Bloch 矢量的分量与期望值的关系为
    $
      P_i = expval(sigma_i) = 2/hbar expval(S_i) \
    $
    因为$abs(vb(P)) = 1$，所以只需要知道$P_x$，$P_z$和$sgn(P_y)$就可以确定$P_y$，从而确定Bloch矢量。

    态矢量可以表示为
    $
      ket(psi) = alpha ket(arrow.t) + beta ket(arrow.b) \
    $
    其中
    $
      P_x & = alpha^* beta + alpha beta^* = sin theta cos phi \
      P_y & = - i (alpha^* beta - alpha beta^*) = sin theta sin phi \
      P_z & = abs(alpha)^2 - abs(beta)^2 = cos theta \
    $
    从而
    $
      abs(alpha)^2 + abs(beta)^2 = 1, abs(alpha)^2 - abs(beta)^2 = cos theta \
      => abs(alpha) = sqrt((1 + cos theta)/2) = cos theta/2, abs(beta) = sqrt((1 - cos theta)/2) = sin theta/2 \
    $
    $alpha, beta$差一个全局相位，可以取为实数，从而
    $
      alpha = cos theta/2, beta = e^(i phi) sin theta/2 \
    $
    态矢量为
    $
      ket(psi) = cos theta/2 ket(arrow.t) + e^(i phi) sin theta/2 ket(arrow.b) \
    $
    而$phi$由$sgn(P_y)$决定。

  + 对任意自旋$1/2$系统，类似地可以得到
    $
      rho = 1/2 (I + vb(P) dot vb(sigma)) , vb(P) = mat(P_x; P_y; P_z) , P_i = expval(expval(sigma_i)) = 2/hbar expval(expval(S_i)) \
    $
    从而在$ket(arrow.t), ket(arrow.b)$基矢下，密度矩阵为
    $
      rho &= 1/2 mat(1 + P_z, P_x - i P_y; P_x + i P_y, 1 - P_z) \
      &= 1/2 mat(1 + 2/hbar expval(expval(S_z)), 2/hbar expval(expval(S_x)) - (2i)/hbar expval(expval(S_y)); 2/hbar expval(expval(S_x)) + (2i)/hbar expval(expval(S_y)), 1 - 2/hbar expval(expval(S_z))) \
    $

]
