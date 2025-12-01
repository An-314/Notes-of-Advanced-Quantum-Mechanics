#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第7次作业],
  author: "AnZrew",
  time: "2025年11月",
)

#exercise(subname: [讲义例题1])[
  $j=1/2$，利用$hat(J)_y$的矩阵形式，证明：
  $
    d^((1/2)) (phi) = mat(cos(phi/2), - sin(phi/2); sin(phi/2), cos(phi/2))
  $
]

#proof[
  $j=1/2$，利用$hat(J)_y$的矩阵形式
  $
    hat(J)_y = hbar/2 mat(0, - i; i, 0)
  $
  得到
  $
    hat(J)_y^2 = (1/4) hbar^2 mat(1, 0; 0, 1)
  $
  利用Taylor展开，可以证明
  $
    e^(- i/hbar beta hat(J)_y) & = sum_(n=0)^oo (- i/hbar beta hat(J)_y)^n/n! \
    & = sum_(n=0)^oo (- i beta/(2))^(2n)/(2n)! mat(1, 0; 0, 1) + sum_(n=0)^oo (- i beta/(2))^(2n+1)/(2n+1)! mat(0, -1; 1, 0) \
    & = mat(cos(beta/2), -sin(beta/2); sin(beta/2), cos(beta/2))\
    &= cos (beta/2) I_2 - 2 i sin (beta/2) hat(J)_y/hbar
  $
  这就是$d^((1/2))_(m'm)(phi) = braket(j m', e^(- i/hbar phi hat(J)_y), j m)$的矩阵形式。
  $
    d^((1/2)) (phi) = mat(cos(phi/2), - sin(phi/2); sin(phi/2), cos(phi/2))
  $
]

#exercise(subname: [讲义例题2])[
  $j=1$，利用$hat(J)_y$的矩阵形式，证明：
  $
    d^((1)) (phi) = mat(
      1/2(1 + cos phi), -1/(sqrt(2)) sin phi, 1/2 (1 - cos phi);
      1/(sqrt(2)) sin phi, cos phi, -1/(sqrt(2)) sin phi;
      1/2 (1 - cos phi), 1/(sqrt(2)) sin phi, 1/2 (1 + cos phi)
    )
  $
]

#proof[
  $j=1$，利用$hat(J)_y$的矩阵形式
  $
    hat(J)_y/hbar = 1/sqrt(2) mat(0, - i, 0; i, 0, - i; 0, i, 0)
  $
  得到
  $
    (hat(J)_y/ hbar)^2 = 1/2 mat(1, 0, - i; 0, 2, 0; i, 0, 1)
  $
  和
  $
    (hat(J)_y/ hbar)^3 = 1/sqrt(2) mat(0, - i, 0; i, 0, - i; 0, i, 0) = hat(J)_y/ hbar
  $
  利用Taylor展开，可以证明
  $
    e^(- i/hbar beta hat(J)_y) & = sum_(n=0)^oo (- i/hbar beta hat(J)_y)^n/n! \
    & = sum_(n=0)^oo (- i beta)^(2n)/(2n)! (hat(J)_y/ hbar)^(2n) + sum_(n=0)^oo (- i beta)^(2n+1)/(2n+1)! (hat(J)_y/ hbar)^(2n+1) \
    & = sum_(n=0)^oo (- i beta)^(2n)/(2n)! (hat(J)_y^2/hbar^2) + sum_(n=0)^oo (- i beta)^(2n+1)/(2n+1)! (hat(J)_y/hbar) \
    & = cos(beta) I_3 - i sin(beta) hat(J)_y/hbar + (1 - cos(beta)) hat(J)_y^2/hbar^2\
    & = mat(
      1/2(1 + cos beta), -1/(sqrt(2)) sin beta, 1/2 (1 - cos beta);
      1/(sqrt(2)) sin beta, cos beta, -1/(sqrt(2)) sin beta;
      1/2 (1 - cos beta), 1/(sqrt(2)) sin beta, 1/2 (1 + cos beta)
    )
  $
  这就是$d^((1))_(m'm)(phi) = braket(j m', e^(- i/hbar phi hat(J)_y), j m)$的矩阵形式。
  $
    d^((1)) (phi) = mat(
      1/2(1 + cos phi), -1/(sqrt(2)) sin phi, 1/2 (1 - cos phi);
      1/(sqrt(2)) sin phi, cos phi, -1/(sqrt(2)) sin phi;
      1/2 (1 - cos phi), 1/(sqrt(2)) sin phi, 1/2 (1 + cos phi)
    )
  $

]

#exercise(subname: [Sakurai 3.17])[
  The wave function of a particle subjected to a spherically symmetrical potential $V(r)$ is given by
  $
    psi(vb(r)) = (x + y + 3 z) f(r)
  $
  + Is $psi$ an eigenfunction of $hat(vb(L))^2$? If so, what is the $l$-value? If not, what are the possible values of $l$ that we may obtain when $hat(vb(L))^2$ is measured?
  + What are the probabilities for the particle to be found in various $m_l$ states?
  + Suppose it is known somehow that $psi(vb(x))$ is an energy eigenfunction with eigenvalue $E$. Indicate how we may find $V(r)$.
]
