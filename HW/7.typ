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

#exercise(subname: [Sakurai 3.17])[
  The wave function of a particle subjected to a spherically symmetrical potential $V(r)$ is given by
  $
    psi(vb(r)) = (x + y + 3 z) f(r)
  $
  + Is $psi$ an eigenfunction of $hat(vb(L))^2$? If so, what is the $l$-value? If not, what are the possible values of $l$ that we may obtain when $hat(vb(L))^2$ is measured?
  + What are the probabilities for the particle to be found in various $m_l$ states?
  + Suppose it is known somehow that $psi(vb(x))$ is an energy eigenfunction with eigenvalue $E$. Indicate how we may find $V(r)$.
]
