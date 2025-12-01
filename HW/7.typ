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

#solution[
  + 考虑到球谐函数
    $
       Y_(1,0) & = sqrt(3/(4 pi)) cos theta = sqrt(3/(4 pi)) z/r \
       Y_(1,1) & = - sqrt(3/(8 pi)) sin theta e^(i phi) = - sqrt(3/(8 pi)) (x + i y)/r \
      Y_(1,-1) & = sqrt(3/(8 pi)) sin theta e^(- i phi) = sqrt(3/(8 pi)) (x - i y)/r
    $
    从而
    $
      z & = sqrt((4 pi)/3) r Y_(1,0) \
      x & = sqrt((2 pi)/3) r (- Y_(1,1) + Y_(1,-1)) \
      y & = i sqrt((2 pi)/3) r (Y_(1,1) + Y_(1,-1))
    $
    可以得到
    $
      psi(vb(r)) &= r (sqrt((2 pi)/3) (- Y_(1,1) + Y_(1,-1)) + i sqrt((2 pi)/3) (Y_(1,1) + Y_(1,-1)) + 3 sqrt((4 pi)/3) Y_(1,0)) f(r) \
      &= r f(r) (sqrt((2 pi)/3) (-1+i) Y_(1,1) + sqrt((2 pi)/3) (1+i) Y_(1,-1) + 3 sqrt((4 pi)/3) Y_(1,0))
    $
    $hat(vb(L))^2$的本征函数满足
    $
      hat(vb(L))^2 ket(l\, m) & = hbar^2 l (l+1) ket(l\, m) \
                   ket(l\, m) & = Y_(l\, m)(theta, phi)
    $
    也有
    $
      hat(vb(L))^2 psi(vb(r)) & = hat(vb(L))^2 (c_1 ket(1\,1) + c_(-1) ket(1\, -1) + c_0 ket(1\,0)) \
                              & = 2 hbar^2 (c_1 ket(1\,1) + c_(-1) ket(1\, -1) + c_0 ket(1\,0)) = 2 hbar^2 psi(vb(r))
    $
    从而可知$psi(vb(r))$是$hat(vb(L))^2$的本征函数，对应$l=1$。
  + 由上可知，$psi(vb(r))$在$hat(vb(L))^2$的$l=1$子空间中，可以写成
    $
      psi(vb(r)) = c_1 ket(1\,1) + c_(-1) ket(1\, -1) + c_0 ket(1\,0)
    $
    其中
    $
      c_1 & = sqrt((2 pi)/3) (-1+i) ,
            c_(-1) & = sqrt((2 pi)/3) (1+i) ,
                     c_0 & = 3 sqrt((4 pi)/3)
    $
    可以得到
    $
      abs(c_1)^2 & = (4 pi)/3 ,
                   abs(c_(-1))^2 & = (4 pi)/3 ,
                                   abs(c_0)^2 & = 12 pi
    $
    从而归一化则为概率
    $
      P(m_l=1) = ((4pi)/3)/( (4pi)/3 + (4pi)/3 + 12pi) = 1/11 \
      P(m_l=-1) = ((4pi)/3)/( (4pi)/3 + (4pi)/3 + 12pi) = 1/11 \
      P(m_l=0) = (12pi)/( (4pi)/3 + (4pi)/3 + 12pi) = 9/11
    $
  + 最后由本征方程反求$V(r)$
    $
      (- hbar^2/(2m) laplacian + V(r)) psi(vb(r)) = E psi(vb(r))
    $
    考虑到
    $
      laplacian (psi phi) = phi laplacian psi + 2 grad psi dot grad phi + psi laplacian phi
    $
    以及
    $
           grad f(r) & = f'(r) vb(r)/r \
      laplacian f(r) & = f''(r) + 2/r f'(r)
    $
    可得到
    $
      laplacian((x + y + 3 z) f(r)) & = (x + y + 3 z) laplacian f(r) + 2 grad (x + y + 3 z) dot grad f(r) + f(r) laplacian (x + y + 3 z) \
      & = (x + y + 3 z) (f''(r) + 2/r f'(r)) + 2 (vu(i) + vu(j) + 3 vu(k)) dot (f'(r) vb(r)/r)\
      & = (x + y + 3 z) (f''(r) + 4/r f'(r))
    $
    代回Schrödinger方程，得到
    $
      V(r) & = E + (hbar^2/(2m)) ((f''(r))/f(r) + 4/r (f'(r))/f(r)) \
    $
]
