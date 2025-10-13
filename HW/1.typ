#import "@preview/scripst:1.1.1": *

#show: scripst.with(
  title: [高等量子力学第1次作业],
  author: "AnZrew",
  time: "2025年9月",
)

#exercise(subname: [Steck 1.3])[
  Let $hat(Q)$ be a linear operator (not necessarily Hermitian) with eigenvalue $q$ and eigenvector $ket(q)$:
  $
    hat(Q) ket(q) = q ket(q)
  $
  Remembering that linear operators transform states on the right,
  $
    braket(phi', hat(Q), phi) := braket(phi', hat(Q) phi)
  $
  show _explicitly_ that $hat(Q)$ effectively ‘‘acts to the left’’ as
  $
    braket(q, hat(Q), phi) = q braket(q, phi) "(这个对吗？)"
  $
  $
    braket(q, hat(Q)^dagger, phi) = q^* braket(q, phi)
  $
  using _only_ the assumptions spelled out in this problem, along with the conjugation rule $braket(a, b) = braket(b, a)^*$ and the definition of the adjoint (Hermitian conjugate) of an operator.
]

#proof[
  Steck书上对$hat(Q)$的Hermite共轭的定义为
  $
    braket(psi, hat(Q) phi) = braket(hat(Q)^dagger psi, phi)
  $
  由于
  $
    braket(phi', hat(Q), phi) := braket(phi', hat(Q) phi)
  $
  结合$braket(a, b) = braket(b, a)^*$和算符的Hermite共轭的定义，我们有
  $
    braket(phi', hat(Q), phi) = braket(phi', hat(Q) phi) = braket(hat(Q)^dagger phi', phi) = braket(phi, hat(Q)^dagger phi')^* = braket(phi, hat(Q)^dagger, phi')^*
  $
  从而
  $
    braket(q, hat(Q)^dagger, phi) = braket(phi, hat(Q), q)^* = braket(phi, hat(Q) q)^* = q^* braket(phi, hat(Q))^* = q^* braket(q, phi)
  $
]

#exercise(subname: [Steck 1.8])[
  + Suppose that a Hermitian operator $hat(Q)$ has eigenvalues $q_n$ and eigenvectors $ket(q_n)$,
    $
      hat(Q) ket(q_n) = q_n ket(q_n)
    $
    for $n in ZZ^+$. If the eigenvectors $ket(q_n)$ form a complete set, prove that $hat(Q)$ may always be written in the form
    $
      hat(Q) = sum_n^oo q_n ketbra(q_n)
    $
    _Note_: to prove this equivalence you must consider the action of both expressions here on an arbitrary vector, not just an eigenvector.
  + Show that $hat(Q)$ has the same form if it is not necessarily *Hermitian*, but is normal, which means that $hat(Q) hat(Q)^dagger = hat(Q)^dagger hat(Q)$.
]

#proof[
  1. 设$ket(psi)$为Hilbert空间中的任意向量，并且Hermite算符$hat(Q)$的本征矢${ket(q_n)}$构成完备归一正交基底，因此
    $
      hat(Q) ket(psi) & = hat(Q) (sum_n ketbra(q_n)) ket(psi) = sum_n hat(Q) ket(q_n) braket(q_n, psi) \
                      & = sum_n q_n ket(q_n) braket(q_n, psi) = (sum_n q_n ketbra(q_n)) ket(psi)
    $
    由于$ket(psi)$是任意向量，因此
    $
      hat(Q) = sum_n q_n ketbra(q_n)
    $
  2. 对于正规算符$hat(Q)$，$hat(Q) hat(Q)^dagger = hat(Q)^dagger hat(Q)$，$hat(Q)$的本征方程为
    $
      hat(Q) ket(q_n) = q_n ket(q_n)
    $
    我们先证明$hat(Q)$的特征矢也是其伴随算符$hat(Q)^dagger$的特征矢：由
    $
      norm((hat(Q) - q_n hat(I))ket(q_n))^2 & = braket(q_n, (hat(Q) - q_n hat(I))^dagger (hat(Q) - q_n hat(I)), q_n) \
                                            & = braket(q_n, (hat(Q) - q_n hat(I))(hat(Q) - q_n hat(I))^dagger, q_n) \
                                            & = norm((hat(Q)^dagger - q_n^* hat(I)) ket(q_n))^2 = 0
    $
    从而有
    $
      hat(Q)^dagger ket(q_n) = q_n^* ket(q_n)
    $
    下面我们证明不同本征值对应的本征向量正交：考虑
    $
      braket(q_m, hat(Q), q_n) = braket(q_m, hat(Q) q_n) = q_n braket(q_m, q_n)\
      braket(q_m, hat(Q), q_n) = braket(hat(Q)^dagger q_m, q_n) = q_m braket(q_m, q_n)
    $
    因此
    $
      (q_n - q_m) braket(q_m, q_n) = 0
    $
    若$q_n != q_m$，则$braket(q_m, q_n) = 0$，即不同本征值对应的本征向量正交。如果我们可以从中取出一组完备归一正交基底${ket(n)}$，则
    $
      hat(Q) ket(psi) & = hat(Q) (sum_n ketbra(n)) ket(psi) = sum_n hat(Q) ket(n) braket(n, psi) = sum_n q_n ket(n) braket(n, psi) = (sum_n q_n ketbra(n)) ket(psi) \
      hat(Q) &= sum_n q_n ketbra(n)
    $
]

#exercise(subname: [Sakurai 1.4])[
  Using the rules of bra-ket algebra, prove or evaluate the following:
  + $tr(X Y) = tr(Y X)$, where X and Y are operators.
  + $(X Y)^dagger = Y^dagger X^dagger$, where X and Y are operators.
  + $exp(i f(A)) = ?$ in ket-bra form, where A is a Hermitian operator whose eigenvalues are known.
  + $sum_a' psi_a'^* (vb(x)') psi_a'(vb(x)'')$, where $psi_a' (vb(x)') = braket(vb(x'), a')$
]

#solution[
  1. 对于任意算符 $hat(X)$ 和 $hat(Y)$，取Hilbert空间的一组完备正交归一基底${ket(i)}$，则
    $
      tr(hat(X) hat(Y)) & = sum_i braket(i, hat(X) hat(Y), i) \
                        & = sum_i bra(i) hat(X) hat(I) hat(Y) ket(i) \
                        & = sum_i bra(i) hat(X) (sum_j ketbra(j)) hat(Y) ket(i) \
                        & = sum_i sum_j bra(i) hat(X) ket(j) bra(j) hat(Y) ket(i) \
                        & = sum_j sum_i bra(j) hat(Y) ket(i) bra(i) hat(X) ket(j) \
                        & = sum_j bra(j) hat(Y) (sum_i ketbra(i)) hat(X) ket(j) \
                        & = sum_j bra(j) hat(Y) hat(I) hat(X) ket(j) \
                        & = sum_j braket(j, hat(Y) hat(X), j) \
                        & = tr(hat(Y) hat(X))
    $
  2. 设$ket(alpha)$为任意右矢，并且
    $
      hat(X) hat(Y) ket(alpha) = ket(beta)
    $
    考虑与$ket(beta)$对偶的左矢$bra(beta)$，有
    $
      hat(X) hat(Y) ket(alpha) <-> bra(alpha) (hat(X) hat(Y))^dagger \
      hat(X) (hat(Y) ket(alpha)) <-> (bra(alpha) hat(Y)^dagger) hat(X)^dagger = bra(alpha) (hat(Y)^dagger hat(X)^dagger)
    $
    从而
    $
      (hat(X) hat(Y))^dagger = hat(Y)^dagger hat(X)^dagger
    $
  3. 算符$hat(A)$的本征值${a_n}$和本征矢${ket(a_n)}$满足
    $
      hat(A) ket(a_n) = a_n ket(a_n)
    $
    且${ket(a_n)}$构成完备归一正交基底，因此
    $
      hat(A) = sum_n a_n ketbra(a_n)
    $
    对于
    $
      g(hat(A)) = exp(i f(hat(A))) &=^"Taylor" sum_(k=0)^oo (g^((k))(0))/k! hat(A)^k = sum_(k=0)^oo (g^((k))(0))/k! (sum_n a_n ketbra(a_n))^k\
      &= sum_(k=0)^oo (g^((k))(0))/k! sum_n a_n^k ketbra(a_n) = sum_n (sum_(k=0)^oo (g^((k))(0))/k! a_n^k) ketbra(a_n) \
      &= sum_n g(a_n) ketbra(a_n)\
      &= sum_n exp(i f(a_n)) ketbra(a_n)
    $
  + $
      sum_a' psi_a'^* (vb(x)') psi_a'(vb(x)'') & = sum_a' braket(vb(x'), a')^* braket(vb(x''), a') \
                                               & = sum_a' braket(a', vb(x')) braket(vb(x''), a') \
                                               & = braket(vb(x''), (sum_a' ketbra(a')) vb(x')) \
                                               & = braket(vb(x''), vb(x')) \
                                               & = delta(vb(x'') - vb(x'))
    $
]
