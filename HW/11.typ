#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第11次作业],
  author: "AnZrew",
  time: "2026年1月",
)

#exercise(subname: [Steck 20.1])[
  + Consider the localized, bosonic Fock state of $N$ particles,
    $
      ket(Psi(vb(r)_1, vb(r)_2, ..., vb(r)_N)) = 1/sqrt(N!) vb(hat(psi))^(dagger)(vb(r)_1) vb(hat(psi))^(dagger)(vb(r)_2) ... vb(hat(psi))^(dagger)(vb(r)_N) ket(0)
    $
    Use the normal-mode expansion for the field operators $hat(phi)(vb(r))$, along with the bosonic commutation relations for the mode operators $a_k$ and $a_k^dagger$, to show that this Fock state is indeed an eigenstate of the total number operator:
    $
      hat(N) ket(Psi) := sum_k a_k^dagger a_k ket(Psi) = N ket(Psi)
    $
  + Indicate how the argument changes but ultimately goes through for the localized, fermionic Fock state of the same form.
]

#solution[
  + 对于Boson体系，利用场算符的模展开
    $
      hat(psi)(vb(r)) = sum_k phi_k (vb(r)) hat(a)_k,\
      hat(psi)^(dagger)(vb(r)) = sum_k phi_k^*(vb(r)) hat(a)_k^dagger
    $
    下证$[hat(N), hat(psi)^(dagger)(vb(r))] = hat(psi)^(dagger)(vb(r))$
    #proof[
      考虑到
      $
        [hat(a)_k^dagger hat(a)_k, hat(a)_k'^dagger] = [hat(a)_k^dagger, hat(a)_k'^dagger] hat(a)_k + hat(a)_k^dagger [hat(a)_k, hat(a)_k'^dagger] = delta_(k k') hat(a)_k^dagger
      $
      对$k$求和得到
      $
        [hat(N), hat(a)_k'^dagger] = sum_k [hat(a)_k^dagger hat(a)_k, hat(a)_k'^dagger] = sum_k delta_(k k') hat(a)_k^dagger = hat(a)_k'^dagger
      $
      从而对于场算符有
      $
        [hat(N), hat(psi)^(dagger)(vb(r))] & = [hat(N), sum_k phi_k^*(vb(r)) hat(a)_k^dagger] = sum_k phi_k^*(vb(r)) [hat(N), hat(a)_k^dagger] \
        & = sum_k phi_k^*(vb(r)) hat(a)_k^dagger = hat(psi)^(dagger)(vb(r))
      $
    ]
    回到原题，得到
    $
      hat(N) hat(psi)^(dagger)(vb(r)) = hat(psi)^(dagger)(vb(r)) + hat(psi)^(dagger)(vb(r)) hat(N)
    $
    从而
    $
      hat(N) ket(Psi) & = 1/sqrt(N!) hat(N) hat(psi)^(dagger)(vb(r)_1) hat(psi)^(dagger)(vb(r)_2) ... hat(psi)^(dagger)(vb(r)_N) ket(0) \
      & = 1/sqrt(N!) (hat(psi)^(dagger)(vb(r)_1) hat(N) + hat(psi)^(dagger)(vb(r)_1)) hat(psi)^(dagger)(vb(r)_2) ... hat(psi)^(dagger)(vb(r)_N) ket(0) \
      & = 1/sqrt(N!) hat(psi)^(dagger)(vb(r)_1) (hat(N) hat(psi)^(dagger)(vb(r)_2) ... hat(psi)^(dagger)(vb(r)_N) ket(0) + hat(psi)^(dagger)(vb(r)_2) ... hat(psi)^(dagger)(vb(r)_N) ket(0)) + ... + 1/sqrt(N!) hat(psi)^(dagger)(vb(r)_1) hat(psi)^(dagger)(vb(r)_2) ... hat(psi)^(dagger)(vb(r)_N) ket(0) \
      & = (N - 1)/sqrt(N!) hat(psi)^(dagger)(vb(r)_1) hat(psi)^(dagger)(vb(r)_2) ... hat(psi)^(dagger)(vb(r)_N) ket(0) + 1/sqrt(N!) hat(psi)^(dagger)(vb(r)_1) hat(psi)^(dagger)(vb(r)_2) ... hat(psi)^(dagger)(vb(r)_N) ket(0) \
      & = N ket(Psi)
    $
  + 而对于Fermion
    $
      {hat(a)_k, hat(a)_k'^dagger} = delta_(k k'), {hat(a)_k, hat(a)_k'} = 0, {hat(a)_k^dagger, hat(a)_k'^dagger} = 0\
    $
    给出
    $
      hat(a)_k hat(a)_k'^dagger = delta_(k k') - hat(a)_k'^dagger hat(a)_k, hat(a)_k'^dagger hat(a)_k^dagger = - hat(a)_k^dagger hat(a)_k'^dagger\
    $
    数算符仍定义为
    $
      hat(N) = sum_k hat(a)_k^dagger hat(a)_k
    $
    下证$[hat(N), hat(psi)^(dagger)(vb(r))] = hat(psi)^(dagger)(vb(r))$
    #proof[
      考虑到
      $
        [hat(a)_k^dagger hat(a)_k, hat(a)_k'^dagger] &= hat(a)_k^dagger hat(a)_k hat(a)_k'^dagger - hat(a)_k'^dagger hat(a)_k^dagger hat(a)_k\
        &= delta_(k k') hat(a)_k^dagger - hat(a)_k'^dagger hat(a)_k^dagger hat(a)_k + hat(a)_k'^dagger hat(a)_k^dagger hat(a)_k\
        &= delta_(k k') hat(a)_k^dagger
      $
      对$k$求和得到
      $
        [hat(N), hat(a)_k'^dagger] = sum_k [hat(a)_k^dagger hat(a)_k, hat(a)_k'^dagger] = sum_k (- delta_(k k') hat(a)_k^dagger) = hat(a)_k'^dagger
      $
      从而对于场算符有
      $
        [hat(N), hat(psi)^(dagger)(vb(r))] & = [hat(N), sum_k phi_k^*(vb(r)) hat(a)_k^dagger] = sum_k phi_k^*(vb(r)) [hat(N), hat(a)_k^dagger] \
        & = sum_k phi_k^*(vb(r)) (hat(a)_k^dagger) = hat(psi)^(dagger)(vb(r))
      $
    ]
    同上可以给出
    $
      hat(N) ket(Psi) = N ket(Psi)
    $
]
