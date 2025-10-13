#import "@preview/scripst:1.1.1": *

#show: scripst.with(
  title: [高等量子力学第2次作业],
  author: "AnZrew",
  time: "2025年10月",
)

#exercise(subname: [Steck 5.13])[
  Recall that a coherent state can be defined by $ket(alpha) = D(alpha) ket(0)$, where $ket(0)$ is the harmonic-oscillator ground state, and in the case of purely imaginary $alpha$,
  $
    D(alpha) = e^(i sqrt(2) abs(alpha) x) , i alpha in RR
  $
  Also remember that the coherent state has the expansion
  $
    ket(alpha) = sum_(n=0)^oo (alpha^n)/(sqrt(n!)) e^(-abs(alpha)^2/2) ket(n)
  $
  in the energy basis.

  Use the above coherent-state expressions to derive an expression for the $n$th moment $braket(0, x^n, 0)$ of the   harmonic-oscillator ground state.9 You should not be setting up any integrals in this problem; work   entirely in Dirac notation.
]

#exercise(subname: [Sakurai 2.7])[
  Consider a particle in three dimensions whose Hamiltonian is given by
  $
    hat(H) = hat(vb(p))^2/(2m) + hat(V)(vb(x))
  $
  By calculating $[x dot p, H]$, obtain
  $
    dv(, t) expval(vb(x) dot hat(vb(p))) = expval(hat(vb(p))^2/m - vb(x) dot grad hat(V)(vb(x)))
  $
  In order for us to identify the preceding relation with the quantum-mechanical ana­logue of the virial theorem, it is essential that the left-hand side vanish. Under what condition would this happen?
]

#exercise(subname: [Sakurai 2.12])[
  Consider a particle subject to a one-dimensional simple harmonic oscillator potential. Suppose that at $t = 0$ the state vector is given by
  $
    exp((-i p a)/hbar) ket(0)
  $
  where $p$  is the momentum operator and a is some number with dimension of length. Using the Heisenberg picture, evaluate the expectation value $expval(x)$  for $t >= 0$.
]

#exercise(subname: [Sakurai 2.13])[
  + Write down the wave function (in coordinate space) for the state specified in Problem 2.12 at $t=0$. You may use
    $
      braket(x', 0) = pi^(-1/4) x_0^(-1/2) exp(-1/2 (x'/x_0)^2) , x_0 = (hbar/(m omega))^(1/2)
    $
  + Obtain a simple expression for the probability that the state is found in the ground state at $t = 0$. Does this probability change for $t > 0$?
]
