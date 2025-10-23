#import "@preview/scripst:1.1.1": *

#show: scripst.with(
  title: [高等量子力学第3次作业],
  author: "AnZrew",
  time: "2025年10月",
)

#exercise(subname: [Sakurai 5.24])[
  Consider a particle bound in a simple harmonic-oscillator potential. Initially ($t < 0$), it is in the ground state. At $t = 0$ a perturbation of the form
  $
    H'(x, t) = A x^2 e^(-t/tau)
  $
  is switched on. Using time-dependent perturbation theory, calculate the probability that after a sufficiently long time ($t >> tau$), the system will have made a transition to a given excited state. Consider all final states.
]

#exercise(subname: [Sakurai 5.25])[
  The unperturbed Hamiltonian of a two-state system is represented by
  $
    H_0 = mat(E_1^0, 0; 0, E_2^0)
  $
  There is, in addition, a V(t) time-dependent  perturbation
  $
    V(t) = mat(0, lambda cos omega t; lambda cos omega t, 0), lambda in RR
  $
  + At t = 0 the system is known to be in the first state, represented by
    $
      mat(1; 0)
    $
    Using time-dependent perturbation theory and assuming that $E_1^0- E_2^0$ is not   close to $±hbar omega$, derive an expression for the probability that the system is found in the second state represented by
    $
      mat(0; 1)
    $
    as a function of t(t > 0).
  + Why is this procedure not valid when $E_1^0- E_2^0$ is close to $±hbar omega$?
]

#exercise(subname: [])[
  The delta function well @hw3.1 supports a single bound state @hw3.2. Calculate the geometric phase change when $alpha$ gradually increases from $alpha_1$ to $alpha_2$. If the increase occurs at a constant rate $(dv(alpha, t)=c)$, what is the dynamic phase change for this process?
  $
    V(x) = - alpha delta(x)
  $<hw3.1>
  $
    psi(x) = sqrt(m alpha)/hbar exp(-(m alpha abs(x))/hbar^2), E = - (m alpha^2)/(2 hbar^2)
  $<hw3.2>
]
