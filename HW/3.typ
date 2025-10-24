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

#solution[
  微扰
  $
    hat(H) = hat(H)_0 + hat(V)(t)\
    hat(H)_0 = p^2/(2m) + (1/2) m omega^2 x^2\
    hat(V)(t) = A x^2 e^(-t/tau) Theta(t)
  $
  在相互作用绘景中
  $
    i hbar pdv(, t) ket(psi^"I" (t)) = hat(V)^"I" (t) ket(psi^"I" (t))\
    hat(V)^"I" (t) = e^(i/hbar hat(H)_0 t) hat(V)(t) e^(-i/hbar hat(H)_0 t)
  $
  初态为$ket(0)$，一阶微扰下转移到终态$ket(n)$的概率振幅为
  $
    cal(A)_"n0"^((1)) (t) = (- i)/hbar integral_0^t dd(t_1) e^(i/hbar (E_n - E_0) t_1) braket(n, hat(V)(t_1), 0)\
  $
  计算矩阵元
  $
    cal(A)_"n0"^((1)) (t->oo) &= (- i A)/hbar integral_0^oo dd(t_1) e^(i/hbar (E_n - E_0) t_1) e^(-t_1/tau) braket(n, x^2, 0)\
    &= - (i A)/hbar braket(n, x^2, 0) tau/(1 - i/hbar (E_n - E_0) tau) \
    &= - (i A)/hbar braket(n, x^2, 0) tau/(1 - i omega_(n 0) tau)
  $
  其中
  $
    omega_(n 0) = (E_n - E_0)/hbar = n omega
  $
  故一阶跃迁概率
  $
    P_(0->n)^((1)) = abs(cal(A)_"n0"^((1)) (t->oo))^2 = (A^2 tau^2)/(hbar^2) abs(braket(n, x^2, 0))^2 / (1 + ((E_n - E_0)^2 tau^2)/hbar^2) = (A^2 tau^2)/(hbar^2) abs(braket(n, x^2, 0))^2 / (1 + n^2 omega^2 tau^2)
  $
  而
  $
    hat(x) = sqrt(hbar/(2 m omega)) (hat(a) + hat(a)^dagger)\
    hat(x)^2 = (hbar/(2 m omega)) (hat(a)^2 + (hat(a)^dagger)^2 + hat(a) hat(a)^dagger + hat(a)^dagger hat(a))\
  $
  作用在基态上
  $
    hat(x)^2 ket(0) = (sqrt(2) hbar)/(2 m omega) ket(2) + hbar/(2 m omega) ket(0) \
  $
  从而有
  $
    braket(n, x^2, 0) = (sqrt(2) hbar)/(2 m omega) delta_"n2" + (hbar)/(2 m omega) delta_"n0" \
  $
  对所有激发态而言
  $
    P_(0->2)^((1)) = (A^2 tau^2)/hbar^2 (hbar^2/(2 m^2 omega^2))/(1 + 4 omega^2 tau^2) = (A^2 tau^2)/(2 m^2 omega^2 (1 + 4 omega^2 tau^2))\
    P_(0->n(n!=2))^((1)) = 0\
  $
  所有激发态中，只有$n = 2$态有非零跃迁概率。
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
  + At $t = 0$ the system is known to be in the first state, represented by
    $
      mat(1; 0)
    $
    Using time-dependent perturbation theory and assuming that $E_1^0- E_2^0$ is not   close to $±hbar omega$, derive an expression for the probability that the system is found in the second state represented by
    $
      mat(0; 1)
    $
    as a function of $t(t > 0)$.
  + Why is this procedure not valid when $E_1^0- E_2^0$ is close to $±hbar omega$?
]

#solution[
  + 相互作用绘景一阶跃迁振幅为
  $
    cal(A)_"21"^((1)) (t) = - i/hbar integral_0^t dd(t_1) e^(i/hbar (E_2^0 - E_1^0) t_1) braket(2, hat(V)(t_1), 1)\
  $
  计算矩阵元
  $
    braket(2, hat(V)(t_1), 1) = lambda cos(omega t) \
  $
  因此
  $
    cal(A)_"21"^((1)) (t) &= - (i lambda)/hbar integral_0^t dd(t_1) e^(i/hbar (E_2^0 - E_1^0) t_1) cos(omega t_1) \
    &= - (i lambda)/(2 hbar) integral_0^t dd(t_1) (e^(i/(hbar) (E_2^0 - E_1^0 + hbar omega) t_1) + e^(i/(hbar) (E_2^0 - E_1^0 - hbar omega) t_1)) \
    &= (lambda)/(2 hbar) [(e^(i(omega_(21) + omega) t) - 1)/(omega_(21) + omega) + (e^(i(omega_(21) - omega) t) - 1)/(omega_(21) - omega)]\
    &= lambda/(2 hbar) ((2i e^(i (omega_(21) + omega) t/2) sin((omega_(21) + omega) t/2))/(omega_(21) + omega) + (2i e^(i (omega_(21) - omega) t/2) sin((omega_(21) - omega) t/2))/(omega_(21) - omega))\
  $
  其中
  $
    omega_(21) = (E_2^0 - E_1^0)/hbar
  $
  由此得到跃迁概率
  $
    P_(1->2)^((1)) (t) &= abs(cal(A)_"21"^((1)) (t))^2 \
    &= (lambda^2)/(4 hbar^2) abs((2i e^(i (omega_(21) + omega) t/2) sin((omega_(21) + omega) t/2))/(omega_(21) + omega) + (2i e^(i (omega_(21) - omega) t/2) sin((omega_(21) - omega) t/2))/(omega_(21) - omega))^2\
    &= (lambda^2)/(hbar^2) ((sin^2((omega_(21) + omega) t)/2)/(omega_(21) + omega)^2 + (sin^2((omega_(21) - omega) t)/2)/(omega_(21) - omega)^2 + (2 sin ((omega_(21) + omega) t)/2 sin ((omega_(21) - omega) t)/2 cos omega t)/((omega_(21) + omega)(omega_(21) - omega)))\
  $
  + 在$E_1^0 - E_2^0$ 接近 $±hbar omega$ 时，分母 $omega_(21) ± omega$ 接近零，导致跃迁概率发散，说明微扰论失效。
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
