#import "@preview/scripst:1.1.1": *
#import "@preview/physica:0.9.7": *

#show: scripst.with(
  title: [高等量子力学第10次作业],
  author: "AnZrew",
  time: "2025年12月",
)

#exercise[
  For $H = vb(P)^2/2m + V(vb(X))$, write down the standard path integral with $integral dd(vb(p))$ and $integral dd(vb(x))$. How to obtain an exactly equivalent path integral with $integral dd(vb(x))$ only?
]

#exercise[
  For a more general Hamiltonian, this cannot be done exactly. Usually the saddle point approximation is used, i.e. substitute for $vb(P) = (...)$ using the classical equation of motion (this is only exact in the case considered above). Try this for a relativistic particle $H = sqrt((vb(P)c)^2 + (m c^2)^2)$.
]
