# Numerical Methods for Poroelasticity

Codes to solve the Poroelasticity problem (Biot's Consolidation Theory) numerically.

## Mathematical Model

The one-dimension model for poroelasticity is gave by

```math
-E\frac{\partial^2u}{\partial x^2} + \frac{\partial p}{\partial x} = U
```

```math
\frac{\partial}{\partial t}\left(\frac{\partial u}{\partial x}\right) - K\frac{\partial^2p}{\partial x^2} = P
```

We will assume the following bondaries conditions ond the left side

```math
E\frac{\partial u(0,t)}{\partial x} = 0
```

```math
p(0,t) = 0
```

on the right side

```math
u(L,t) = 0
```

```math
K\frac{\partial p(L,t)}{\partial x} = 0
```

## Numerical Model

Using the Finite Volume Method (FVM), we obtain the following numerical model for one-dimension poroelasticity

```math
u_i = \frac{u_{i-1}+u_{i+1}}{2(\Delta x)^2} + \frac{p_{i-1}-p_{i+1}}{8E} + \frac{(\Delta x)^2}{2E}U(x_i,t_i)
```

```math
p_i = \frac{E \Delta x}{4EK \Delta t + (\Delta x)^2}(u_{i-1}-u_{i+1}) + (u^o_{i-1}-u^o_{i+1}) + \frac{4EK \Delta t + (\Delta x)^2}{E \Delta x}(p_{i-1} + p_{i+1}) + \frac{\Delta x}{2E}(p^o_{i-1}-2p^o_i+p^o_{i+1}) + 4\Delta x \Delta tP(x_i,t_i))
```

