# Numerical Methods for Poroelasticity

Codes to solve the Poroelasticity problem (Biot's Consolidation Theory) numerically.

## Mathematical Model

The one-dimension model for poroelasticity is gave by the displacement equation

```math
-E\frac{\partial^2u}{\partial x^2} + \frac{\partial p}{\partial x} = U
```

and the pressure equation

```math
\frac{\partial}{\partial t}\left(\frac{\partial u}{\partial x}\right) - K\frac{\partial^2p}{\partial x^2} = P
```

where $E$ denotes the elastic modulus, $K$ represents hydraulic conductivity, $U$ signifies the density of the force applied to the body, and $P$ denotes the injection or extraction force of the fluid within the porous medium. The displacement is defined by \(u(x,t)\), and the pressure by $p(x,t)$, where $x$ is the spatial variable, and $t$ is the temporal variable. We consider a spatial domain $[0,L]$ and a temporal domain $[0,tf]$.

For the boundary conditions, we will assume free drainage without variation in displacement on the left boundary

```math
E\frac{\partial u(0,t)}{\partial x} = 0
```

```math
p(0,t) = 0
```

and stiffness without pressure variation on the right boundary

```math
u(L,t) = 0
```

```math
K\frac{\partial p(L,t)}{\partial x} = 0
```

## Numerical Model

To enhance system stability during the numerical simulation without compromising results a smoothing term is coupled on the left side of the pressure equation

```math
-\frac{\left(\Delta x\right)^2}{4E}\frac{\partial}{\partial t}\left(\frac{\partial^2p}{\partial x^2}\right)
```

Using the Finite Volume Method (FVM) and Implicit Euler Method, we obtain the following numerical model for one-dimension poroelasticity

```math
u_i = \frac{u_{i-1}+u_{i+1}}{2(\Delta x)^2} + \frac{p_{i-1}-p_{i+1}}{8E} + \frac{(\Delta x)^2}{2E}U(x_i,t_i)
```

```math
p_i = \frac{E \Delta x}{4EK \Delta t + (\Delta x)^2}(u_{i-1}-u_{i+1}) + (u^o_{i-1}-u^o_{i+1}) + \frac{4EK \Delta t + (\Delta x)^2}{E \Delta x}(p_{i-1} + p_{i+1}) + \frac{\Delta x}{2E}(p^o_{i-1}-2p^o_i+p^o_{i+1}) + 4\Delta x \Delta tP(x_i,t_i))
```

To implement the boundary conditions, we will a fictitious control volume on the left boundary

```math
u_0 = u_1 
```

```math
p_0 = -p_1
```

and on the right condition

```math
u_{N+1} = -u_{N}
```

```math
p_{N+1} = p_{N}
```


