def gauss_seidel_solver(u, p, j, u0, p0, N, U, P, c1, c2, c3, c4, teta=1):
  """
  This function update the solution vector of the equation
  :param u: solution vector for displacement variable
  :param p: solution vector for pressure variable
  :param j: time step
  :param u0: solution vector of previous time step for displacement variable (initial condition for the first time step)
  :param p0: solution vector of previous time step for pressure variable (initial condition for the first time step)
  :param N: números de volumes de controle
  :param U: vector of density of the force applied to the body
  :param P: vector of injection or extraction force of the fluid
  :param c1: coefficient 1
  :param c2: coefficient 2
  :param c3: coefficient 3
  :param c4: coefficient 4
  :param teta: euler method param (implicit -> teta=1)
  :return: updated vectors for displacement and pressure variables
  """

  # Loop to update displacement variables
  for i in range(N + 2):
    # Left control volume
    if i == 0:
      # Left boundary condition
      u[i] = u[i + 1]
    # Right control volume
    elif i == N + 1:
      # Right boundary condition
      u[i] = - u[i - 1]
    # Middle control volumes
    else:
      u[i] = 1 / (2 * c1) * (c1 * (u[i - 1] + u[i + 1]) + c2 * (p[i - 1] - p[i + 1]) + U[j][i])

  # Loop to update pressure variables
  for i in range(N + 2):
    # Left control volume
    if i == 0:
      # Left boundary condition
      p[i] = - p[i + 1]
    # Right control volume
    elif i == N + 1:
      # Right boundary condition
      p[i] = p[i - 1]
    # Middle control volumes
    else:
      p[i] = 1 / (2 * c3) * (c2 * (u[i - 1] - u[i + 1]) + c3 * (p[i - 1] + p[i + 1]) + c2 * (
              u0[i + 1] - u0[i - 1]) + c4 * (p0[i - 1] - 2 * p0[i] + p0[i + 1]) + (
              teta * P[j][i] + (1 - teta) * P[j - 1][i]) * dt)

  return u, p
