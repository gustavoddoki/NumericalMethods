def compute_results(N, m, E, K, solver):
  """
  This code compute the solution vectors for the one-dimension poroelasticity problem
  :param N: quantity of control volumes
  :param m: quantity of time steps
  :param E: elasticity modulus (Pa)
  :param K: hydraulic conductivity (m/s)
  :param solver: solver used to compute the results for each time step
  :return: solution vector of displacement variable, solution vector of pressure variable, sinaliize if the solution converged
  """

  # Flag to sinalize if the solution converges
  flag_converge = True

  # Lenght of each control volume
  dx = L / N

  # Interval between each time step
  dt = tf / m

  # Discretize the domain
  x = np.arange(- dx / 2, L + dx, dx) # Espacial (Método dos volumes fictícios)
  t = np.linspace(0, tf, m + 1)       # Temporal

  # Initial conditions
  u0 = np.cos(np.pi * x / (2 / L)) # Deslocamento
  p0 = np.sin(np.pi * x / (2 / L)) # Pressão

  # Empty lists to store U and P values for each position and time
  U = []
  P = []

  # Compute the value of U and P for each position and time
  c = 0
  while c <= m:
    U.append((E * np.pi / (2 * L) + 1) * np.pi / (2 * L) * np.cos(np.pi / (2 * L) * x) * np.exp(- t[c]))
    P.append((1 + K * np.pi / (2 * L)) * np.pi / (2 * L) * np.sin(np.pi / (2 * L) * x) * np.exp(- t[c]))
    c += 1

  # First time step
  c = 1

  # Set the estimative of the solutiions of the first time step as the initial conditiions
  u = u0.copy()
  p = p0.copy()


  # Function to compute the pause criteria used in this code (can be other one)
  def compute_pause_criteria(Ru, Rp, Rui, Rpi):
    """
    This function computes the pause criteria
    :param Ru: residual of the solution vector of the displacement variable
    :param Rp: residual of the solution vector of the pressure variable
    :param Rui: residual of the estimative vector of the displacement variable
    :param Rup: residual of the estimative vector of the pressure variable
    :return: pause criteria (ratio between the sum of infinity norms of the solution vectors and the sum of infinity norms of estimative vectors)

    # Infinity norm of the residuals of the solution vectors
    Ru_norm = np.linalg.norm(Ru, np.inf)
    Rp_norm = np.linalg.norm(Rp, np.inf)

    # Infinity norm of the resifuals of the estimative vectors
    Rui_norm = np.linalg.norm(Rui, np.inf)
    Rpi_norm = np.linalg.norm(Rpi, np.inf)

    return (Ru_norm + Rp_norm) / (Rui_norm + Rpi_norm)


  # While the time step is lower or equal the last time step
  while c <= m:

    # Compute the residuals of the estimative vectors
    Rui, Rpi = calcular_residuo(u, p, c, u0, p0, N, E, K, U, P)

    # Infinity loop that breaks when the pause criteria satisfies the tolerance or the solution not converges
    while True:

      # Compute the solutions
      u, p = solver(u, p, c, u0, p0, N, E, K, U, P, dx, dt)

      # Compute the residuals of the solutions vectors
      Ru, Rp = calcular_residuo(u, p, c, u0, p0, N, E, K, U, P)

      # Pause criteria
      pause_criteria = compute_pause_criteria(Ru, Rp, Rui, Rpi)

      # Check if the pause criteria satisfies the tolerance or the pause criteria explodes
      if pause_criteria < tol or pause_criteria > 1e50:
        # Break the loop
        break

    # Check if pause criteria explodes (when the solution doesn't converge)
    if pause_criteria > 1e50:
      # Break the loop and set flag_converge as False
      flag_converge = False
      break

    # Set the initial conditions of the next time step as the solution of this time step
    u0 = u.copy()
    p0 = p.copy()

    # Next time step
    c += 1

  return u, p, flag_converge
