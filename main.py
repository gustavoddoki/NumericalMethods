from ComputeResiduals 


def main():

  N = 128
  m = 256

  

  E = 1
  K = 1

  u, p, converge = compute_results(N, m, E, K, solver, tol=1e-5)

  epoch = 0
  
  A, B, m, n, w, z = 1.585791797112397, 2.3089208021724614, 1.0218350908968998, -1.5703177005933322, 0.8910888059894608, -1.4271235393157982
  print(A, B, m, n, w, z)

  x_domain = np.linspace(0, L, num=size_x, dtype=np.float64)
  t_domain = np.linspace(0, t_f, num=size_t, dtype=np.float64).reshape((size_t, 1))
  x_domain = np.tile(x_domain, (size_t, 1))
  
  batch_vars, random_vars, real_solution = create_sample(L, t_f, E, K, size_x, size_t, batch_size, x_domain, t_domain, [A, B, m, n, w, z])
  
  u_real = real_solution[0]
  p_real = real_solution[1]


  trunk_input = tf.concat([x, t], axis=-1)
  
  inicio = time.time()
  u_model, p_model = model({'branch_input_u0': u0, 'branch_input_p0': p0, 'branch_input_U': U, 'branch_input_P': P, 'trunk_input': trunk_input})
  fim = time.time()
  
  erro_u = np.mean(np.abs(u_model - u_real) / u_real)
  erro_p = np.mean(np.abs(p_model - p_real) / p_real)
  
  print('DeepONet', erro_u, erro_p, fim-inicio)

  dict_dados = {'u': [], 'p': [], 'u_real': [], 'p_real': []}
  for i, j, ii, jj in zip(u_model, p_model, u_real, p_real):
      dict_dados['u'].append(float(i))
      dict_dados['p'].append(float(j))
      dict_dados['u_real'].append(float(ii))
      dict_dados['p_real'].append(float(jj))

  pd.DataFrame.from_dict(dict_dados).to_excel(f'teste/model_resultados_finetuned_E1K1.xlsx')

  for NN, mm in zip([32, 64, 128], [64, 128, 256]):

      # Comprimento de cada volume de controle
      dx = L / NN

      # Intervalo a cada passo de tempo
      dt = t_f / mm

      # Domínio discretizado
      x = np.arange(- dx / 2, L + dx, dx) # Espacial (Método dos volumes fictícios)
      t = np.linspace(0, t_f, mm + 1)      # Temporal

      U_batch2, P_batch2 = [], []
      for time_num in t:
          U_batch, P_batch = calculate_source_terms(L, E, K, x, time_num, A, B, m, n, w, z)
          U_batch2.append(U_batch)
          P_batch2.append(P_batch)
      u0_batch, p0_batch = calculate_initial_condition(L, x, A, B, n, z)
      u_real, p_real = calculate_real_solution(L, x, t[-1], A, B, m, n, w, z)

      dict_dados = {}

      inicio = time.time()
      u, p, converge = calcular_resultado(NN, mm, E, K, atualizar_variaveis_GS, t_f, L, tol, U_batch2, P_batch2, u0_batch, p0_batch)
      fim = time.time()

      erro_u = np.mean(np.abs(u[1:-1] - u_real[1:-1]) / u_real[1:-1])
      erro_p = np.mean(np.abs(p[1:-1] - p_real[1:-1]) / p_real[1:-1])

      dict_dados = {'u': u, 'p': p, 'u_real': u_real, 'p_real': p_real}
      pd.DataFrame.from_dict(dict_dados).to_excel(f'teste/GS_N{NN}_m{mm}_E1K1.xlsx')

      print('GS', NN, mm, erro_u, erro_p, fim-inicio)
