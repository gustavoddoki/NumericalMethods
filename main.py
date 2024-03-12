from SolverGaussSeidel import gauss_seidel_solver 


def main():

  N = 128
  m = 256

  E = 1
  K = 1

  u, p, converge = compute_results(N, m, E, K, gauss_seidel_solver, tol=1e-5)

  dict_dados = {'u': u, 'p': p}
  pd.DataFrame.from_dict(dict_dados).to_excel(f'GS.xlsx')

  return None


if __name__ == "__main__":
  main()
