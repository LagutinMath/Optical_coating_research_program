from deposition_simulation import simulation, Design, SetUpParameters

def simulation_witness_monitoring(des_th: Design, term_algs: list, set_up_pars: SetUpParameters) -> list:
    witness_des = [des_th.design_on_witness(j) for j in range(1, des_th.N, 2)]
    sim_list = [simulation(des, term_algs, set_up_pars) for des in witness_des]
    error_d = []
    for each in sim_list:
        error_d += each.errors_d[1:]
    return error_d
