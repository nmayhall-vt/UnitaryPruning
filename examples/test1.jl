include("../src/read_in.jl")
@time up.compute_expectation_value_recurse(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=1e-7)
@time up.compute_expectation_value_iter(ref_state, ham_ops, ham_par, ansatz_ops, ansatz_par, thresh=1e-7)
