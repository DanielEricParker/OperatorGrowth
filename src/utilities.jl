random_pauli_string(n) = reduce(*,rand(["I","X","Y","Z"],n))
random_pauli(n) = PS(random_pauli_string(n))

function r_statistic(eigvals)
	#must be sorted
	eigenvalues = sort(eigvals)
	rs = Array{Float64}(undef, length(eigenvalues)-2)

	for i in 1:length(rs)
		delta_na = eigenvalues[i+1]-eigenvalues[i]
		delta_nb = eigenvalues[i+2]-eigenvalues[i+1]
		r_min = min(delta_nb,delta_na)
		r_max = max(delta_nb,delta_na)
		r = (abs(r_min - r_max) < 10e-14) ? 0.0 : r_min/r_max
		#just in case, check we're in the bounds
		if r < -0.001 || r > 1.001
			error("out of bounds", i," ", r_min," ", r_max," ", r)
		end
		rs[i] = r
	end
	return mean(rs)
end