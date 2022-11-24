using OperatorGrowth
using Test
using LinearAlgebra
using KrylovKit
using SparseArrays


@testset "Running all tests..." begin
	# include("test_pauli_mult.jl")
	include("test_Lanczos.jl")
end