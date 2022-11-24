@testset "Test Dense Lanczos" begin

    using KrylovKit
    H = OperatorMap(Dict("XX"=>1,"Z"=>-1.05,"X"=>0.5))
    L = Liouvillian(H)
    O0 = OperatorMap(Dict("XX"=>1.05,"Z"=>1))

    #set up Lanczos
    iter = LanczosIterator(L,O0);

    #run Lanczos
    N=22
    n=0
    F = initialize(iter)
    while n < N
        expand!(iter, F)
        n+=1
    end
    bs = F.βs

    #compare to correct values
    Xiangyu_values = vec([4.15880082  4.31401627  4.50365026  4.79068279  5.15486467  5.58710085 6.06482781  6.68292298  7.41393735  7.90678234  8.0922167   8.40878809 8.84958688  9.2086087   9.59364746 10.07842675 10.53190668 10.97745372 11.33212599 11.6183719])
        # check = foldl( (bool,index) -> bool && isapprox(betas[index],Xiangyu_values[index]), 1:min(layers,length(Xiangyu_values)); init=true)
        # println("same values? ", check)

    for n in 1:20
        @test abs(bs[n]-Xiangyu_values[n]) < 1e-8
    end
end