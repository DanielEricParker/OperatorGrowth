struct Basis{T}
    name :: String
    ind :: Dict{T,Int64}
    vec :: Dict{Int64,T}
    partitions :: Array{Int64,1}
end

Base.length(b :: Basis) = length(b.ind)

function partitions(b :: Basis)
    pa = b.partitions
    ranges = [n==1 ? (1:pa[1]) : (pa[n-1]+1:pa[n]) for n in eachindex(pa)]
    return ranges
end

function list_paulis(nmax; has_id=true) 
    ops_lists = [BPS[] for _ in 1:nmax]
    has_id && push!(ops_lists[1],BPS(false,false,Vertex(0,0)))
    ra = 0:((1<<nmax)-1)
    # ys_coefficient = ComplexF64[im,-1,-im,1]
    ys_coefficient = [(true,false),(false,true),(true,true),(false,false)]
    for w in ra, v in ra
        if (v & 1) | (w & 1) !=0
            onesMask = ~(v|w)
            stop = (63 - leading_ones(onesMask))
            δϵ = ys_coefficient[mod1(count_ones(v&w),4)]
            P = BPS(δϵ...,Vertex(v,w))
            push!(ops_lists[stop+1],P)
        end
    end
    partitions = cumsum(length.(ops_lists))
    return reduce(vcat,ops_lists), partitions
end

function pauli_basis(nmax;has_id=false)
    name = "Pauli Basis $(nmax)"
    lp,partitions = list_paulis(nmax; has_id)
    ra = 1:length(lp)
    indices = Dict{BPS,Int64}(lp .=> ra)
    vectors = Dict{Int64,BPS}(ra .=> lp)
    return Basis(name,indices,vectors, partitions)
end


function build_graph(L :: Liouvillian, basis)
    Is = Int64[]
    Js = Int64[]
    Ks = ComplexF64[]
    for (col,BPScol) in basis.vec
        o = PS(BPScol)
        O = OperatorMap(Dict(o.V => o.z))
        LO = L * O
        for (vert,z) in LO.terms
            z, BPSrow = split(PS(z,vert))
            row = get(basis.ind,BPSrow,0)
            if col > row #no diagonal for L
                if row != 0
                    push!(Is,row)
                    push!(Js,col)
                    push!(Ks,z)

                    push!(Is,col)
                    push!(Js,row)
                    push!(Ks,conj(z))
                end
            end
        end
    end
    return sparse(Is,Js,Ks,length(basis),length(basis))
end

function Base.vec(O :: OperatorMap, b :: Basis)
    v = zeros(ComplexF64,length(b))
    for (vect,z) in O.terms
        z,bps = split(PS(z,vect))
        v[b.ind[bps]] = z
    end
    return v
end