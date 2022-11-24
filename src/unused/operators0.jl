#unused 

#Array Version
struct Operator
    terms :: Array{PS,1}
end

#######################################################################
## Commutation ---- Array Version #####################################
#######################################################################

function commutator!(buffer, V1 :: Vertex, V2 :: Vertex)
    d1 = count_ones(V1.w & V2.v) % 2
    d2 = count_ones(V2.w & V1.v) % 2
    if d1 != d2
        ϵ = d1 == 1
        newVertex = Vertex(xor(V1.v,V2.v), xor(V1.w,V2.w))
        push!(buffer, PSInternal(false, ϵ, newVertex))
    end
    return buffer
end
commutator(V1 :: Vertex, V2 :: Vertex) = commutator!(PSInternal[],V1,V2)

function commutator!(buffer :: Array{PS,1}, V1 :: Vertex, V2 :: Vertex, offset, λ)
    if offset < 0
        V2 = Vertex(V2.v >>> offset, V2.w >>> offset)
    elseif offset > 0
        V1 = Vertex(V1.v << offset, V1.w << offset)
    end
    d1 = count_ones(V1.w & V2.v) % 2
	d2 = count_ones(V2.w & V1.v) % 2
    if d1 != d2
        ϵ = (d1==1)
        newVertex = rightAlign(Vertex(xor(V1.v,V2.v), xor(V1.w,V2.w)))
        z = ifelse(ϵ,-λ,λ)
        P = PS(z,newVertex)
        push!(buffer,P)
    end
end

function convolve!(O :: Operator, P1 :: PS, P2 :: PS)
    buffer = O.terms
    e1 = edges(P1.V)
    e2 = edges(P2.V)
    min_offset = e2[1] - e1[2]
    max_offset = e2[2] - e1[1]
    λ = P1.z * P2.z
    for offset in min_offset:max_offset
        commutator!(buffer, P1.V, P2.V, offset, λ)
    end
    return O
end
convolve(P1 :: PS, P2 :: PS) = convolve!(Operator(PS[]),P1,P2)