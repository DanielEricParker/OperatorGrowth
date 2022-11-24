#######################################################################
## Operator ###########################################################
#######################################################################

struct OperatorMap
    terms :: Dict{Vertex,ComplexF64}
    OperatorMap(terms :: Dict{Vertex,ComplexF64}) = new(terms)
    OperatorMap() = new(Dict{Vertex,ComplexF64}())
end

function Base.:+(O :: OperatorMap, P :: PS)
    Od = O.terms
    Od[P.V] = haskey(Od,P.V) ? Od[P.V] + P.z : P.z
    return O
end

function OperatorMap(M :: Dict{String,T}) where T <: Number
    O = OperatorMap(Dict{Vertex,ComplexF64}())
    for (string,z) in M
        O += z*PS(string)
    end
    return O
end

Base.length(O :: OperatorMap) = length(O.terms)
Base.iszero(O :: OperatorMap) = (length(O) == 0)

#######################################################################
## Commutation ---- Dictionary Version ################################
#######################################################################

function commutator!(O :: OperatorMap, V1 :: Vertex, V2 :: Vertex, offset, λ)
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
        z = 2*ifelse(ϵ,-λ,λ) #λ * im or λ
        O += PS(z,newVertex)
    end
    return O
end

function convolve!(O :: OperatorMap, P1 :: PS, P2 :: PS)
    e1 = edges(P1.V)
    e2 = edges(P2.V)
    min_offset = e2[1] - e1[2]
    max_offset = e2[2] - e1[1]
    λ = P1.z * P2.z
    for offset in min_offset:max_offset
        commutator!(O, P1.V, P2.V, offset, λ)
    end
    return O
end
convolve(P1 :: PS, P2 :: PS) = convolve!(OperatorMap(Dict()),P1,P2)

#######################################################################
## Operator Algebra ###################################################
#######################################################################

#The following methods are required for KrylovKit
#Base.:*(α, v)
#Base.similar(v)
#LinearAlgebra.mul!(w, v, α)
#LinearAlgebra.rmul!(v, α)
#LinearAlgebra.axpy!(α, v, w)
#LinearAlgebra.axpby!(α, v, β, w)
#LinearAlgebra.dot(v,w)
#LinearAlgebra.norm(v)

function Base.:*(α,v :: OperatorMap)
    w = deepcopy(v)
    map!(x->α*x,values(w.terms))
    return w
end

Base.similar(:: OperatorMap) = OperatorMap()
#preallocating w buys us nothing here =(
function LinearAlgebra.mul!(w :: OperatorMap, v :: OperatorMap, α)
    dw = w.terms
    empty!(dw)
    for (vert,z) in v.terms
        dw[vert] = α*z
    end
    return w
end

function LinearAlgebra.rmul!(v :: OperatorMap, α :: Bool)
    if !α
        empty!(v.terms)
    end
    return v
end
LinearAlgebra.rmul!(v :: OperatorMap, α) = (map!(x->x*α,values(v.terms)); v)

#w <- α*v + w
function LinearAlgebra.axpy!(α, v :: OperatorMap, w ::OperatorMap)
    dv, dw = v.terms, w.terms
    for (vert, zv) in dv
        dw[vert] = haskey(dw,vert) ? dw[vert] + α * zv : α*zv
    end
    return w 
end

function LinearAlgebra.axpby!(α, v :: OperatorMap, β, w :: OperatorMap) #w <- α*v + β*w
    dv,dw = v.terms, w.terms
    for (vert, zw) in dw
        dw[vert] = haskey(dv,vert) ? α*dv[vert] + β*zw : β*zw
    end
    for (vert, zv) in dv
        if !haskey(dw,vert)
            dw[vert] = α*zv
        end
    end
    return w
end

function LinearAlgebra.dot(v :: OperatorMap,w::OperatorMap)
    dv,dw = v.terms, w.terms
    (length(dv) == 0 || length(dw) == 0) && return zero(Float64)
    return sum(dot(zv,get(dw,vert,zero(ComplexF64))) for (vert,zv) in dv)
end

LinearAlgebra.norm(v :: OperatorMap) = sqrt(sum(abs2.(values(v.terms))))

function trim!(v :: OperatorMap, ϵ)
    for (vert,z) in v.terms
        if abs(z) < ϵ
            delete!(v.terms, vert)
        end
    end
    return v
end