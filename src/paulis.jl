#######################################################################
## Vertex #############################################################
#######################################################################

struct Vertex
    v :: UInt64
    w :: UInt64
end

Base.zero(::Type{Vertex}) = Vertex(0,0)
Base.iszero(V :: Vertex) = (V.v == 0) && (V.w == 0)
#custom hash for speed. 10x speedup
Base.hash(V :: Vertex) = hash(V.v ⊻ hash(V.w))

function edges(V :: Vertex)
    onesMask = V.v|V.w
    return ifelse(onesMask==0,(0,0),(trailing_zeros(onesMask),63-leading_zeros(onesMask)))
end
function rightAlign(V :: Vertex)
    start = trailing_zeros(V.v|V.w)
    return ifelse(start==0,V,Vertex(V.v>>>start,V.w>>>start))
end

function vertexStringHelper(V :: Vertex)
    (start, stop) = edges(V)
    numYs = count_ones(V.v & V.w)
    vs = digits(V.v,base=2,pad=stop+1)
	ws = digits(V.w,base=2,pad=stop+1)

	names = ["I" "Z"; "X" "Y"]
    arr = [names[vs[i]+1,ws[i]+1] for i in start+1:stop+1]
	opnames = foldr(*,arr)
	leading = start > 0 ? (start > 1 ? "I^$(start)" : "I") : ""
    return numYs,leading*opnames
end

function Base.show(io :: IO, V :: Vertex) 
    if iszero(V)
        print(io,"I")
    else
        numYs, vertexstring = vertexStringHelper(V)
        print(io,"$(-im^numYs)*"*vertexstring)
    end
end

#######################################################################
## Pauli String #######################################################
#######################################################################

struct PS
    z :: ComplexF64
    V :: Vertex
end


function PS(name :: String)
    names = map(String,split(name,""))
    numYs = 0
    v = 0
    w = 0
    lettersToXZs = Dict("I" => (0,0), "X"=>(1,0), "Z"=>(0,1),"Y"=>(1,1))
    for (i,name) in enumerate(names)
        !haskey(lettersToXZs,name) && error("Invalid Pauli letter $(name).")
        (vi,wi) = lettersToXZs[name]
        v = v | (vi << (i-1))
        w = w | (wi << (i-1))
    end
    numYs = count_ones(v&w)%4
    return PS((im)^numYs,Vertex(v,w))
end

PS(P :: Pair{Vertex,ComplexF64}) = PS(P[2],P[1])

Base.zero(::Type{PS}) = PS(0,zero(Vertex))
Base.iszero(P :: PS) = iszero(P.V) || P.z == zero(ComplexF64)

Base.:*(λ, P :: PS) = PS(λ*P.z,P.V)
Base.isapprox(P1 :: PS, P2 :: PS; kwargs...) = P1.V == P2.V && isapprox(P1.z,P2.z; kwargs...)

function canonicalForm(P :: PS)
    (start, stop) = edges(P)
    if start > 0
        return PS(z,Vertex(P.V.v >> start, P.V.w >> start))
    else
        return P
    end
end

function Base.show(io :: IO, P :: PS)
    if iszero(P)
        print(io,"I")
    else
        numYs, vertexstring = vertexStringHelper(P.V)
        print(io,"$(P.z * (-im)^numYs)*"*vertexstring)
    end
end

function showVerbose(P :: PS)
	(start,stop) = edges(P.V)
	print("PS[$(P.z), $(digits(P.V.v,base=2,pad=stop+1)), $(digits(P.V.w,base=2,pad=stop+1)), $(start), $(stop)]\n")
end

#######################################################################
#################### Base Pauli Strings ###############################
#######################################################################


#######################################################################
#################### Base Pauli Strings ###############################
#######################################################################

@doc raw"""
`BPS(δ, epsilon, v, w, len)``

Type for storing Hermitian basis vectors of Pauli matrices. The representation used is
``
i^{\δ} * (-1)^{\epsilon} * \prod_{i=1}^{|v|} Z_i^{v_i} X_i^{w_i}.
``
Also, `start` is the first non-trivial site and `stop` is the last non-trivial site.
These are indexed in ``[0,63]``. For the operator "I^64", one has the convention `start = 0, stop = 0`.

Ref.
    Jeroen Dehaene and Bart De Moor
    Clifford group, stabilizer states, and linear and quadratic operations
    over GF(2)
    Phys. Rev. A 68, 042318 – Published 20 October 2003
"""
struct BPS
    δ :: Bool
    ϵ :: Bool
    V :: Vertex
end

Base.zero(::Type{BPS}) = BPS(false,false,zero(Vertex))
Base.iszero(P :: BPS) = iszero(P.V)
Base.hash(BP :: BPS) = hash((BP.δ << 1 ⊻ BP.ϵ) ⊻ hash(BP.V))
#ordering on vertices to allow sorting

function showVerbose(P :: BPS)
	(start,stop) = edges(P.V)
	print("BPS[$(P.δ), $(P.ϵ), $(digits(P.V.v,base=2,pad=stop+1)), $(digits(P.V.w,base=2,pad=stop+1))]\n")
end

function PS(BP :: BPS)
    # z = ifelse(BP.δ,0.0+1.0im,1.0+0.0im)
    # z = ifelse(BP.ϵ,-z,z)
    z = (im)^BP.δ * (-1)^BP.ϵ
    return PS(z,BP.V)
end

# """Given the number of Ys `numYs`, return ``(i)^δ*(-1)^ϵ = (-i)^numY``."""
# function numYsToδϵ(numYs)
#     δ = Bool(numYs & 1)
#     ϵ = Bool(xor(numYs & 1, (numYs>>1) & 1))
#     return δ,ϵ
# end

# P(z,V) = z' * BPS(δ,ϵ,V) s.t. BPS(δ,ϵ,V) is Hermitian
# M = numY(V)
# P = z * [(-iY)^M (XZ)^...] = (z * (-i)^M) * [(Y)^M (XZ)^...]
function Base.split(P :: PS)
    z = P.z
    numYs = mod(count_ones(P.V.v & P.V.w), 4)
    if numYs == 1
        z = complex(z.im,-z.re) #z * (-i)^1 = y - i x
        δ,ϵ = true, false #i^1
    elseif numYs == 2 #z*(-i)^2 = z*(-1) = -x - i y
        z = complex(-z.re,-z.im)
        δ,ϵ = false, true #i^1 = -1
    elseif numYs == 3 #z*(-i)^3 = z*i = -y + i x
        z = complex(-z.im,z.re)
        δ,ϵ = true, true #i^3 = -i
    else #z*(-i)^4 = z
        z = z
        δ,ϵ = false, false
    end
    return z, BPS(δ,ϵ,P.V)
end

Base.show(io :: IO, P :: BPS) = show(io, PS(P))