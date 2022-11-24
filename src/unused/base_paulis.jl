#unused 


#######################################################################
#################### Base Pauli Strings ###############################
#######################################################################

@doc raw"""
`BPS(δ, epsilon, v, w, len)``

Type for storing strings of Pauli matrices. The representation used is
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
struct PauliOp
    δ :: Bool
    ϵ :: Bool
    v :: UInt64
    w :: UInt64

    PauliOp(prefactor,δ,ϵ,start,stop,v,w) = new(prefactor,δ,ϵ,start,stop,v,w)
	#user-friendly constructor
	function PauliOp(name :: String)
		#icky conversion from SubStrings to Strings is necessary
		return NamesToPauliOp(map(String,split(name,"")))
	end
    function PauliOp(v,w)
        stop = count_ones(v|w)-1
        new(1.0,false,false,0,stop,v,w)
    end
end


# #######################################################################
# ## Internal Pauli String ##############################################
# #######################################################################

# struct PSInternal
#     δ :: Bool
#     ϵ :: Bool
#     V :: Vertex
# end

# Base.zero(::Type{PSInternal}) = PSInternal(false,false,zero(Vertex))
# Base.iszero(P :: PSInternal) = iszero(P.V)
# Base.show(io :: IO, P :: PSInternal) = show(io,PS(P))

# function PS(P :: PSInternal)
#     z = ifelse(P.δ,complex(0,1.0),complex(1.0,0))
#     z = ifelse(P.ϵ,-z,z)
#     return PS(z,P.V)
# end