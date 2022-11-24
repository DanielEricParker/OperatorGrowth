#######################################################################
## Liouvillian ########################################################
#######################################################################

struct Liouvillian
    H :: OperatorMap
end


function LinearAlgebra.mul!(LO :: OperatorMap, L :: Liouvillian, O ::OperatorMap)
    for h in L.H.terms, o in O.terms
        convolve!(LO,PS(h),PS(o))
    end
    return LO
end
Base.:*(L :: Liouvillian, O :: OperatorMap) = mul!(OperatorMap(),L,O)
(L :: Liouvillian)(O :: OperatorMap) = L*O