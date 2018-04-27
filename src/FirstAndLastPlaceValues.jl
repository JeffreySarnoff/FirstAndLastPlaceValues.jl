module FirstAndLastPlaceValues

export ufp, ulp, ulps

import Base: IEEEFloat, prevfloat

prevfloat(x::T, n::Int) where {T<:AbstractFloat} = -nextfloat(-x, n)

# nominal `ulp` values for IEEEFloat Types
ulp(::Type{T}) where {T<:IEEEFloat} =
    inv(ldexp(1.0, precision(T)))


# actual `ufp` (unit first place) values for IEEEFloats 
function ufp(x::T) where {T<:IEEEFloat}
    p = precision(T) - 1
    u = reinterpret(Unsigned, x)
    u = (u >> p) << p
    return reinterpret(T, u)
end

# acutal `ulp` (unit last place) values for IEEEFloats
ulp(x::T) where {T<:IEEEFloat} = ufp(x) * ulp(T)

#=
     let  n = binade_ulps( refval, obsval ) 
          where obsval is in the same binad as refval
     then nextfloat(min(refval, obsval), n) == max(refval, obsval)
=#     
function binade_ulps(reference_value::T, investigative_value::T) where {T<:AbstractFloat}
    distance = reference_value - investigative_value
    separation = signbit(distance) ? abs(distance)/2 : distance
    reference_ulp = ulp(reference_value)
    relative_ulps  = round(Int, separation / reference_ulp)
    return relative_ulps 
end



# provision for BigFloat
# CAUTION: BigFloat can misrepresent a few trailing digits

ulp(x::BigFloat) = eps(x)/2
ufp(x::BigFloat) = ulp(x) * frexp(1.0, precision(x))


end # FirstAndLastPlaceValues
