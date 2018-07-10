module FirstAndLastPlaceValues

export ufp, ulp, ulps,
       rreu, 𝐮           # relative rounding error unit: 𝐮 (\bfu), 𝐮(T) == 2.0^(-precision(T))
       subnmin, η        # minimum positive subnormal:   η (\eta), η(T)

import Base: IEEEFloat, prevfloat

prevfloat(x::T, n::Int) where {T<:AbstractFloat} = -nextfloat(-x, n)

𝐮(::Type{Float64}) = ldexp(1.0, -precision(Float64))
𝐮(::Type{Float32}) = ldexp(1.0, -precision(Float32))
𝐮(::Type{Float16}) = ldexp(1.0, -precision(Float16))

η(::Type{Float64}) =  ldexp(0.5, -1073)
η(::Type{Float32}) =  ldexp(0.5, -148)
η(::Type{Float16}) =  ldexp(0.5, -23)

const rreu = 𝐮
const subnmin = η

# nominal `ulp` values for IEEEFloat Types
const NominalSignificand = 0.5 # or 1.0 (per paper)

ulp(::Type{T}) where {T<:IEEEFloat} =
    inv(ldexp(NominalSignificand, precision(T)))


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
