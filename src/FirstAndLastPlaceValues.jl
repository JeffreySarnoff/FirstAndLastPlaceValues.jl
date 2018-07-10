module FirstAndLastPlaceValues

export ufp, ulp

import Base: IEEEFloat, prevfloat

if !hasmethod(prevfloat, (Float64, Int))
    prevfloat(x::T, n::Int) where {T<:AbstractFloat} = -nextfloat(-x, n)
end

# ufp is "unit first place"

ufp(x::T) where {T<:IEEEFloat} = x !== zero(T) ? ldexp(one(T), exponent(x)) : x
ufp(x::Float64) = x !== 0.0 ? ldexp(1.0, exponent(x)) : x
ufp(x::Float32) = x !== 0.0 ? ldexp(1.0f0, exponent(x)) : x
ufp(x::Float16) = x !== 0.0 ? ldexp(one(Float16), exponent(x)) : x


# ulp is "unit last place"
# ulp(T) == ulp(one(T)) == ldexp(one(T), 1-precision(T))

const ulp_Float64 = ldexp(1.0, -52)
const ulp_Float32 = ldexp(1.0f0, -23)
const ulp_Float16 = ldexp(one(Float16), -10)

# acutal `ulp` (unit last place) values for IEEEFloats

ulp(x::Float64) = x !== 0.0 ? ulp_Float64 * ldexp(1.0, exponent(x)) : 0.0
ulp(x::Float32) = x !== 0.0f0 ? ulp_Float32 * ldexp(1.0f0, exponent(x)) : 0.0f0
ulp(x::Float16) = x !== Float16(0.0) ? ulp_Float16 * ldexp(one(Float16), exponent(x)) : zero(Float16)


#=

# 𝐮 is the relative rounding error unit
# 𝐮 is half the distance between 1.0 and its successor
# 𝐮 == (nextfloat(1.0) - 1.0) / 2 == eps(1.0) / 2
# 𝐮(T) == ldexp(one(T), -precision(T))

const 𝐮_Float64 = ldexp(1.0, -precision(Float64))
const 𝐮_Float32 = ldexp(1.0f0, -precision(Float32))
const 𝐮_Float16 = ldexp(one(Float16), -precision(Float16))

# 𝛈 is the smallest subnormal
# 𝛈 == nextfloat(0.0)

const 𝛈_Float64 =  ldexp(1.0, -1074)
const 𝛈_Float32 =  ldexp(1.0f0, -149)
const 𝛈_Float16 =  ldexp(one(Float16), -24)



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

=#
end # FirstAndLastPlaceValues
