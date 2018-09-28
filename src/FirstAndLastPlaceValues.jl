module FirstAndLastPlaceValues

export ufp, ulp, rre, ε⁻, ε⁺

#=

ACCURATE FLOATING-POINT SUMMATION PART I: FAITHFUL ROUNDING
SIEGFRIED M. RUMP, TAKESHI OGITA, AND SHIN’ICHI OISHI
    http://oishi.info.waseda.ac.jp/~oishi/papers/RuOgOi07I.pdf


Fast quadruple-double floating point format
    Naoya Yamanaka and Shin’ichi Oishi
    2014-Jan-01

ERROR ESTIMATION OF FLOATING-POINT SUMMATION AND DOT PRODUCT
SIEGFRIED M. RUMP

On the definition of ulp (x)
Jean-Michel Muller
=#

two(::Type{T}) where {T<:Real} = one(T) + one(T)
two(::Type{Float64}) = 2.0
two(::Type{Float32}) = 2.0f0
two(::Type{Float16}) = Float16(2.0)

#=
  epsilon (ε)
  ε⁻(T) = (one(T) - prevfloat(one(T)))
  ε⁺(T) = (nextfloat(one(T)) - one(T))
=#
ε⁻(::Type{T}) where {T<:AbstractFloat} = (one(T) - prevfloat(one(T)))
ε⁺(::Type{T}) where {T<:AbstractFloat} = (nextfloat(one(T)) - one(T))

ε⁻(::Type{Float64}) = 1.1102230246251565e-16
ε⁻(::Type{Float32}) = 5.9604645f-8
ε⁻(::Type{Float16}) = Float16(0.0004883)

ε⁺(::Type{Float64}) = 2.220446049250313e-16
ε⁺(::Type{Float32}) = 1.1920929f-7
ε⁺(::Type{Float16}) = Float16(0.000977)


#=
  relative rounding error (rre)
  rre(T) = (nextfloat(one(T)) - one(T)) / (one(T) + one(T))
  rre(T) = eps(T) / 2
=#

rre(::Type{Float64}) = 1.1102230246251565e-16
rre(::Type{Float32}) = 1.1920929f-7
rre(::Type{Float16}) = Float16(0.000977)

rre(::Type{T}) where {T<:AbstractFloat} = 
    (nextfloat(one(T)) - one(T)) / two(T)

#=
inv2rre(T)
=#
inv2rre(::Type{T}) where {T<:AbstractFloat} =
    inv(two(T))*rre(T) + one(T)


#=
the unsigned constants used in `ufp(x::T)`
     Float64: ((~UInt64(0)) >> 52) << 52
     Float32: ((~UInt32(0)) >> 23) << 23
     Float16: ((~UInt16(0)) >> 10) << 10

the shifts: Base.trailing_ones(Base.significand_mask)
=#
     
ufp(x::Float64) =
    reinterpret(Float64, 
        reinterpret(UInt64, x) & 0xfff0000000000000)  
        
ufp(x::Float32) =
    reinterpret(Float32, 
        reinterpret(UInt32, x) & 0xff800000)

ufp(x::Float16) =
    reinterpret(Float16, 
        reinterpret(UInt16, x) & 0xfc00)

function ufp(x::T) where {T<:AbstractFloat}
    q = inv2rre(T) * x
    return q - (one(T) - rre(T))
end



ulp(::Type{T}) where {T<:AbstractFloat} = two(T) * ε⁻(T)

ulp(::Type{Float64}) = 2.220446049250313e-16
ulp(::Type{Float32}) = 1.1920929f-7
ulp(::Type{Float16}) = Float16(0.000977)

ulp(x::T) where {T<:AbstractFloat} = two(T) * ε⁻(T) * ufp(x)
# ulp(x::T) where {T<:AbstractFloat} = ε⁺(T) * ufp(x)

ulp(x::Float64) = 2.220446049250313e-16 * ufp(x)
ulp(x::Float32) = 1.1920929f-7 * ufp(x)
ulp(x::Float16) = Float16(0.000977) * ufp(x)







# ufp is "unit first place"
# ufp(x::T) where {T<:IEEEFloat} = x !== zero(T) ? ldexp(one(T), exponent(x)) : x

ufp(x::Float64) = x !== 0.0 ? ldexp(1.0, exponent(x)) : 0.0
ufp(x::Float32) = x !== 0.0 ? ldexp(1.0f0, exponent(x)) : 0.0f0
ufp(x::Float16) = x !== 0.0 ? ldexp(one(Float16), exponent(x)) : zero(Float16)


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


if !hasmethod(prevfloat, (Float64, Int))
    Base.prevfloat(x::T, n::Int) where {T<:AbstractFloat} = -nextfloat(-x, n)
end

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
