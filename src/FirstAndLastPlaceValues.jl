module FirstAndLastPlaceValues

export ufp, ulp, uls, rre, ε⁻, ε⁺, eta

#=
References

Accurate Floating-Point Summation Part I: Faithful Rounding
Siegfried M. Rump, Takeshi Ogita and Shin'ichi Oishi
http://oishi.info.waseda.ac.jp/~oishi/papers/RuOgOi07I.pdf

Fast quadruple-double floating point format
Naoya Yamanaka and Shin’ichi Oishi
2014-Jan-01

Error Estimation of Floating-Point Summation and Dot Product
Siegfried M. Rump

On the definition of ulp(x)
Jean-Michel Muller

Computing predecessor and successor in rounding to nearest
Siegfried Rump, Paul Zimmermann, Sylvie Boldo, Guillaume Melquiond

On various ways to split a floating-point number
Claude-Pierre Jeannerod, Jean-Michel Muller, Paul Zimmermann
=#

two(::Type{T}) where {T<:Real} = one(T) + one(T)
two(::Type{Float64}) = 2.0
two(::Type{Float32}) = 2.0f0
two(::Type{Float16}) = Float16(2.0)

#=
  Nicolas Fabiano, Jean-Michel Muller, Joris Picot.
  Algorithms for triple-word arithmetic.
  IEEE Transactions on Computers, Institute of Electrical and Electronics Engineers, 2019, 68 (11), pp.1573-1583.
  10.1109/TC.2019.2918451. hal-01869009v2
  Sylvie Boldo, Mioara Joldes, Jean-Michel Muller, Valentina Popescu.
  Formal Verification of a FloatingPoint Expansion Renormalization Algorithm.
  8th International Conference on Interactive Theorem Proving (ITP’2017), Sep 2017, Brasilia, Brazil. 
  hal-01512417f
  ACCURATE FLOATING-POINT SUMMATION PART I: FAITHFUL ROUNDING
  SIEGFRIED M. RUMP, TAKESHI OGITA, AND SHIN’ICHI OISHI
  in SIAM J. Sci. Comput. Volume 31, Issue 1, pp. 189-224 (2008)
=#

#=
  ufp(x)  unit in the first place
              the weight of the most significant bit
  ulp(x)  unit in the last place
              the weight of the least significant bit
  uls(x)  unit in the last significant place
              the weight of the rightmost nonzero bit
              the largest power of 2 that divides x
              the largest k s.t. x / 2^k is an integer
  
  reu(T)  roundoff error unit
             2^(-precision)
             ulp(1) / 2
             
  blp(x)  bit in the last place:              ulp(x) == 2.0^blp(x)
  bfp(x)  bit in the first place:             ufp(x) == 2.0^bfp(x)
  bls(x)  bit in the last significant place:  uls(x) == 2.0^bls(x)
=#
#=
julia> x = cbrt(17*pi)        #  3.765878158314341
julia> ufp(x),ulp(x),uls(x)   #  (2.0, 4.440892098500626e-16, 8.881784197001252e-16)
julia> bfp(x),blp(x),bls(x)   #  (1, -51, -50)
julia> x = 1/cbrt(17*pi)      #  0.26554231389355776
julia> ufp(x),ulp(x),uls(x)   #  (0.25, 5.551115123125783e-17, 4.440892098500626e-16)
julia> bfp(x),blp(x),bls(x)   #  (-2, -54, -51)
=#

#=
    The relative rounding error unit (reu) , 
    the distance from 1.0 to the next smaller floating-point number, is denoted by eps, and 
    the underflow unit by eta, that is the smallest positive (subnormal) floating-point number.
    
For IEEE 754 double precision we have eps = 2^(−53) and eta = 2^(−1074).
      The smallest positive normalized floating-point number is ets = (1/2)(eps^(-1))(eta) = (1/2)(2^(53-1074))
=#
#=
const BitPrecision = 1 + trailing_ones(Base.significand_mask(Float64))  # 53 = 1 implicit bit + 52 significand bits

const reu = 1.0 - prevfloat(1.0)  # 2.0^(-53)          # 1.0 - prevfloat(1.0)
const eta = nextfloat(0.0)        # 2.0^(-1074)        # smallest positive subnormal == nextfloat(0.0)
const ets = floatmin(Float64)     # 2.0^(53-1074) / 2  # smallest positive normal    == floatmin(Float64)
=#

explicit_precision(::Type{T}) where {T} = trailing_ones(Base.significand_mask(T))
implicit_precision(::Type{T}) where {T} = one(T) + explicit_precision(T)

reu(::Type{T}) where {T} = one(T) - prevfloat(one(T))
eta(::Type{T}) where {T} = nextfloat(zero(T))
ets(::Type{T}) where {T} = floatmin(T)

ufp(x) = 2.0^floor(Int,log2(abs(x)))
ulp(x::T) where {T} = ufp(x) * 2.0^(-explicit_precision(T)+1)
uls(x) = ulp(x) * 2.0^(trailing_zeros(reinterpret(UInt64,x)))

bfp(x) = round(Int, log2(ufp(x)))
blp(x) = round(Int, log2(ulp(x)))
bls(x) = round(Int, log2(uls(x)))

# for x <= sqrt(floatmax(T)) !! bounds need testing
fast_ufp(x::Float64) = reinterpret(Float64, reinterpret(UInt64, x) & 0xfff0000000000000)
fast_ufp(x::Float32) = reinterpret(Float32, reinterpret(UInt32, x) & 0xfff0_0000)
fast_ufp(x::Float16) = reinterpret(Float16, reinterpret(UInt16, x) & 0xfff0)

      

#=
eta is the smallest positive (nonzero, denormal)
eta(T) = nextfloat(zero(T))
=#
eta(::Type{T}) where {T<:AbstractFloat} = nextfloat(zero(T))

eta(::Type{Float64}) = 5.0e-324
eta(::Type{Float32}) = 1.0f-45
eta(::Type{Float16}) = Float16(6.0e-8)


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
ufp(x) unit in the first place of x

ufp(zero(T)) ≝ zero(T)

ufp(x::T) ≝ 2^floor(Int, log2(abs(x)))
          ≡ ldexp(one(T), floor(Int,log2(abs(x))))


the unsigned constants used in `ufp(x::T)`
     Float64: ((~UInt64(0)) >> 52) << 52
     Float32: ((~UInt32(0)) >> 23) << 23
     Float16: ((~UInt16(0)) >> 10) << 10

the shifts: Base.trailing_ones(Base.significand_mask)
=#

ulp(::Type{T}) where {T<:AbstractFloat} = two(T) * ε⁻(T)

ulp(::Type{Float64}) = 2.220446049250313e-16
ulp(::Type{Float32}) = 1.1920929f-7
ulp(::Type{Float16}) = Float16(0.000977)

ulp(x::T) where {T<:AbstractFloat} = two(T) * ε⁻(T) * ufp(x)
# ulp(x::T) where {T<:AbstractFloat} = ε⁺(T) * ufp(x)

ulp(x::Float64) = 2.220446049250313e-16 * ufp(x)
ulp(x::Float32) = 1.1920929f-7 * ufp(x)
ulp(x::Float16) = Float16(0.000977) * ufp(x)


end # FirstAndLastPlaceValues
