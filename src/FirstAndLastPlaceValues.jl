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
beta == base == 2

ufp(x::Real)  unit in the first place
              ufp(0) := 0                              iszero(x)
              ufp(x) := 2^floor(log2(abs(x)))          !iszero(x)

              ufp(x) is the weight of the most significant bit in the representation of x.
                                          first nonzero `digit` of x in its base beta representation.

              !iszero(x):  ufp(x) <= abs(x) < 2 * ufp(x)
              this definition is indpendent of some floating-point format.

      ufp(fp::T) where T = signbit(fp) ? ufp_abs(abs(fp)) : ufp_abs(fp)
      function ufp_abs(fp::T) where T
            twopow_fractionbits =  # lookup table T => 2.0^(precision(T) - 1)
            fp_scaled = fp * twopow_fractionbits
            succ(fp_scaled) - fp_scaled
      end                      

      function ufp(fp, precision)
          scaled_fp = abs(fp) * 2^(precision - 1)    # ldexp(precision - 1, abs(fp))
          succ(scaled_fp) - scaled_fp
      end

      function ufp(fp, precision)
          scaled_fp = abs(fp) * 2^(precision - 1)    # ldexp(precision - 1, abs(fp))
          succ(scaled_fp) - scaled_fp
      end

      The smallest postive subnormal floating-point number ris 2^Emin, denoted by subrealmin
      The smallest positive normalized floating-point number is eta / (2 * eps)

      const predecessor_one = 1 - subrealmin  # prevfloat(one(T))
                                              # = 1 - 2^(-precision)
      const phi = 2^(precision-1) + 1         # 2^fracbits + 1
        Emin <= -1 < p <= Emax, |fp| < 2^(Emax - p + 1), roundToZero OR roundDown
        ufp(x) == ufp(|x|)        fl(ufp(fp)) = ufp(fp)
        2^(precision - 1) * ufp(fp) <= ufp(fp) <= 2^precision* ufp(fp)

      function ufp(fp, p, 2)
         q = phi * abs(fp)
         q - predecessor_one * q      #   ufp(fp) = q * (1 - predecessor_one)
                                      #   ufp(fp) = abs(fp) * phi * (1 - predecessor_one)
                                      #   ufp(fp) / abs(fp) = phi * (1 - predecessor_one)
                                      #   ufp(fp) / abs(fp) = (2^fracbits + 1) * (1 - predecessor_one)
                                      #   ufp(fp) / abs(fp) = (1 + 2^fracbits) * (1 - predecessor_one)
                                      #   ufp(fp) / abs(fp) = (1 - predecessor_one) + (2^fracbits * (1 - predecessor_one))
                                      #   ufp(fp) / abs(fp) = (1 - predecessor_one) + (2^fracbits - (2^fracbits * predecessor_one)))

# Unit in the rst place (ufp) in precision-p base-b page 5

function ufp(fp)
    feature(setround,0)            # rounding RoundToZero
    p  = precision(fp)             # precision p
                                   # flbeta(m, e) genertes m*2^e
    p1 = 1 - subrealmin(flbeta)    # predecessor of 1
    phi = flbeta(1,p-1) + 1        # beta^(p-1) + 1
    q = phi*abs(f)                 # result in the normalized range
    q - p1*q                       # ufp(fp)
end

function ulp(fp) in roundUp
    fp = abs(fp)           # where abs(fp) < mrealmax, te largest representable floating-point number
    succ(fp) - fp          # fp + subrealmin  # succesor(abs(fp)) == nextfloat(abs(fp))
end
     for fp >0,  ulp(fp) = succ(fp) - fp      
     for fp <0,  ulp(fp) = pred(fp) - fp
     for fp !=0, ulp(fp) = succ(abs(fp)) - abs(fp)
                 ulp(fp) = aulp(abs(fp))            aulp(fp) = succ(fp) - fp

      eta = 1 / (2 * eps) 
      The distance from 1.0 to the next smaller floating-point number, is denoted by eps
      eps = 1 / (2 * eta)

      For IEEE 754 double precision we have eps = 2^(−53) and eta = 2^(−1074).


ulp(x)  unit in the last place
            the weight of the least significant bit

        from Muller on the definition of ulp(x)
        Conclusion (definition 5)
          If x is a real number that lies between two finite consecutive FP numbers
          a and b, without being equal to one of them, then ulp (x) = |b - a|, 
          otherwise ulp (x) is the distance between the two finite FP numbers nearest x. 
          Moreover, ulp (NaN) is NaN.

          x::Real              fl(x) != x
          a::MachineFloat      isfinite(a)    a < x
          b::MachineFloat      isfinite(b)    x < b
                                              a < x < b
                                              of all fl(_) a, b are two machine floats nearest to x
                                              R(a) preceeds x , R(b) suceeds x
                                              R(a) < x < R(b)
                                              a <= fl(x) <= b
                                                                a == fl(x) && fl(x) != b
                                                                a != fl(x) && fl(x) == b

          Define L as the largest finite FP number, and L- as its predecessor
          If x is larger than L then KahanUlp(x) = L - L-
                                    | X - x | < (1/2) KahanUlp(x)  ==>  X = roundNearest(x)
          HOWEVER IEEE 754 says 
            an infinitely precise reesult with magnitude at least
            2^emax * (2 - 2^(-n))
            shall round to infinity with no change in sign.
            With this converntion, if X is finite, in radix 2
                X = roundNearesst(c) ==> |X - x| <= (1/2)KahanUlp(x)
                                         |X - x| <  (1/2)HarrisonUlp(x) ==> X = roundNearest(x)
                HarrisonUlp(x) = b - a where a != b and a,b are the closest straddling points (a < x < b)
                                (assumes unbounded exponent range)
                KahanUlp (x) is the width of the interval whose endpoints are the two finite representable
                    numbers nearest x (even if x is not contained within that interval).




          and !exists m::MachineFloat such that  a < m < b
        

uls(x)  unit in the last significant place
              the weight of the rightmost nonzero bit
              the largest power of 2 that divides x
              the largest k s.t. x / 2^k is an integer
  
blp(x)  bit in the last place:              ulp(x) == 2.0^blp(x)
bfp(x)  bit in the first place:             ufp(x) == 2.0^bfp(x)
bls(x)  bit in the last significant place:  uls(x) == 2.0^bls(x)

the relative rounding error unit (u)
define u as half the distance between 1 and its successor
u = 2^(1-precision) / 2 = 2^(1 - precision - 1) = 2^(1-1 - precision) = 2^(-precision)

    reu(T)  roundoff error unit === eps(T)
    c is c[orrectly rounded]

    where c is finite nonnegative and nextfloat(c) is finite and c in F
      succ(1) = 1 + 2u
      succ(c) >= max( c(1+2u), c+eta )
      pred(c) >= min( c(1-2u), c-eta )

    where c is finite negative and prevfloat(c) finite and c in F
      pred(1) = 1 - u
      pred(c) >= min( c(1+2u), c+eta )
      succ(c) >= max( c(1-2u), c-eta )

    The relative rounding error unit (reu) , 

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
fast_ufp(x::Float64) = reinterpret(Float64, reinterpret(UInt64, x) & 0xfff0_000_000_000_000)
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
