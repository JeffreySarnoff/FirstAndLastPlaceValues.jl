module FirstAndLastPlaceValues

export ufp, ulp, rre, ε⁻, ε⁺, eta

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
=#

two(::Type{T}) where {T<:Real} = one(T) + one(T)
two(::Type{Float64}) = 2.0
two(::Type{Float32}) = 2.0f0
two(::Type{Float16}) = Float16(2.0)


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


end # FirstAndLastPlaceValues
