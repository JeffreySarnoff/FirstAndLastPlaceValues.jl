two(::Type{T}) where {T<:Number} = one{T} + one{T}

fracbits(::Type{T}) where {T<:Real} = precision(T) - 1

implicitbit(::Type{Float16})  = one(UInt16)  << precision(Float16)
implicitbit(::Type{Float32})  = one(UInt32)  << precision(Float32)
implicitbit(::Type{Float64})  = one(UInt64)  << precision(Float64)
implicitbit(::Type{Float128}) = one(UInt128) << precision(Float128)

significant_bits(::Type{T}) where {T} = precision(T)
significant_mask(::Type{T}) where {T} = Base.significand_mask(T) | implicitbit(T)

# twopow_eplicit(T) = 2^fracbits(T)  
# 2^(count of explicitly stored significand bits)
# 2^trailing_significand_bits
twopow_explicit(::Type{T}) where {T<:Real} = trunc(Int64, two(T)^fracbits(T))

function trailing_ones(Base.significand_mask(T)) where {T<:Real}

#=
  Float16  => 2^Base.significand_bits(Float16)
  Float32  => 2^Base.significand_bits(Float32)
  Float64  => 2^Base.significand_bits(Float64)
  Float128 => 2^Base.significand_bits(Float128)


const twopow_fracbits = Dict([ 
    Float16  => trailing_ones(Base.significand_mask(Float16)),
    Float32  => trailing_ones(Base.significand_mask(Float32)),
    Float64  => trailing_ones(Base.significand_mask(Float64)),
    Float128 => trailing_ones(Base.significand_mask(Float128)) ])

function ufp_abs(fp::T) where {T<:Reak}
      twopow_fractionbits =  # lookup table T => 2.0^(precision(T) - 1)
      fp_scaled = fp * twopow_fractionbits
      succ(fp_scaled) - fp_scaled
end                      

@inline function ufp(fp::T) where {T<:Real}
    ufp_abs(abs(fp))
end  
  
