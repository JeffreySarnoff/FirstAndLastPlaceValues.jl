

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
  
