module Uncertain

export
  UncertainNumber


#
# Define the type
#
type UncertainNumber{T<:FloatingPoint} <: Real
    value::T
    uncertainty::T
end

UncertainNumber{T<:FloatingPoint}(v::T) = UncertainNumber{T}(v, zero(T))
UncertainNumber{T<:Real}(v::T) = UncertainNumber{FloatingPoint}(float(v), zero(FloatingPoint))
function UncertainNumber{T<:FloatingPoint, S<:FloatingPoint}(v::T, u::S)
    R = promote_type(T, S)
    UncertainNumber{R}(convert(R,v), convert(R, u))
end
UncertainNumber{T<:Real, S<:Real}(v::T, u::S) = UncertainNumber(float(v), float(u))


#
# Implement basic operations
#
import Base.+, Base.-, Base.*, Base./, Base.^

+(a::UncertainNumber, b::UncertainNumber) = UncertainNumber(a.value + b.value, hypot(a.uncertainty, b.uncertainty))

-(a::UncertainNumber, b::UncertainNumber) = UncertainNumber(a.value - b.value, hypot(a.uncertainty, b.uncertainty))

function *(a::UncertainNumber, b::UncertainNumber)
    val = a.value * b.value
    UncertainNumber(val, val * hypot(a.uncertainty/a.value,
                    b.uncertainty/b.value))
end

function /(a::UncertainNumber, b::UncertainNumber)
    val = a.value / b.value
    UncertainNumber(val, val * hypot(a.uncertainty/a.value,
                    b.uncertainty/b.value))
end

# XXX these seem to be ignored, with ^ induced from * instead, which is not correct
^(a::UncertainNumber, b::FloatingPoint) = UncertainNumber(a.value^b, b*a.value^(b-1)*a.uncertainty)

^(a::UncertainNumber, b::UncertainNumber) = UncertainNumber(a.value^b.value, hypot( b.value * a.value^(b.value-1) * a.uncertainty,
        a.value ^ b.value * log(a.value) * b.uncertainty) )

#
# Trig functions
#
import Math.sin, Math.cos, Math.tan, Math.cot, Math.sec, Math.csc,
       Math.asin, Math.acos, Math.atan, Math.acot, Math.asec, Math.acsc,
       Math.sinh, Math.cosh, Math.tanh, Math.coth, Math.sech, Math.csch,
       Math.asinh, Math.acosh, Math.atanh, Math.acoth, Math.asech, Math.acsch,
       Math.sind, Math.cosd, Math.tand, Math.cotd, Math.secd, Math.cscd,
       Math.asind, Math.acosd, Math.atand, Math.acotd, Math.asecd, Math.acscd

# TODO (this list contains other functions that Math exports, and hasn't been filtered at all):
#       atan2, radians2degrees, degrees2radians,
#       log, log2, log10, log1p, logb, exp, exp2, expm1, 
#       cbrt, sqrt, square, erf, erfc, ceil, floor, trunc, round, 
#       lgamma, hypot, gamma, lfact, max, min, ilogb, ldexp, frexp,
#       clamp, modf, ^, 
#       airy, airyai, airyprime, airyaiprime, airybi, airybiprime,
#       besselj0, besselj1, besselj, bessely0, bessely1, bessely,
#       hankelh1, hankelh2, besseli, besselk, besselh,
#       beta, lbeta, eta, zeta, psigamma, digamma

# Warning: these fail for zeros of df
for (f, df) in (
                (:sin, :cos),
                (:cos, :sin),
                (:tan, x->sec(x)^2),
                (:cot, x->csc(x)^2),
                (:sec, x->sec(x)*tan(x)),
                (:csc, x->csc(x)*cot(x)),

                (:asin, x->1.0/sqrt(1-x^2)),
                (:acos, x->1.0/sqrt(1-x^2)),
                (:atan, x->1.0/(1+x^2)),
                (:acot, x->1.0/(1+x^2)),
                (:asec, x->1.0/(x*sqrt(x^2-1))),
                (:acsc, x->1.0/(x*sqrt(x^2-1))),

                (:sinh, :cosh),
                (:cosh, :sinh),
                (:tanh, x->sech(x)^2),
                (:coth, x->csch(x)^2),
                (:sech, x->sech(x)*tanh(x)),
                (:csch, x->csch(x)*coth(x)),

                (:asinh, x->1.0/sqrt(x^2+1)),
                (:acosh, x->1.0/sqrt(x^2-1)),
                (:atanh, x->1.0/(1-x^2)),
                (:acoth, x->1.0/(1-x^2)),
                (:asech, x->1.0/(x*sqrt(1-x^2))),
                (:acsch, x->1.0/(x*sqrt(1+x^2))),

                (:sind, :cosd),
                (:cosd, :sind),
                (:tand, x->secd(x)^2),
                (:cotd, x->cscd(x)^2),
                (:secd, x->secd(x)*tand(x)),
                (:cscd, x->cscd(x)*cotd(x)),

                (:asind, x->1.0/sqrt(1-x^2)),
                (:acosd, x->1.0/sqrt(1-x^2)),
                (:atand, x->1.0/(1+x^2)),
                (:acotd, x->1.0/(1+x^2)),
                (:asecd, x->1.0/(x*sqrt(x^2-1))),
                (:acscd, x->1.0/(x*sqrt(x^2-1))),

                )
    @eval begin
        ($f)(y::UncertainNumber) = UncertainNumber($f(y.value), abs($df(y.value)) * y.uncertainty)
    end
end

#
# Handle conversions and promotions
#
import Base.promote_rule, Base.convert

convert{T<:FloatingPoint}(::Type{UncertainNumber{T}}, x::UncertainNumber) = UncertainNumber(convert(T, x.value), convert(T, x.uncertainty))
convert{T<:FloatingPoint}(::Type{UncertainNumber{T}}, x::Rational) = UncertainNumber(convert(T,x))
convert{T<:FloatingPoint}(::Type{UncertainNumber{T}}, x::Real) = UncertainNumber(x)

#convert{T<:FloatingPoint}(::Type{T}, x::UncertainNumber) = convert(T, x.value)

promote_rule{T<:FloatingPoint}(::Type{UncertainNumber{T}}, ::Type{T}) = UncertainNumber{T}
promote_rule{T<:FloatingPoint, S<:Real}(::Type{UncertainNumber{T}}, ::Type{S}) = UncertainNumber{promote_type(T,S)}
promote_rule{T<:FloatingPoint, S<:FloatingPoint}(::Type{UncertainNumber{T}}, ::Type{UncertainNumber{S}}) = UncertainNumber{promote_type(T,S)}


end
