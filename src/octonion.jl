struct Octonion{T<:Real} <: Number
  s::T
  v1::T
  v2::T
  v3::T
  v4::T
  v5::T
  v6::T
  v7::T
end

Octonion{T}(x::Real) where {T<:Real} = Octonion(convert(T, x))
Octonion{T}(o::Octonion) where {T<:Real} =
  Octonion{T}(o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7)

Octonion(s::Real, v1::Real, v2::Real, v3::Real, v4::Real, v5::Real, v6::Real, v7::Real) =
  Octonion(promote(s, v1, v2, v3, v4, v5, v6, v7)...)
Octonion(x::Real) = Octonion(x, zero(x), zero(x), zero(x), zero(x), zero(x), zero(x), zero(x))

const OctonionF16 = Octonion{Float16}
const OctonionF32 = Octonion{Float32}
const OctonionF64 = Octonion{Float64}

Base.promote_rule(::Type{Octonion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}
Base.promote_rule(::Type{Octonion{T}}, ::Type{Octonion{S}}) where {T <: Real, S <: Real} = Octonion{promote_type(T, S)}

octo(p, v1, v2, v3, v4, v5, v6, v7) = Octonion(p, v1, v2, v3, v4, v5, v6, v7)
octo(x) = Octonion(x)

Base.real(o::Octonion) = o.s
imag_part(o::Octonion) = (o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7)

Base.:/(o::Octonion, x::Real) = Octonion(o.s / x, o.v1 / x, o.v2 / x, o.v3 / x, o.v4 / x, o.v5 / x, o.v6 / x, o.v7 / x)
Base.:*(o::Octonion, x::Real) = Octonion(o.s * x, o.v1 * x, o.v2 * x, o.v3 * x, o.v4 * x, o.v5 * x, o.v6 * x, o.v7 * x)
Base.:*(x::Real, o::Octonion) = o * x

Base.conj(o::Octonion) = Octonion(o.s, -o.v1, -o.v2, -o.v3, -o.v4, -o.v5, -o.v6, -o.v7)
Base.abs(o::Octonion) = _hypot((o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7))
Base.float(q::Octonion{T}) where T = convert(Octonion{float(T)}, q)
abs_imag(o::Octonion) = _hypot((o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7))
Base.abs2(o::Octonion) = RealDot.realdot(o, o)
function Base.inv(o::Octonion)
    if isinf(o)
        return octo(
            copysign(zero(o.s), o.s),
            flipsign(-zero(o.v1), o.v1),
            flipsign(-zero(o.v2), o.v2),
            flipsign(-zero(o.v3), o.v3),
            flipsign(-zero(o.v4), o.v4),
            flipsign(-zero(o.v5), o.v5),
            flipsign(-zero(o.v6), o.v6),
            flipsign(-zero(o.v7), o.v7),
        )
    end
    a = max(abs(o.s), maximum(abs, imag_part(o)))
    p = o / a
    io = conj(p) / (a * abs2(p))
    return io
end

Base.isreal(o::Octonion) = iszero(o.v1) & iszero(o.v2) & iszero(o.v3) & iszero(o.v4) & iszero(o.v5) & iszero(o.v6) & iszero(o.v7)
Base.isfinite(o::Octonion) = isfinite(real(o)) & isfinite(o.v1) & isfinite(o.v2) & isfinite(o.v3) & isfinite(o.v4) & isfinite(o.v5) & isfinite(o.v6) & isfinite(o.v7)
Base.iszero(o::Octonion) = iszero(real(o)) & iszero(o.v1) & iszero(o.v2) & iszero(o.v3) & iszero(o.v4) & iszero(o.v5) & iszero(o.v6) & iszero(o.v7)
Base.isnan(o::Octonion) = isnan(real(o)) | isnan(o.v1) | isnan(o.v2) | isnan(o.v3) | isnan(o.v4) | isnan(o.v5) | isnan(o.v6) | isnan(o.v7)
Base.isinf(o::Octonion) = isinf(real(o)) | isinf(o.v1) | isinf(o.v2) | isinf(o.v3) | isinf(o.v4) | isinf(o.v5) | isinf(o.v6) | isinf(o.v7)

Base.:-(o::Octonion) = Octonion(-o.s, -o.v1, -o.v2, -o.v3, -o.v4, -o.v5, -o.v6, -o.v7)

Base.:+(o::Octonion, w::Octonion) = Octonion(o.s + w.s,
                                            o.v1 + w.v1,
                                            o.v2 + w.v2,
                                            o.v3 + w.v3,
                                            o.v4 + w.v4,
                                            o.v5 + w.v5,
                                            o.v6 + w.v6,
                                            o.v7 + w.v7)

Base.:-(o::Octonion, w::Octonion) = Octonion(o.s - w.s,
                                            o.v1 - w.v1,
                                            o.v2 - w.v2,
                                            o.v3 - w.v3,
                                            o.v4 - w.v4,
                                            o.v5 - w.v5,
                                            o.v6 - w.v6,
                                            o.v7 - w.v7)

function Base.:*(o::Octonion, w::Octonion)
    s  = ((o.s * w.s - o.v4 * w.v4) - (o.v2 * w.v2 + o.v6 * w.v6)) - ((o.v1 * w.v1 + o.v5 * w.v5) + (o.v3 * w.v3 + o.v7 * w.v7))
    v1 = ((o.s * w.v1 + o.v1 * w.s) + (o.v6 * w.v5 - o.v5 * w.v6)) + ((o.v2 * w.v3 - o.v3 * w.v2) + (o.v7 * w.v4 - o.v4 * w.v7))
    v2 = ((o.s * w.v2 + o.v2 * w.s) + (o.v4 * w.v6 - o.v6 * w.v4)) + ((o.v3 * w.v1 - o.v1 * w.v3) + (o.v7 * w.v5 - o.v5 * w.v7))
    v3 = ((o.s * w.v3 + o.v3 * w.s) + (o.v5 * w.v4 - o.v4 * w.v5)) + ((o.v1 * w.v2 - o.v2 * w.v1) + (o.v7 * w.v6 - o.v6 * w.v7))
    v4 = ((o.s * w.v4 + o.v4 * w.s) + (o.v3 * w.v5 - o.v5 * w.v3)) + ((o.v1 * w.v7 - o.v7 * w.v1) + (o.v6 * w.v2 - o.v2 * w.v6))
    v5 = ((o.s * w.v5 + o.v5 * w.s) + (o.v2 * w.v7 - o.v7 * w.v2)) + ((o.v1 * w.v6 - o.v6 * w.v1) + (o.v4 * w.v3 - o.v3 * w.v4))
    v6 = ((o.s * w.v6 + o.v6 * w.s) + (o.v3 * w.v7 - o.v7 * w.v3)) + ((o.v2 * w.v4 - o.v4 * w.v2) + (o.v5 * w.v1 - o.v1 * w.v5))
    v7 = ((o.s * w.v7 + o.v7 * w.s) + (o.v5 * w.v2 - o.v2 * w.v5)) + ((o.v4 * w.v1 - o.v1 * w.v4) + (o.v6 * w.v3 - o.v3 * w.v6))
    return Octonion(s, v1, v2, v3, v4, v5, v6, v7)
end

function Base.:/(o::Octonion{T}, w::Octonion{T}) where T
    # handle over/underflow while matching the behavior of /(a::Complex, b::Complex)
    a = max(abs(real(w)), maximum(abs, imag_part(w)))
    if isinf(w)
        if isfinite(o)
            return octo(
                zero(T)*sign(o.s)*sign(w.s),
                -zero(T)*sign(o.v1)*sign(w.v1),
                -zero(T)*sign(o.v2)*sign(w.v2),
                -zero(T)*sign(o.v3)*sign(w.v3),
                -zero(T)*sign(o.v4)*sign(w.v4),
                -zero(T)*sign(o.v5)*sign(w.v5),
                -zero(T)*sign(o.v6)*sign(w.v6),
                -zero(T)*sign(o.v7)*sign(w.v7),
            )
        end
        return octo(T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN))
    end
    p = w / a
    return (o * conj(p)) / RealDot.realdot(w, p)
end

Base.:(==)(q::Octonion, w::Octonion) = (q.s == w.s) & (q.v1 == w.v1) & (q.v2 == w.v2) & (q.v3 == w.v3) &
                                 (q.v4 == w.v4) & (q.v5 == w.v5) & (q.v6 == w.v6) & (q.v7 == w.v7)
function Base.isequal(q::Octonion, w::Octonion)
    return (isequal(q.s, w.s) & isequal(q.v1, w.v1) & isequal(q.v2, w.v2) &
            isequal(q.v3, w.v3) & isequal(q.v4, w.v4) & isequal(q.v5, w.v5) &
            isequal(q.v6, w.v6) & isequal(q.v7, w.v7))
end

"""
    extend_analytic(f, o::Octonion)

Evaluate the extension of the complex analytic function `f` to the octonions at `o`.

Given ``o = s + a u``, where ``s`` is the real part, ``u`` is a pure unit octonion,
and ``a \\ge 0`` is the magnitude of the imaginary part of ``o``,
```math
f(o) = \\Re(f(z)) + \\Im(f(z)) u,
```
is the extension of `f` to the octonions, where ``z = s + a i`` is a complex analog to
``o``.

See [^DentoniSce1973] and [^ColomboSabadini2020] for details.

[^DentoniSce1973]: Dentoni, P. and Sce M. "Funzioni regolari nell'algebra di Cayley."
                   Rendiconti del Seminario matematico della Università di Padova 50 (1973): 251-267.
                   Translation: [^ColomboSabadini2020]
[^ColomboSabadini2020]: Colombo, F., Sabadini, I., Struppa, D.C. (2020).
                        Regular Functions in the Cayley Algebra.
                        In: Michele Sce's Works in Hypercomplex Analysis.
                        doi: [10.1007/978-3-030-50216-4_6](https://doi.org/10.1007/978-3-030-50216-4_6)
"""
function extend_analytic(f, o::Octonion)
    # Adapted from Quaternions.jl
    a = abs_imag(o)
    s = o.s
    z = complex(s, a)
    w = f(z)
    wr, wi = reim(w)
    scale = wi / a
    # o == real(o), so f(real(o)) may be real or complex, i.e. wi may be nonzero.
    # we choose to embed complex numbers in the octonions by identifying the first
    # imaginary octonion basis with the complex imaginary basis.
    wi_octo = a > 0 ? map(x -> x * scale, imag_part(o)) : ntuple(_ -> zero(scale), Val(7))
    return octo(wr, wi_octo...)
end

for f in (:sqrt, :exp, :exp2, :exp10, :expm1, :log2, :log10, :log1p,
    :sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh,
    :csc, :sec, :cot, :acsc, :asec, :acot, :csch, :sech, :coth, :acsch, :asech, :acoth,
    :sinpi, :cospi,
)
    @eval Base.$f(o::Octonion) = extend_analytic($f, o)
end

for f in (@static(VERSION ≥ v"1.6" ? (:sincos, :sincospi) : (:sincos,)))
    @eval begin
        function Base.$f(o::Octonion)
            a = abs_imag(o)
            z = complex(o.s, a)
            s, c = $f(z)
            sr, si = reim(s)
            cr, ci = reim(c)
            sscale = si / a
            cscale = ci / a
            ov = imag_part(o)
            si_octo = a > 0 ? map(x -> x * sscale, ov) : ntuple(_ -> zero(sscale), Val(7))
            ci_octo = a > 0 ? map(x -> x * cscale, ov) : si_octo
            return octo(sr, si_octo...), octo(cr, ci_octo...)
        end
    end
end

function Base.log(o::Octonion)
  a = abs(o)
  o = o / a
  s = o.s
  M = abs_imag(o)
  th = atan(M, s)
  if M > 0
    M = th / M
    return Octonion(log(a),
                     o.v1 * M,
                     o.v2 * M,
                     o.v3 * M,
                     o.v4 * M,
                     o.v5 * M,
                     o.v6 * M,
                     o.v7 * M)
  else
    z = zero(th)
    return Octonion(log(a), ifelse(iszero(a), z, th), z, z, z, z, z, z)
  end
end

Base.:^(o::Octonion, w::Octonion) = exp(w * log(o))

octorand(rng::AbstractRNG = Random.GLOBAL_RNG) = octo(randn(rng), randn(rng), randn(rng), randn(rng), randn(rng), randn(rng), randn(rng), randn(rng))

function Base.rand(rng::AbstractRNG, ::Random.SamplerType{Octonion{T}}) where {T<:Real}
  Octonion{T}(rand(rng, T), rand(rng, T), rand(rng, T), rand(rng, T),
              rand(rng, T), rand(rng, T), rand(rng, T), rand(rng, T))
end

function Base.randn(rng::AbstractRNG, ::Type{Octonion{T}}) where {T<:AbstractFloat}
  Octonion{T}(
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
      randn(rng, T) * INV_SQRT_EIGHT,
  )
end

function RealDot.realdot(o::Octonion, w::Octonion)
    return ((o.s * w.s + o.v4 * w.v4) + (o.v2 * w.v2 + o.v6 * w.v6)) +
           ((o.v1 * w.v1 + o.v5 * w.v5) + (o.v3 * w.v3 + o.v7 * w.v7))
end

# copied from https://github.com/JuliaLang/julia/blob/v1.9.0-beta3/base/math.jl
# to work around 3+arg hypot being slow on <v1.9
# https://github.com/JuliaLang/julia/issues/44336
function _hypot(x::NTuple{N,<:Number}) where {N}
    maxabs = maximum(abs, x)
    if isnan(maxabs) && any(isinf, x)
        return typeof(maxabs)(Inf)
    elseif (iszero(maxabs) || isinf(maxabs))
        return maxabs
    else
        return maxabs * sqrt(sum(y -> abs2(y / maxabs), x))
    end
end
