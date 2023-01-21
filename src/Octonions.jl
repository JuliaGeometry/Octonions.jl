module Octonions

import Base: +, -, *, /, ^, ==
import Base: abs, abs2, conj, exp, inv, isreal, isfinite, isinf, iszero, isnan, log, real, sqrt
import Base: promote_rule, float
import Base: rand, randn
using Quaternions: Quaternion
using Random

Base.@irrational INV_SQRT_EIGHT 0.3535533905932737622004 sqrt(big(0.125))

include("octonion.jl")

export Octonion, OctonionF16, OctonionF32, OctonionF64
export imag_part, octo, octorand

end # module
