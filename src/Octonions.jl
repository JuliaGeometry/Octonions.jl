module Octonions

using Random
using RealDot: RealDot

Base.@irrational INV_SQRT_EIGHT 0.3535533905932737622004 sqrt(big(0.125))

include("octonion.jl")

export Octonion, OctonionF16, OctonionF32, OctonionF64
export imag_part, octo, octorand

end # module
