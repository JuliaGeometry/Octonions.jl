module Octonions

using Random
using RealDot: RealDot

include("octonion.jl")

export Octonion, OctonionF16, OctonionF32, OctonionF64
export imag_part, octo, octorand

end # module
