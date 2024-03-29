using LinearAlgebra
using Octonions
using Quaternions: Quaternions, Quaternion, QuaternionF64
using Random
using RealDot: realdot
using Test

_octo(a::Real) = octo(a)
function _octo(c::Complex{T}) where T
    z = zero(T)
    return octo(reim(c)..., z, z, z, z, z, z)
end
function _octo(q::Quaternion{T}) where T
    z = zero(T)
    return octo(real(q), Quaternions.imag_part(q)..., z, z, z, z)
end

@testset "Octonion" begin
    @testset "type aliases" begin
        @test OctonionF16 === Octonion{Float16}
        @test OctonionF32 === Octonion{Float32}
        @test OctonionF64 === Octonion{Float64}
    end

    @testset "Constructors" begin
        @testset "from coefficients" begin
            cs = [(1, 2.0, 3.0f0, 4//1, 5, 6, 7, 8), (1//1, 2.0f0, 3.0f0, 4, 5, 6, 7, 8)]
            @testset for coef in cs, T in (Float32, Float64, Int)
                q = @inferred Octonion{T}(coef...)
                @test q isa Octonion{T}
                @test q === Octonion{T}(convert.(T, coef)...)
                q2 = @inferred Octonion(convert.(T, coef)...)
                @test Octonion(convert.(T, coef)...) === q
            end
        end
        @testset "from real" begin
            @testset for x in (-1//1, 1.0, 2.0), T in (Float32, Float64, Int, Rational{Int})
                coef = T.((x, zeros(7)...))
                @test @inferred(Octonion{T}(x)) === Octonion{T}(coef...)
                @test @inferred(Octonion(T(x))) === Octonion{T}(coef...)
            end
        end
        @testset "from octonion" begin
            os = (
                Octonion(1, 2, 3, 4, 5, 6, 7, 8), OctonionF64(0, 1, 0, 0, 0, 0, 0, 0)
            )
            @testset for o in os, T in (Float32, Float64)
                coef = T.((o.s, o.v1, o.v2, o.v3, o.v4, o.v5, o.v6, o.v7))
                @test @inferred(Octonion{T}(o)) === Octonion{T}(coef...)
                @test @inferred(Octonion(o)) === o
            end
        end
    end

    @testset "==" begin
        @test Octonion(1.0, 2, 3, 4, 5, 6, 7, 8) == Octonion(1, 2, 3, 4, 5, 6, 7, 8)
        @test Octonion(1.0, 2, 3, 4, 5, 6, 7, 8) != Octonion(1, 2, 3, 4, 1, 2, 3, 4)
    end

    @testset "isequal" begin
        @test isequal(Octonion(1:8...), Octonion(1.0:8.0...))
        x = ntuple(identity, 8)
        o = Octonion(x...)
        for i in 1:8
            x2 = Base.setindex(x, 0, i)
            o2 = Octonion(x2...)
            @test !isequal(o, o2)
        end
        o = Octonion(NaN, -0.0, Inf, -Inf, 0.0, NaN, Inf, -Inf)
        @test isequal(o, o)
        @test !isequal(o, Octonion(NaN, 0.0, Inf, -Inf, 0.0, NaN, Inf, -Inf))
        @test !isequal(o, Octonion(NaN, -0.0, Inf, -Inf, -0.0, NaN, Inf, -Inf))
    end

    @testset "convert" begin
        @test convert(Octonion{Float64}, 1) === Octonion(1.0)
        @test convert(Octonion{Float64}, Octonion(1, 2, 3, 4, 5, 6, 7, 8)) ===
            Octonion(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
        @test convert(Octonion{Float64}, Octonion(0, 0, 0, 0, 1, 0, 0, 0)) ===
            Octonion(0.0, 0.0, 0.0, 0.0, 1, 0, 0, 0)
    end

    @testset "promote" begin
        @test promote(Octonion(1.0, 2:8...), 1.0) === (Octonion(1.0, 2:8...), Octonion(1.0))
        @test promote(Octonion(1.0f0, 2:8...), 2.0) ===
            (Octonion(1.0, 2:8...), Octonion(2.0))
        @test promote(Octonion(1.0f0), Octonion(2.0)) === (Octonion(1.0), Octonion(2.0))

        @test Octonion(1) == 1.0
    end

    @testset "shorthands" begin
        @test octo(1) === Octonion(1)
        @test octo(1, 2, 3, 4, 5, 6, 7, 8) === Octonion(1, 2, 3, 4, 5, 6, 7, 8)
        @test octo(Octonion(1, 2, 3, 4, 5, 6, 7, 8)) === Octonion(1, 2, 3, 4, 5, 6, 7, 8)
    end

    @testset "random generation" begin
        @testset "octorand" begin
            @test_deprecated octorand()
            o = octorand()
            @test o isa Octonion
        end

        @testset "rand($H)" for H in (OctonionF32, OctonionF64)
            rng = Random.MersenneTwister(42)
            o = rand(rng, H)
            @test o isa H

            os = rand(rng, H, 1000)
            @test eltype(os) === H
            @test length(os) == 1000
            xs = map(os) do o
                return [real(o); imag_part(o)...]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0.5; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 12; atol=0.01))
        end

        @testset "randn($H)" for H in (OctonionF32, OctonionF64)
            rng = Random.MersenneTwister(42)
            o = randn(rng, H)
            @test o isa H

            os = randn(rng, H, 10000)
            @test eltype(os) === H
            @test length(os) == 10000
            xs = map(os) do o
                return [real(o); imag_part(o)...]
            end
            xs_mean = sum(xs) / length(xs)
            xs_var = sum(x -> abs2.(x .- xs_mean), xs) / (length(xs) - 1)
            @test all(isapprox.(xs_mean, 0; atol=0.1))
            @test all(isapprox.(xs_var, 1 / 8; atol=0.1))
        end
    end

    @testset "basic" begin
        q = randn(OctonionF64)
        qnorm = sign(q)
        @test real(q) === q.s
        @test_throws MethodError imag(q)
        @test imag_part(q) === (q.v1, q.v2, q.v3, q.v4, q.v5, q.v6, q.v7)
        @test conj(q) ===
            Octonion(q.s, -q.v1, -q.v2, -q.v3, -q.v4, -q.v5, -q.v6, -q.v7)
        @test conj(qnorm) === Octonion(
            qnorm.s,
            -qnorm.v1,
            -qnorm.v2,
            -qnorm.v3,
            -qnorm.v4,
            -qnorm.v5,
            -qnorm.v6,
            -qnorm.v7,
        )
        @test conj(conj(q)) === q
        @test conj(conj(qnorm)) === qnorm
        @test float(Octonion(1:8...)) === Octonion(1.0:8.0...)
        @test Octonions.abs_imag(q) ≈
            abs(Octonion(0, q.v1, q.v2, q.v3, q.v4, q.v5, q.v6, q.v7))
    end

    @testset "abs/abs_imag don't over/underflow" begin
        for x in [1e-300, 1e300, -1e-300, -1e300]
            for i in 1:8
                z = Base.setindex(ntuple(zero, 8), x, i)
                o = octo(z...)
                @test abs(o) == abs(x)
                i == 1 || @test Octonions.abs_imag(o) == abs(x)
            end
        end
        @test isnan(abs(octo(fill(NaN, 8)...)))
        @test abs(octo(NaN, Inf, fill(NaN, 6)...)) == Inf
        @test abs(octo(-Inf, fill(NaN, 7)...)) == Inf
        @test abs(octo(0.0)) == 0.0
        @test abs(octo(Inf)) == Inf
        @test abs(octo(1, -Inf, 3:8...)) == Inf
        @test isnan(Octonions.abs_imag(octo(0, fill(NaN, 7)...)))
        @test Octonions.abs_imag(octo(0, Inf, fill(NaN, 6)...)) == Inf
        @test Octonions.abs_imag(octo(0, NaN, -Inf, fill(NaN, 5)...)) == Inf
        @test Octonions.abs_imag(octo(0.0)) == 0.0
        @test Octonions.abs_imag(octo(0.0, 0.0, Inf, fill(0, 5)...)) == Inf
        @test Octonions.abs_imag(octo(0, 1, -Inf, 3:7...)) == Inf
    end

    @testset "algebraic properties" begin
        for _ in 1:100, T in (Float32, Float64, Int32, Int64)
            if T <: Integer
                q, q1, q2, q3 = [octo(rand((-T(100)):T(100), 8)...) for _ in 1:4]
                c1, c2 = [complex(rand((-T(100)):T(100), 2)...) for _ in 1:2]
            else
                q, q1, q2, q3 = randn(Octonion{T}, 4)
                c1, c2 = randn(Complex{T}, 2)
            end

            # skewfield
            test_group(q1, q2, q3, +, zero(q), -)
            # test all group properties but associativity
            test_neutral(q1, one(q1), *)
            test_inverse(q1, one(q1), *, inv)
            test_multiplicative(q1, q2, *, norm)

            # complex embedding
            test_multiplicative(c1, c2, *, _octo)
            test_multiplicative(c1, c2, +, _octo)
        end
    end

    @testset "inv does not under/overflow" begin
        x = 1e-300
        y = inv(x)
        for i in 1:8
            z = zeros(8)
            z[i] = x
            z2 = vcat(0.0, fill(-0.0, 7))
            z2[i] = i == 1 ? y : -y
            o = octo(z...)
            @test inv(octo(z...)) == octo(z2...)

            z[i] = y
            z2[i] = i == 1 ? x : -x
            @test inv(octo(z...)) == octo(z2...)
        end
        @test isequal(inv(octo(-Inf, 1, -2, 3, -4, 5, -6, 7)), octo(-0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0))
        @test isequal(inv(octo(1, -2, Inf, 3, -4, 5, -6, 7)), octo(0.0, 0.0, -0.0, -0.0, 0.0, -0.0, 0.0, -0.0))
    end

    @testset "isreal" begin
        @test isreal(octo(1))
        @test !isreal(octo(1, 1, 0, 0, 0, 0, 0, 0))
        @test !isreal(octo(1, 0, 1, 0, 0, 0, 0, 0))
        @test !isreal(octo(1, 0, 0, 1, 0, 0, 0, 0))
        @test !isreal(octo(1, 0, 0, 0, 1, 0, 0, 0))
        @test !isreal(octo(1, 0, 0, 0, 0, 1, 0, 0))
        @test !isreal(octo(1, 0, 0, 0, 0, 0, 1, 0))
        @test !isreal(octo(1, 0, 0, 0, 0, 0, 0, 1))
    end

    @testset "iszero" begin
        @test iszero(octo(0))
        @test !iszero(octo(1))
        @test !iszero(octo(0, 1, 0, 0, 0, 0, 0, 0))
        @test !iszero(octo(0, 0, 1, 0, 0, 0, 0, 0))
        @test !iszero(octo(0, 0, 0, 1, 0, 0, 0, 0))
        @test !iszero(octo(0, 0, 0, 0, 1, 0, 0, 0))
        @test !iszero(octo(0, 0, 0, 0, 0, 1, 0, 0))
        @test !iszero(octo(0, 0, 0, 0, 0, 0, 1, 0))
        @test !iszero(octo(0, 0, 0, 0, 0, 0, 0, 1))
    end

    @testset "isone" begin
        @test isone(octo(1))
        @test !isone(octo(-1))
        @test !isone(octo(0, 1, 0, 0, 0, 0, 0, 0))
        @test !isone(octo(1, 1, 0, 0, 0, 0, 0, 0))
        @test !isone(octo(1, 0, 1, 0, 0, 0, 0, 0))
        @test !isone(octo(1, 0, 0, 1, 0, 0, 0, 0))
        @test !isone(octo(1, 0, 0, 0, 1, 0, 0, 0))
        @test !isone(octo(1, 0, 0, 0, 0, 1, 0, 0))
        @test !isone(octo(1, 0, 0, 0, 0, 0, 1, 0))
        @test !isone(octo(1, 0, 0, 0, 0, 0, 0, 1))
    end

    @testset "isfinite" begin
        @test isfinite(octo(1:8...))
        for val in (Inf, -Inf, NaN)
            @test !isfinite(octo(val, 0, 0, 0, 0, 0, 0, 0))
            @test !isfinite(octo(0, val, 0, 0, 0, 0, 0, 0))
            @test !isfinite(octo(0, 0, val, 0, 0, 0, 0, 0))
            @test !isfinite(octo(0, 0, 0, val, 0, 0, 0, 0))
            @test !isfinite(octo(0, 0, 0, 0, val, 0, 0, 0))
            @test !isfinite(octo(0, 0, 0, 0, 0, val, 0, 0))
            @test !isfinite(octo(0, 0, 0, 0, 0, 0, val, 0))
            @test !isfinite(octo(0, 0, 0, 0, 0, 0, 0, val))
        end
    end

    @testset "isinf" begin
        @test !isinf(octo(1, 2, 3, 4, 5, 6, 7, 8))
        @test !isinf(octo(1, 2, 3, 4, 5, 6, 7, NaN))
        for inf in (Inf, -Inf)
            @test isinf(octo(inf, 0, 0, 0, 0, 0, 0, 0))
            @test isinf(octo(0, inf, 0, 0, 0, 0, 0, 0))
            @test isinf(octo(0, 0, inf, 0, 0, 0, 0, 0))
            @test isinf(octo(0, 0, 0, inf, 0, 0, 0, 0))
            @test isinf(octo(0, 0, 0, 0, inf, 0, 0, 0))
            @test isinf(octo(0, 0, 0, 0, 0, inf, 0, 0))
            @test isinf(octo(0, 0, 0, 0, 0, 0, inf, 0))
            @test isinf(octo(0, 0, 0, 0, 0, 0, 0, inf))
        end
    end

    @testset "isnan" begin
        @test !isnan(octo(1, 2, 3, 4, 5, 6, 7, 8))
        @test !isnan(octo(1, 2, 3, 4, 5, 6, 7, Inf))
        @test !isnan(octo(1, 2, 3, 4, 5, 6, 7, -Inf))
        @test isnan(octo(NaN, 2, 3, 4, 5, 6, 7, 8))
        @test isnan(octo(1, NaN, 3, 4, 5, 6, 7, 8))
        @test isnan(octo(1, 2, NaN, 4, 5, 6, 7, 8))
        @test isnan(octo(1, 2, 3, NaN, 5, 6, 7, 8))
        @test isnan(octo(1, 2, 3, 4, NaN, 6, 7, 8))
        @test isnan(octo(1, 2, 3, 4, 5, NaN, 7, 8))
        @test isnan(octo(1, 2, 3, 4, 5, 6, NaN, 8))
        @test isnan(octo(1, 2, 3, 4, 5, 6, 7, NaN))
    end

    @testset "*" begin
        # verify basic correctness
        q0 = Octonion(1,0,0,0,0,0,0,0)
        q1 = Octonion(0,1,0,0,0,0,0,0)
        q2 = Octonion(0,0,1,0,0,0,0,0)
        q3 = Octonion(0,0,0,1,0,0,0,0)
        q4 = Octonion(0,0,0,0,1,0,0,0)
        q5 = Octonion(0,0,0,0,0,1,0,0)
        q6 = Octonion(0,0,0,0,0,0,1,0)
        q7 = Octonion(0,0,0,0,0,0,0,1)
        @test q0 * q0 == q0
        @test q0 * q1 == q1
        @test q0 * q2 == q2
        @test q0 * q3 == q3
        @test q0 * q4 == q4
        @test q0 * q5 == q5
        @test q0 * q6 == q6
        @test q0 * q7 == q7
        @test q1 * q0 == q1
        @test q1 * q1 == -q0
        @test q1 * q2 == q3
        @test q1 * q3 == -q2
        @test q1 * q4 == -q7
        @test q1 * q5 == -q6
        @test q1 * q6 == q5
        @test q1 * q7 == q4
        @test q2 * q0 == q2
        @test q2 * q1 == -q3
        @test q2 * q2 == -q0
        @test q2 * q3 == q1
        @test q2 * q4 == q6
        @test q2 * q5 == -q7
        @test q2 * q6 == -q4
        @test q2 * q7 == q5
        @test q3 * q0 == q3
        @test q3 * q1 == q2
        @test q3 * q2 == -q1
        @test q3 * q3 == -q0
        @test q3 * q4 == -q5
        @test q3 * q5 == q4
        @test q3 * q6 == -q7
        @test q3 * q7 == q6
        @test q4 * q0 == q4
        @test q4 * q1 == q7
        @test q4 * q2 == -q6
        @test q4 * q3 == q5
        @test q4 * q4 == -q0
        @test q4 * q5 == -q3
        @test q4 * q6 == q2
        @test q4 * q7 == -q1
        @test q5 * q0 == q5
        @test q5 * q1 == q6
        @test q5 * q2 == q7
        @test q5 * q3 == -q4
        @test q5 * q4 == q3
        @test q5 * q5 == -q0
        @test q5 * q6 == -q1
        @test q5 * q7 == -q2
        @test q6 * q0 == q6
        @test q6 * q1 == -q5
        @test q6 * q2 == q4
        @test q6 * q3 == q7
        @test q6 * q4 == -q2
        @test q6 * q5 == q1
        @test q6 * q6 == -q0
        @test q6 * q7 == -q3
        @test q7 * q0 == q7
        @test q7 * q1 == -q4
        @test q7 * q2 == -q5
        @test q7 * q3 == -q6
        @test q7 * q4 == q1
        @test q7 * q5 == q2
        @test q7 * q6 == q3
        @test q7 * q7 == -q0

        @testset "* same between quaternion and Octonion" begin
            q1 = Octonion(1,0,0,0,0,0,0,0)
            qi = Octonion(0,1,0,0,0,0,0,0)
            qj = Octonion(0,0,1,0,0,0,0,0)
            qk = Octonion(0,0,0,1,0,0,0,0)
            @test q1 * q1 == q1
            @test q1 * qi == qi
            @test q1 * qj == qj
            @test q1 * qk == qk
            @test qi * q1 == qi
            @test qi * qi == -q1
            @test qi * qj == qk
            @test qi * qk == -qj
            @test qj * q1 == qj
            @test qj * qi == -qk
            @test qj * qj == -q1
            @test qj * qk == qi
            @test qk * q1 == qk
            @test qk * qi == qj
            @test qk * qj == -qi
            @test qk * qk == -q1
        end
    end

    @testset "abs2" for _ in 1:100, T in (Float16, Float32, Float64)
        o = rand(Octonion{T})
        @test abs2(o) == o'*o
    end

    @testset "/" begin
        for _ in 1:100
            o, o2 = randn(OctonionF64, 2)
            x = randn()
            @test o / o ≈ o \ o ≈ one(o)
            @test o / o2 ≈ o * inv(o2)
            @test o2 \ o ≈ inv(o2) * o
            @test o / x ≈ x \ o ≈ inv(x) * o
        end

        @testset "no overflow/underflow" begin
            @testset for x in [1e-300, 1e300, -1e-300, -1e300]
                @test octo(x) / octo(x) == octo(1)
                @testset for i in 2:8
                    z = Base.setindex(ntuple(zero, 8), x, i)
                    z2 = Base.setindex(ntuple(zero, 7), -1, i - 1)
                    @test octo(x) / octo(z...) == octo(0, z2...)
                end
                @test octo(0, x, zeros(6)...) / octo(x, 0, 0, 0, zeros(4)...) == octo(0, 1, 0, 0, zeros(4)...)
                @test octo(0, x, zeros(6)...) / octo(0, x, 0, 0, zeros(4)...) == octo(1, 0, 0, 0, zeros(4)...)
                @test octo(0, x, zeros(6)...) / octo(0, 0, x, 0, zeros(4)...) == octo(0, 0, 0, -1, zeros(4)...)
                @test octo(0, x, zeros(6)...) / octo(0, 0, 0, x, zeros(4)...) == octo(0, 0, 1, 0, zeros(4)...)
            end
            @testset for T in [Float32, Float64]
                o = one(T)
                z = zero(T)
                inf = T(Inf)
                nan = T(NaN)
                @testset for s in [1, -1], t in [1, -1]
                    @test isequal(octo(o) / octo(s*inf), octo(s*z, fill(-z, 7)...))
                    @test isequal(octo(o) / octo(s*inf, t*o, z, t*z, z, t*z, z, t*z), octo(s*z, -t*z, -z, -t*z, -z, -t*z, -z, -t*z))
                    @test isequal(octo(o) / octo(s*inf, t*nan, t*z, z, t*z, z, t*z, z), octo(s*z, nan, -t*z, -z, -t*z, -z, -t*z, -z))
                    @test isequal(octo(o) / octo(s*inf, t*inf, t*z, z, t*z, z, t*z, z), octo(s*z, -t*z, -t*z, -z, -t*z, -z, -t*z, -z))
                end
                @test isequal(octo(inf) / octo(inf, 1:7...), octo(fill(nan, 8)...))
                @test isequal(octo(inf) / octo(inf, 1, 2, -inf, 4:7...), octo(fill(nan, 8)...))
            end
        end
    end

    @testset "^" begin
        @testset "^(::Octonion, ::Real)" begin
            for _ in 1:100
                o = randn(OctonionF64)
                @test @inferred(o^2.0) ≈ o * o
                @test o^1.0 ≈ o
                @test o^-1.0 ≈ inv(o)
                @test o^1.3 ≈ exp(1.3 * log(o))
                @test o^7.8 ≈ exp(7.8 * log(o))
                @test o^1.3f0 ≈ exp(1.3f0 * log(o))
                @test o^7.8f0 ≈ exp(7.8f0 * log(o))
            end
        end
        @testset "^(::Octonion, ::Octonion)" begin
            @test octo(Float64(ℯ))^octo(0, 0, 0, 0, 0, 0, π / 2, 0) ≈
                octo(0, 0, 0, 0, 0, 0, 1, 0)
            z = (3.5 + 2.3im)^(0.2 + 1.7im)
            @test octo(3.5, 0, 0, 0, 0, 0, 2.3, 0)^octo(0.2, 0, 0, 0, 0, 0, 1.7, 0) ≈
                octo(real(z), 0, 0, 0, 0, 0, imag(z), 0)
            for _ in 1:100
                q, p = randn(OctonionF64, 2)
                @test @inferred(q^p) ≈ exp(p * log(q))
            end
        end
    end

    @testset "non-analytic functions" begin
        unary_funs = [conj, abs, abs2, norm, sign]
        # since every octonion is conjugate to a quaternion,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            for _ in 1:100
                o1, o2 = randn(OctonionF64, 2)
                q = randn(QuaternionF64)
                o = _octo(q)
                @test @inferred(fun(o)) ≈ _octo(fun(q))
                @test o2 * fun(o1) * inv(o2) ≈ fun(o2 * o1 * inv(o2))
            end
        end
    end

    @testset "analytic functions" begin
        # all complex analytic functions can be extended to the octonions
        #! format: off
        unary_funs = [
            sqrt, inv, exp, exp2, exp10, expm1, log, log2, log10, log1p,
            sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh,
            csc, sec, cot, acsc, asec, acot, csch, sech, coth, acsch, asech, acoth,
            sinpi, cospi,
        ]
        #! format: on
        # since every octonion is conjugate to a quaternion,
        # one can establish correctness as follows:
        @testset for fun in unary_funs
            for _ in 1:100
                o1, o2 = randn(OctonionF64, 2)
                q = randn(QuaternionF64)
                o = _octo(q)
                @test @inferred(fun(o)) ≈ _octo(fun(q))
                @test o2 * fun(o1) * inv(o2) ≈ fun(o2 * o1 * inv(o2))
            end
        end

        @testset "identities" begin
            for _ in 1:100
                o = randn(OctonionF64)
                @test inv(o) * o ≈ o * inv(o) ≈ one(o)
                @test sqrt(o) * sqrt(o) ≈ o
                @test exp(log(o)) ≈ o
                @test exp(zero(o)) === one(o)
                @test log(one(o)) === zero(o)
                @test exp2(log2(o)) ≈ o
                @test exp10(log10(o)) ≈ o
                @test expm1(log1p(o)) ≈ o
                @test sinpi(o) ≈ sin(π * o)
                @test cospi(o) ≈ cos(π * o)
                @test all(sincos(o) .≈ (sin(o), cos(o)))
                @test all(sincos(zero(o)) .≈ (sin(zero(o)), cos(zero(o))))
                if VERSION ≥ v"1.6"
                    @test all(sincospi(o) .≈ (sinpi(o), cospi(o)))
                    @test all(sincospi(zero(o)) .≈ (sinpi(zero(o)), cospi(zero(o))))
                end
                @test tan(o) ≈ cos(o) \ sin(o) ≈ sin(o) / cos(o)
                @test tanh(o) ≈ cosh(o) \ sinh(o) ≈ sinh(o) / cosh(o)
                @testset for (f, finv) in [
                    (sin, csc),
                    (cos, sec),
                    (tan, cot),
                    (sinh, csch),
                    (cosh, sech),
                    (tanh, coth),
                ]
                    @test f(o) ≈ inv(finv(o))
                end
                @testset for (f, finv) in [
                    (asin, acsc),
                    (acos, asec),
                    (atan, acot),
                    (asinh, acsch),
                    (acosh, asech),
                    (atanh, acoth),
                ]
                    @test f(o) ≈ finv(inv(o))
                end
            end
        end

        @testset "additional properties" begin
            @testset "log" begin
                @test log(zero(OctonionF64)) === octo(-Inf)
                @test log(one(OctonionF64)) === octo(0.0)
                @test log(-one(OctonionF64)) ≈ _octo(log(complex(-1.0)))
                x = rand()
                @test log(octo(x)) ≈ octo(log(x))
                @test log(octo(-x)) ≈ _octo(log(complex(-x)))
            end

            @testset "exp" begin
                @test exp(octo(0)) === octo(1.0)
                @test exp(octo(2)) === octo(exp(2))
                @test norm(exp(octo(0))) ≈ 1
                @test norm(exp(octo(2))) ≠ 1
                @test exp(octo(0.0)) === octo(1.0)
                for i in 2:8
                    z = setindex!(zeros(8), 2, i)
                    z2 = setindex!(zeros(8), sin(2), i)
                    @test exp(octo(z...)) === octo(cos(2), z2[2:end]...)
                    @test norm(exp(octo(z...))) ≈ 1
                    @test exp(octo(2.0)) === octo(exp(2))
                end
                @test exp(octo(0)) isa OctonionF64
                @test exp(octo(0.0)) isa OctonionF64
                @test exp(octo(0//1)) isa OctonionF64
                @test exp(octo(BigFloat(0))) isa Octonion{BigFloat}
                @test exp(octo(fill(1, 8)...)) ≈ exp(octo(fill(1.0, 8)...))
            end
        end
    end

    @testset "sign" begin
        for _ in 1:100
            q = randn(OctonionF64)
            qnorm = @inferred sign(q)
            @test abs(qnorm) ≈ 1
            @test q ≈ abs(q) * qnorm
            @test sign(qnorm) ≈ qnorm
        end
        @inferred(sign(octo(1:8...)))
    end

    @testset "RealDot with $T" for T in (Float32, Float64)
        for _ in 1:10
            q1 = randn(Octonion{T})
            q2 = randn(Octonion{T})
            # Check real∘dot is equal to realdot.
            @test real(dot(q1,q2)) == @inferred(realdot(q1,q2))
            # Check realdot is commutative.
            @test realdot(q1,q2) == realdot(q2,q1)
            # Check real∘dot is also commutative just in case.
            @test real(dot(q1,q2)) == real(dot(q2,q1))
            # Check the return type of realdot is correct.
            @test realdot(q1,q2) isa T
        end
    end
end
