# multiple dispatch
x = range(0.0, 1.0, length = 20)
x_2 = 1:1:20 # if integers
@show typeof(x)
@show typeof(x_2)
@show supertypes(typeof(x))

using Pkg
Pkg.add("BenchmarkTools")

using BenchmarkTools
N = 3
A = rand(N, N)
x = rand(N)

@btime $A * $x
@btime inv($A)

#exercise 4
using Pkg
Pkg.add("Polynomials")

using Polynomials
p = Polynomial([2, -5, 2], :x)
@show p
p_prime = derivative(p)
@show p(0.1), p_prime(0.1)
@show roots(p)




# BAD
function g2(x::AbstractFloat)
    if x > 0.0   # can't efficiently call with `x::Integer`
        return x + 1.0   # OK - assumes you can promote `Float64` to `AbstractFloat`
    otherwise
        return 0   # BAD! Returns a `Int64`
    end
end

x = BigFloat(1.0)
x2 = BigFloat(-1.0)

@show typeof(g2(x))
@show typeof(g2(x2))  # type unstable

# GOOD
function g3(x) #
    if x > zero(x)   # any type with an additive identity
        return x + one(x)  # more general but less important of a change
    otherwise
        return zero(x)
    end
end

@show typeof(g3(x))
@show typeof(g3(x2));  # type stable
