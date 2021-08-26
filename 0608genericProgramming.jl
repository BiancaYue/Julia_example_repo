using Pkg
Pkg.add("InstantiateFromURL")


using InstantiateFromURL
# optionally add arguments to force installation: instantiate = true, precompile = true
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.8.0")

using Distributions
x = 1
y = Normal()
z = "foo"
@show x, y, z
@show typeof(x), typeof(y), typeof(z)
@show supertype(typeof(x))
typeof(x)|>supertype # |> (x,f) applies a function to the preceding argument, same as supertypeof(typeof(x))

using Base: show_supertypes
show_supertypes(Int64)

# custom type
struct MyType
    a::Float64
end

myval= MyType(2.0)
@show myval
@show typeof(myval)
@show supertype(typeof(myval))



# Distributions

using Distributions
d1 = Normal(1.0, 2.0)
@show d1
show_supertypes(typeof(d1))
rand(d1)

# simulate a stochastic process
function simulateprocess(x0; a = 1.0, b = 1.0, N = 5, d::Sampleable{Univariate,Continuous})
    x = zeros(typeof(x0), N+1)
    for t in 2:N+1
        x[t] = a*x[t-1] + b * rand(d)
    end
    return x
end
@show simulateprocess(0.0, d = Normal(0.5,4))

d1 = Normal(1.0, 2.0)
d2 = Exponential(0.1)
@show d1
@show d2
@show supertype(typeof(d1))
@show supertype(typeof(d2))
@show pdf(d1, 0.1)
@show pdf(d2, 0.1)
@show cdf(d1, 0.1)
@show cdf(d2, 0.1)
@show support(d1)
@show support(d2)
@show minimum(d1)
@show minimum(d2)
@show maximum(d1)
@show maximum(d2);

using Pkg
Pkg.add("StatsPlots")
using StatsPlots
d = Normal(2.0, 1.0)
plot(d)


# complex numbers

x = 4.0 + 1.0im
@show x, typeof(x)

xBig = BigFloat(4.0) + 1.0im
@show xBig, typeof(xBig)
y = Complex(4.0, 1.0)
y == x

@which abs(x)


Pkg.add("Interpolations")
using Interpolations

x = 0.0:0.2:1.0
f(x) = x^2
f_int = LinearInterpolation(x, f.(x))
@show f_int(1.0)

function plotfunctions(f)
    intf(x) = quadgk(f, 0.0, x)[1]  # int_0^x f(x) dx

    x = 0:0.1:1.0
    f_x = f.(x)
    plot(x, f_x, label="f")
    plot!(x, intf.(x), label="int_f")
end

Pkg.add("QuadGK")
using QuadGK

plotfunctions(f_int)
