using InstantiateFromURL
# optionally add arguments to force installation: instantiate = true, precompile = true
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.8.0")

using LinearAlgebra, Statistics
using Pkg
#Pkg.add("ForwardDiff")
#Pkg.add("Zygote")
#Pkg.add("Optim")
#Pkg.add("JuMP")
#Pkg.add("Ipopt")
#Pkg.add("BlackBoxOptim")
#Pkg.add("Roots")
#Pkg.add("NLsolve")
#Pkg.add("LeastSquaresOptim")
using ForwardDiff, Zygote, Optim, JuMP, Ipopt, BlackBoxOptim, Roots, NLsolve, LeastSquaresOptim
using Optim: converged, maximum, maximizer, minimizer, iterations #some extra functions

# Automatic differentiation
function f(x_1, x_2)
    w_1 = x_1
    w_2 = x_2
    w_3 = w_1 * w_2
    w_4 = sin(w_1)
    w_5 = w_3 + w_4
    return w_5
end

# forward-mode automatic differentation
using ForwardDiff
h(x) = sin(x[1]) + x[1] * x[2] + sinh(x[1] * x[2])
x = [1.4 2.2]
@show ForwardDiff.gradient(h,x) # use Automatic differentation, seeds from x

f(x) = sum(sin, x) + prod(tan, x) * sum(sqrt, x)
g = (x) -> ForwardDiff.gradient(f,x); # g() is now the gradient
g(rand(5)) # gradient at a random point

# AD complicated functions with embedded iterations
function squareroot(x)
    z = copy(x) # initial starting point for Newton method
    while abs(z*z -x) > 1e-13
        z = z - (z*z-x)/(2z)
    end
    return z
end

squareroot(2.0)


using ForwardDiff
dsqrt(x) = ForwardDiff.derivative(squareroot, x)
dsqrt(2.0)

# Zygote.jl
using Zygote

h(x, y) = 3x^2 + 2x + 1 + y*x - y
gradient(h, 3.0, 5.0)

D(f) = x-> gradient(f,x)[1] # returns first in tuple
D_sin = D(sin)
D_sin(4.0)

# for functions of one (julia) variable, find the derivative using the ' after a function name

using Statistics
p(x) = mean(abs, x)
p'([1.0, 3.0, -2.0])
squareroot'(2.0)

# zygote supports combinations of vectors and scalars as the function parameters
h(x,n) = (sum(x.^n))^(1/n)
@show gradient(h, [1.0, 4.0, 6.0], 2.0) # x is vector, n is scalar

using Optim, LinearAlgebra
N = 1000000
y = rand(N)
λ = 0.01
obj(x) = sum((x .- y).^2) + λ*norm(x)

x_iv = rand(N)
function g!(G, x)
    G .=  obj'(x)
end

results = optimize(obj, g!, x_iv, LBFGS()) # or ConjugateGradient()
println("minimum = $(results.minimum) with in "*
"$(results.iterations) iterations")


# Optimization

# univariate functions on bounded intervals
using Optim
using Optim: converged, maximum, maximizer, minimizer, iterations

result = optimize(x-> x^2, -2.0, 1.0)
# check if the results converbed, and throw errors otherwise
converged(result) || error("Failed")
xmin = result.minimizer
xmin

f(x) = -x^2
result = maximize(f, -2.0, 1.0)
converged(result) || error("Failed to converge in $(iterations(reuslt)) iterations")
xmin = maximizer(result)
fmax = maximum(result)

# unconstrained multivariate optimization
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x_iv = [0.0, 0.0]
results = optimize(f, x_iv) # i.e. optimize(f, x_iv, NelderMead())

# no derivative, used finite differences to approximate the gradient of f(x)
results = optimize(f, x_iv, LBFGS())
println("minimum = $(results.minimum) with argmin = $(results.minimizer)) in $(results.iterations) iterations ")

    # use autodifferentiation
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x_iv = [0.0, 0.0]
results = optimize(f, x_iv, LBFGS(), autodiff=:forward) # i.e. use ForwardDiff.jl
println("minimum = $(results.minimum) with argmin = $(results.minimizer) in $(results.iterations) iterations")

f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x_iv = [0.0, 0.0]
function g!(G, x)
    G[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
    G[2] = 200.0 * (x[2] - x[1]^2)
end

results = optimize(f, g!, x_iv, LBFGS()) # or ConjugateGradient()
println("minimum = $(results.minimum) with argmin = $(results.minimizer) in "*
"$(results.iterations) iterations")

results = optimize(f, x_iv, SimulatedAnnealing())




# JuMP.jl
# meaning of Ipopt: Interior Point OPTimizer - a nonlinear solver in Julia

using JuMP, Ipopt

# https://jump.dev/JuMP.jl/0.18/quickstart.html

# solve
# max(x[1] + x[2])
# st sqrt(x[1]^2 + x[2]^2) <=1

function squareroot(x) # pretending we don't know sqrt()
    z = x # Initial starting point for Newton’s method
    while abs(z*z - x) > 1e-13
        z = z - (z*z-x)/(2z)
    end
    return z
end


m = Model(with_optimizer(Ipopt.Optimizer))
# need to register user defined functions for AD
JuMP.register(m,:squareroot, 1, squareroot, autodiff=true)

@variable(m, x[1:2], start=0.5) # start is the initial condition
@objective(m, Max, sum(x))
@NLconstraint(m, squareroot(x[1]^2+x[2]^2) <= 1)
@show JuMP.optimize!(m)

    # example: quadratic objective
    # solve
    # min (1-x)^2 + 100(y-x^2)^2)
    # st x + y >= 10

using JuMP,Ipopt
m = Model(with_optimizer(Ipopt.Optimizer)) # settings for the solver
@variable(m, x, start = 0.0)
@variable(m, y, start = 0.0)

@NLobjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

JuMP.optimize!(m)
println("x = ", value(x), " y = ", value(y))

# adding a (linear) constraint
@constraint(m, x + y == 10)
JuMP.optimize!(m)
println("x = ", value(x), " y = ", value(y))


# simple example using JuMP

using JuMP
using Pkg
#Pkg.add("Clp")
using Clp
m = Model(Clp.Optimizer)
@variable(m, 0 <= x <= 2)
@variable(m, 0 <= y <= 30)

@objective(m, Max, 5*x + 3*y)
@constraint(m, 1x+5*y <=3.0)

print(m)

JuMP.optimize!(m)
println("x = ", value(x), " y = ", value(y))
@show value(x)



# BlackBoxOptim.jl

using BlackBoxOptim

function rosenbrock2d(x)
    return(1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

results = bboptimize(rosenbrock2d; SearchRange = (-5.0, 5.0), NumDimensions = 2);


# systems of equations and least squares

# Roots.jl
using Roots
f(x) = sin(4 * (x - 1/4)) + x + x^20 - 1
fzero(f, 0, 1)

# NLsolve.jl
    # solve for multivariate systems of equations and fixed points

using NLsolve

f(x) = [(x[1]+3)*(x[2]^3-7)+18
        sin(x[2]*exp(x[1])-1)] # returns an array

results = nlsolve(f, [ 0.1; 1.2])

results = nlsolve(f, [ 0.1; 1.2], autodiff=:forward)

println("converged=$(NLsolve.converged(results)) at root=$(results.zero) in $(results.iterations) iterations and $(results.f_calls) function calls")


# LeastSquaresOptim.jl
    # nonlinear least squares

# Jacobian matrix: first order partial derivatives
# Hessian matrix definition: a square matrix of second-order partial derivatives of a scalar-valued function, or scalar field. It describes the local curvature of a function of many variables.

using LeastSquaresOptim
function rosenbrock(x)
    [1 - x[1], 100 * (x[2]-x[1]^2)]
end
LeastSquaresOptim.optimize(rosenbrock, zeros(2), Dogleg())


function rosenbrock_f!(out, x)
    out[1] = 1 - x[1]
    out[2] = 100 * (x[2]-x[1]^2)
end
LeastSquaresOptim.optimize!(LeastSquaresProblem(x = zeros(2),
                                f! = rosenbrock_f!, output_length = 2))

# if you want to use gradient
function rosenbrock_g!(J, x)
    J[1, 1] = -1
    J[1, 2] = 0
    J[2, 1] = -200 * x[1]
    J[2, 2] = 100
end
LeastSquaresOptim.optimize!(LeastSquaresProblem(x = zeros(2), f! = rosenbrock_f!, g! = rosenbrock_g!, output_length = 2))


# exercise: autodifferentiation implementation
struct DualNumber{T} <: Real
    val::T
    ϵ::T
end

import Base: +, *, -, ^, exp
+(x::DualNumber, y::DualNumber) = DualNumber(x.val + y.val, x.ϵ + y.ϵ)
+(x::DualNumber, a::Number) = DualNumber(x.val + a, x.ϵ)
+(a::Number, x::DualNumber) = DualNumber(x.val + a, x.ϵ)


f(x, y) = 3.0 + x + y
x = DualNumber(2.0, 1.0) # x -> 2.0+ 1.0ϵ
y = DualNumber(3.0, 0.0) # y = 3.0, no derivative
# seeded calculates both teh function and d/dx gradient
f(x, y)


# Questions:
    # 1. define a function g!(), what's the meaning of !
    # 2. JuMP first example
    # 3. when to use struct?

# NOTES:
    # optimize, JuMP, NLsolve - all three optimization can use autodiff=:forward i.e. forward autodifferentiation
