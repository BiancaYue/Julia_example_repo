# Introductory Examples Exercise 2
# Using only rand(), write a function binomial_rv such that binomial_rv(n,p) generates one draw of Y~Bin（n, p)
# Hint: if U is uniform on (0,1) and p ∈ (o,1) then the expression U<p evaluates to true with probability p
using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.8.0")

# answer
function binomial_rv(n, p)
    count = 0
    U = rand(n)
    for i = 1:n
        if U[i] < p
            count += 1
        end
    end
    return count
end

for j = 1:25
    b = binomial_rv(10, 0.5)
    print("$b,")
end

# exervise 3
# Compute an approximation to π using Monte Carlo
n = 1000000
cnt = 0
for i = 1:n
    u, v = rand(2)
    d = sqrt((u - 0.5)^2 + (v - 0.5)^2)
    if d < 0.5
        cnt += 1
    end
end
print(count)
area_estimate = cnt / n
print(area_estimate * 4)


# exercise 4
# Flip an unbiased coin 10 times.
# If 3 consecutive heads occur one or more times within this sequence, pay one dollar.
# If not, pay nothing.
rand(10)
payoff = 0
count = 0
print("Count = ")
for i in 1:10
    U = rand()
    if U < 0.5
        count +=1
    else
        count = 0
    end
    print(count)
    if count == 3
        payoff = 1
    end
end
println("\npayoff = $payoff")

#exervise 6
# x[t+1]= α * x[t] +ϵ
# plot three simulated time series, one for each of the cases α=0, α = 0.8, α = 0.98
αs = [0.0, 0.8, 0.98, 1.0]
n = 200
p = plot()

for α in αs
    x = zeros(n + 1)
    x[1] = 0.0
    for t in 1:n
        x[t + 1] = α * x[t] + randn()
    end
    plot!(p, x, label = "alpha = $α") # add to plot p
end
p


# plot white noise
using Plots
gr(fmt=:png)
n = 1000
ϵ = randn(n)
plot(1:n, ϵ)

Pkg.add("Distributions")

using Distributions
function plothistogram(distribution, n)
    ϵ = rand(distribution, n)
    histogram(ϵ)
end
lp = Laplace()
plothistogram(lp, 500)
plothistogram(Normal(0,1),10000)

#composing packages
eps()






# Arrays, tuples and other types
# exercise 1
using LinearAlgebra


function compute_asymptotic_var(A, ∑; S0 = ∑ * ∑', tolerance = 1e-6, maxiter = 500)
    V = ∑ * ∑'
    S = S0
    err = tolerance +1
    i = 1
    while err > tolerance && i ≤ maxiter
        next_S = A * S * A' + V
        err = norm(S - next_S)
        S = next_S
        i +=1
    end
    return S
end


A = [0.8 -0.2; -0.1 0.7]
∑ = [0.5 0.4; 0.4 0.6]

maximum(abs, eigvals(A))

our_solution = compute_asymptotic_var(A, ∑)

Pkg.add("QuantEcon")

using QuantEcon

norm(our_solution - solve_discrete_lyapunov(A, ∑ * ∑'))


# good style
function fixedpointmap(f; iv, tolerance=1E-7, maxiter=1000)
    # setup the algorithm
    x_old = iv
    normdiff = Inf
    iter = 1
    while normdiff > tolerance && iter <= maxiter
        x_new = f(x_old) # use the passed in map
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    return (value = x_old, normdiff=normdiff, iter=iter) # A named tuple
end

# define a map and parameters
p = 1.0
β = 0.9
f(v) = p + β * v # note that p and β are used in the function!

sol = fixedpointmap(f, iv=0.8, tolerance=1.0E-8) # don't need to pass
println("Fixed point = $(sol.value), and |f(x) - x| = $(sol.normdiff) in $(sol.iter)"*
        " iterations")

x = 0.0:0.2:1.0
print(x)
typeof(x)
collect(0:3:6)
collect(1:1:5)


using Pkg
Pkg.add("Parameters")

using Parameters
paramgen = @with_kw (α = 0.1, β = 0.2)  # create named tuples with defaults

# creates named tuples, replacing defaults
@show paramgen()
print(typeof(paramgen()))
print(typeof(paramgen)) # calling without arguments gives all defaults
@show paramgen(α = 0.2)
@show paramgen(α = 0.2, β = 0.5);


subtypes(AbstractFloat)

@show typeof(zero)
