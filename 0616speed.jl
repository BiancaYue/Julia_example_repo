using InstantiateFromURL
# optionally add arguments to force installation: instantiate = true, precompile = true
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.8.0")
using LinearAlgebra, Statistics

# what is multiple dispatch
# Answer: When an operation like addition is requested, the Julia compiler inspects the type of data to be acted on and hands it out to the appropriate method.
# This process is called multiple dispatch.

+(1,1)

@which +(1.0 , 1.0)

# Adding methods

# + method, add integer and string
import Base: +
+(x::Integer, y::String) = x + parse(Int,y)
+(100, "10")
100 + "215"

function q(x)  # or q(x::Any)
    println("Default (Any) method invoked")
end

function q(x::Number)
    println("Number method invoked")
end

function q(x::Integer)
    println("Integer method invoked")
end

q("str")
q(3)

q(2+1im)

f(x) = x > 0.0 ? x : 0.0

f(-1)
@code_warntype f(1)


# summary and tips
#Use functions to segregate operations into logically distinct blocks.
# Data types will be determined at function boundaries.
# If types are not supplied then they will be inferred.
# If types are stable and can be inferred effectively your functions will run fast.
