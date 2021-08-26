using InstantiateFromURL
# optionally add arguments to force installation: instantiate = true, precompile = true
github_project("QuantEcon/quantecon-notebooks-julia", version = "0.8.0")

# Constructing and Accessing a DataFrame
using LinearAlgebra, Statistics
using DataFrames
using RDatasets, DataFramesMeta, CategoricalArrays, Query, VegaLite
using GLM

using DataFrames, RDatasets

# note use of missing
commodities = ["crude", "gas", "gold", "silver"]
last_price = [4.2, 11.3, 12.1, missing]
df = DataFrame(commod = commodities, price = last_price)

df.price

DataFrames.describe(df)

nt = (commod = "nickel", price=5.1)
push!(df, nt)

nt = (t = 1, col1 = 3.0)
df2 =  DataFrame([nt])
push!(df2, (t = 2, col1 = 4.0))

df[!, :price]
df[!, :price] *= 2.0

allowmissing!(df2, :col1) # necessary to add in a for col1
push!(df2, (t=3, col1 = missing))
push!(df2, (t=4, col1 = 5.1))

@show mean(df2.col1)
@show mean(skipmissing(df2.col1))


# coalesce return the first value in the arguments which is not equal to missing, if any. otherwise, return missing
df2.col1 .=coalesce.(df2.col1, 0.0) # replace all missing with 0.0

# Manipulating and Transforming DataFrames
using DataFramesMeta
f(x) = x^2
df2 = @transform(df2, col2 = f.(:col1))

# categorical data
using CategoricalArrays
id = [1, 2, 3, 4]
y = ["old", "young", "young", "old"]
y = CategoricalArray(y)
df = DataFrame(id=id, y=y)

levels(df.y)

# visualization, querying, and Plots

x = 3.0
f(x) = x^2
g(x) = log(x)

@show g(f(x))
@show x |> f |> g; # pipes nest function calls

using Query
df = DataFrame(name=["John", "Sally", "Kirk"], age=[23., 42., 59.], children=[3,5,2])

x = @from i in df begin
    @where i.age>50
    @select {i.name, i.children}
    @collect DataFrame
end

using RDatasets, VegaLite
iris = dataset("datasets", "iris")

iris |> @vlplot(
    :point,
    x=:PetalLength,
    y=:PetalWidth,
    color=:Species
)

iris
using GLM

x = randn(100)
y = 0.9 .* x + 0.5 * rand(100)
df = DataFrame(x=x, y=y)
ols = lm(@formula(y~x), df) #R-style notation

using RegressionTables
regtable(ols)
# regtable(ols,  renderSettings = latexOutput()) # for LaTex output

using FixedEffectModels
cigar = dataset("plm", "Cigar")
cigar.StateCategorical =  categorical(cigar.State)
cigar.YearCategorical =  categorical(cigar.Year)
fixedeffectresults = reg(cigar, @formula(Sales ~ NDI + fe(StateCategorical) + fe(YearCategorical)),
                            weights = :Pop, Vcov.cluster(:State))
regtable(fixedeffectresults)
