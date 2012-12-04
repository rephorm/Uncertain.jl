Uncertain.jl
============

Handle error propogation in Julia

Example Usage
-------------

~~~julia

julia> load("uncertain.jl")
julia> import Uncertain.UncertainNumber

julia> a = UncertainNumber(1.0, 0.1)
UncertainNumber{Float64}(1.0,0.1)

julia> b = UncertainNumber(2.0, 0.2)
UncertainNumber{Float64}(2.0,0.2)

julia> 5*a
UncertainNumber{Float64}(5.0,0.5)

julia> a*2
UncertainNumber{Float64}(2.0,0.2)

julia> a*b
UncertainNumber{Float64}(2.0,0.28284271247461906)

julia> cosh(a^2-3*sinh(b*a))
UncertainNumber{Float64}(9773.540317776484,31230.955771261884)
~~~
