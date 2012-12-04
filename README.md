Uncertain.jl
============

Handle error propogation in Julia

NB: See Caveats below!

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

julia> cosh(a^2 * b)
UncertainNumber{Float64}(3.7621956910836314,1.2563812996702035)
~~~

Caveats
-------

Currently, only basic arithmetic functions (+-*/) and trig/hypertrig functions are supported. Feel free to implement more and send a pull request.

Correlations are currently ignored, so, e.g. a+a gives a smaller uncertainty than 2*a.

The ^ operator is also currently incorrect. It looks like it is being induced from * instead of using the one defined in the module.

