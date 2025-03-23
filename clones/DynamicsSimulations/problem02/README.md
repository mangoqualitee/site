---
title: "problem-02"
author: "Vishal Paudel"
date: "2025/01/24"
---

> 2\. $\quad$ 3D. Get good at vectors. Assume that the positions relative to an origin of four random points, which are randomly located in space are given as $\vec{\mathbf{r}}_A$, $\vec{\mathbf{r}}_B$, $\vec{\mathbf{r}}_C$ and $\vec{\mathbf{r}}_D$. Assume force $\vec{\mathbf{F}}$ is given. For each problem below write a single vector formula (one for each problem) that answers the question.  
> 
>     a) The points A and B define an infinite line. So do the points C and D. Find the distance between these two lines. ‘The’ distance means ‘the minimum distance’, that is the length of the shortest line segment connecting the two lines. Either write a formula (or sequence of formulas), or write computer code that gives the answer, or both.  
>     b) Same problem as above, but also find the end points of the shortest line segment.  
>     c) Find the volume of the tetrahedron ABCD (you should reason-out and not quote any formulas for the volume of a tetrahedron, that is, see if you can derive the formula: ’volume = one third base times height’).  
>     d) Assume points A, B and C are fixed to a structure. All three are connected, by massless rods, to a ball and socket at each end, to point D. At point D the force  F is applied. Find the tension in bar AD. Find a formula for the answer, or write computer code to find the answer, or both. The goal is to find a formula for the tension in terms of the positions and the force vector.

# Write formula or code for the length of the shortest line segment

The lines can be represented as:

$$l_1: \vec{p}_1(\lambda_1) = \quad \vec{\mathbf{r}}_A + \lambda_1 (\vec{\mathbf{r}}_B-\vec{\mathbf{r}}_A)$$

$$l_2: \vec{p}_2(\lambda_2) = \quad \vec{\mathbf{r}}_C + \lambda_2 (\vec{\mathbf{r}}_D-\vec{\mathbf{r}}_C)$$

Hence, for two points on either lines, $p_1$, $p_2$, the distance is:

$$D(\lambda_1, \lambda_2) = \quad \left|\left| \mathbf{\vec{p}}_2 - \mathbf{\vec{p}}_1 \right|\right|_n$$

Which can be framed as an optimisation problem, minimum length, $\hat{s}$:

$$\hat{s} = \quad \min_{\lambda_1, \lambda_2} || \left( \vec{\mathbf{r}}_A + \lambda_1 (\vec{\mathbf{r}}_B-\vec{\mathbf{r}}_A) - (\vec{\mathbf{r}}_C + \lambda_2 (\vec{\mathbf{r}}_D-\vec{\mathbf{r}}_C)) \right) ||_n$$

or,

$$\hat{s} \leftarrow D\left\{ \nabla D(\lambda_1, \lambda_2) = 0 \right\}$$

The code corresponding to this in file [./MinimumDistanceBetweenTwoLines/src/MinimumDistanceBetweenTwoLines.jl](./MinimumDistanceBetweenTwoLines/src/MinimumDistanceBetweenTwoLines.jl) is:

```julia
module MinimumDistanceBetweenTwoLines

import Symbolics
import LinearAlgebra

# problem setup
const N = 3
Symbolics.@variables r⃗ₛ[1:N] r⃗ₜ[1:N] r⃗ᵤ[1:N] r⃗ᵥ[1:N]
Symbolics.@variables λ₁ λ₂

p⃗₁ = r⃗ₛ + λ₁ * (r⃗ₜ - r⃗ₛ)
p⃗₂ = r⃗ᵤ + λ₂ * (r⃗ᵥ - r⃗ᵤ)
s⃗ = p⃗₂ - p⃗₁
D² = LinearAlgebra.dot(s⃗, s⃗)
D² = Symbolics.scalarize(D²)
∇D² = Symbolics.gradient(D², [λ₁, λ₂])
eq = ∇D² .~ 0

lambdas_symbolic = Dict([λ₁, λ₂] .=> Symbolics.symbolic_linear_solve(eq, [λ₁, λ₂]))
p̂₁_symbolic, p̂₂_symbolic = Symbolics.substitute.([p⃗₁, p⃗₂], Ref(lambdas_symbolic))
ŝ_symbolic = Symbolics.substitute.(D², Ref(lambdas_symbolic))

#...
```

Note that in the code, we are optimising for distance squared, and the $s$ symbols in the code represent the distance squared.

In the same file, we take specific values:

```julia
#...

# evaluation
subs = Dict(
    r⃗ₛ => [0.0, 0.0, 0.0],
    r⃗ₜ => [1.0, 0.0, 0.0],
    r⃗ᵤ => [0.0, 0.0, 1.0],
    r⃗ᵥ => [0.0, 1.0, 2.0],
    # λ₁ => 3,
    # λ₂ => 4,
)
∇D²_eval = Symbolics.substitute.(∇D², Ref(subs))

eq_eval = Symbolics.substitute.(eq, Ref(subs))

result_eval = Dict(symbol => Symbolics.substitute(expresion, subs) for (symbol, expresion) in lambdas_symbolic)
merged_eval = merge(subs, result_eval)

ŝ_eval = Symbolics.substitute.(D², Ref(merged_eval))
p̂₁_eval, p̂₂_eval = Symbolics.substitute.([p⃗₁, p⃗₂], Ref(merged_eval))

println("=== Substituted results: =====")
println("Λ̂  : ", result_eval)
println("ŝ  = ", ŝ_eval)
println("p̂₁ = ", p̂₁_eval)
println("p̂₂ = ", p̂₂_eval)

println("\n=== Symbolic results: =====")
println("Λ̂  : ", lambdas_symbolic)
println("ŝ  = ", Symbolics.scalarize(ŝ_symbolic))
println("p̂₁ = ", Symbolics.scalarize(p̂₁_symbolic))
println("p̂̂₂ = ", Symbolics.scalarize(p̂₂_symbolic))

end # module MinimumDistanceBetweenTwoLines
```

This produces the following:

```julia
julia> include("src/MinimumDistanceBetweenTwoLines.jl")
WARNING: replacing module MinimumDistanceBetweenTwoLines.
=== Substituted results: =====
Λ̂  : Dict{Symbolics.Num, Float64}(λ₂ => -0.5, λ₁ => -0.0)
ŝ  = 0.5
p̂₁ = [0.0, 0.0, 0.0][Base.OneTo(3)]
p̂₂ = [0.0, -0.5, 0.5][Base.OneTo(3)]

=== Symbolic results: =====
Λ̂  : Dict{Symbolics.Num, SymbolicUtils.BasicSymbolic{Real}}(λ₂ => ((-(2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(⃗rₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3])) / ((-((-2(-⃗rᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)), λ₁ => ((-((-(2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) + 2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3])) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)))
ŝ  = (r⃗ᵤ[1] - r⃗ₛ[1] + ((-r⃗ᵤ[1] + r⃗ᵥ[1])*(((-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) + (((((-(-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - ⃗rₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) + 2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(r⃗ₛ[1] - r⃗ₜ[1])) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)))^2 + (r⃗ᵤ[2] - r⃗ₛ[2] + (((((-(-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) + 2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(r⃗ₛ[2] - r⃗ₜ[2])) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + ((((-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-r⃗ᵤ[2] + r⃗ᵥ[2])) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)))^2 + (r⃗ᵤ[3] - r⃗ₛ[3] + ((((((-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) - 2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-r⃗ₛ[3] + r⃗ₜ[3])) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + ((((-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-r⃗ᵤ[3] + r⃗ᵥ[3])) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)))^2
p̂₁ = Symbolics.Num[r⃗ₛ[1] + ((((((2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(⃗rₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) + 2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-r⃗ₛ[1] + r⃗ₜ[1])) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)), r⃗ₛ[2] + ((((((2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - ⃗rₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) + 2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-r⃗ₛ[2] + r⃗ₜ[2])) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)), r⃗ₛ[3] + ((((((2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)) + 2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) + 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-r⃗ₛ[3] + r⃗ₜ[3])) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2))]
p̂̂₂ = Symbolics.Num[r⃗ᵤ[1] + ((-r⃗ᵤ[1] + r⃗ᵥ[1])*(((-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-⃗rᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)), r⃗ᵤ[2] + ((((-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-r⃗ᵤ[2] + r⃗ᵥ[2])) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2)), r⃗ᵤ[3] + ((((-2(r⃗ᵤ[1] - r⃗ₛ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(r⃗ᵤ[2] - r⃗ₛ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(r⃗ᵤ[3] - r⃗ₛ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))*(-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) + 2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ᵤ[1] - r⃗ₛ[1]) + 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ᵤ[2] - r⃗ₛ[2]) + 2(r⃗ᵤ[3] - r⃗ₛ[3])*(-r⃗ᵤ[3] + r⃗ᵥ[3]))*(-r⃗ᵤ[3] + r⃗ᵥ[3])) / ((-((-2(-r⃗ᵤ[1] + r⃗ᵥ[1])*(r⃗ₛ[1] - r⃗ₜ[1]) - 2(-r⃗ᵤ[2] + r⃗ᵥ[2])*(r⃗ₛ[2] - r⃗ₜ[2]) - 2(-r⃗ᵤ[3] + r⃗ᵥ[3])*(r⃗ₛ[3] - r⃗ₜ[3]))^2)) / (-2((r⃗ₛ[1] - r⃗ₜ[1])^2) - 2((r⃗ₛ[2] - r⃗ₜ[2])^2) - 2((r⃗ₛ[3] - r⃗ₜ[3])^2)) - 2((-r⃗ᵤ[1] + r⃗ᵥ[1])^2) - 2((-r⃗ᵤ[2] + r⃗ᵥ[2])^2) - 2((-r⃗ᵤ[3] + r⃗ᵥ[3])^2))]
Main.MinimumDistanceBetweenTwoLines
```

# Also find the end points of the shortest line segment

The formulation above can be reused, the optimal lambdas $\hat{\Lambda} = (\lambda_1, \lambda_2)$, correspond to the closest points through the equation of the lines.

$$\left(\hat{\lambda}_1, \hat{\lambda}_2\right) = \underset{\lambda_1, \lambda_2}{\texttt{argmin}} \left\{ D(\lambda_1, \lambda_2) \right\}$$

A specific example's numerical output was given in the previous section.

# Find the volume of the tetrahedron

For any 2d shape extended to a point, the volume is $\frac{1}{3} \text{Base-Area} \times \text{height}$. Why? Take a slice, the planar dimensions reduce linearly to zero, so the area reduces quadratically. Integrate the function for the area as a function of the height from base to the tip.

This way, we can see that a tetrahedron is geometric one-sixth of a parallelapiped, therefore, volume $V(\vec{\mathbf{r}}_a, \vec{\mathbf{r}}_b, \vec{\mathbf{r}}_c, \vec{\mathbf{r}}_d)$ is (assuming the vectors are in three dimensions for cross product to work):

$$V(\vec{\mathbf{r}}_a, \vec{\mathbf{r}}_b, \vec{\mathbf{r}}_c, \vec{\mathbf{r}}_d) = \frac{1}{6} \left|(\vec{\mathbf{r}}_c - \vec{\mathbf{r}}_b) \times (\vec{\mathbf{r}}_c - \vec{\mathbf{r}}_a) \cdot (\vec{\mathbf{r}}_d - \vec{\mathbf{r}}_d)\right|$$

# Write formula or code for the tensions in the rods

$$\vec{\mathbf{F}} + T_{ad}\hat{r}_{ad} + T_{bd}\hat{r}_{bd} + T_{cd}\hat{r}_{cd} = \vec{\mathbf{0}}$$

There was an in-class mention about Cramer's rule to solve this. There are three equations and three unknowns, but one could use Gaussian-elimination.

The other way to solve this using vector algebra is to dot the whole equation with a vector which is perpendicular to two of the three tension vectors. For example by $\hat{r}_{ad}$.
