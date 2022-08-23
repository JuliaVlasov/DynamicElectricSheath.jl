```@meta
CurrentModule = DynamicElectricSheath
```

# Example


Set physical and numerical parameters

```@example test
using DynamicElectricSheath
using Plots

physics = Physics()

ν = physics.ν
λ = physics.λ
T = physics.T
μ = physics.μ
xmin, xmax = physics.xmin, physics.xmax
vmin, vmax = physics.vmin, physics.vmax

numerics = Discretization(physics)

Nt = numerics.Nt
dt = numerics.dt
CFL_x = numerics.CFL_x
CFL_v = numerics.CFL_v
Nx = numerics.Nx
dx = numerics.dx
Nv = numerics.Nv
dv = numerics.dv

println("Parameters : ")
println("ν = $ν, λ = $λ, μ = $μ")
println("T = $T, Nt = $Nt, dt = $dt (CFL_x = $CFL_x, CFL_v = $CFL_v)")
println("[xmin,xmax] = [$xmin,$xmax], Nx = $Nx, dx = $dx")
println("[vmin,vmax] = [$vmin,$vmax], Nv = $Nv, dv = $dv")
```

Set mesh grid 

```@example test
xx = LinRange(xmin, xmax, Nx + 1)
vv = LinRange(vmin, vmax, Nv + 1)
vv_plus = vv .* (vv .> 0.0)
vv_minus = vv .* (vv .< 0.0) 
```

Initialize fields

```@example test
EE = E0.(xx)       
fi = fi_0.(xx, vv')  
fe = fe_0.(xx, vv') 
```

```@example test
p = plot(layout=(1,2), size=(800,300))
contourf!(p[1], xx, vv, fi', label="ions")
contourf!(p[2], xx, vv, fe', label="electrons")
```


```@example test
ρ = zeros(Nx + 1)       # charge density
ρi = zeros(Nx + 1)        # ion charge density
ρe = zeros(Nx + 1)        # electron charge density
compute_charge!(ρi, fi, dv)
compute_charge!(ρe, fe, dv)
compute_charge!(ρ, fi .- fe, dv)
plot(xx, ρ)
```

```julia
# Boundary conditions (all 0, never updated afterwards)
fi[begin, :] .= 0.0
fi[end, :] .= 0.0 # speed distribution is almost 0 
fi[:, begin] .= 0.0
fi[:, end] .= 0.0 # non-emmiting wall
fe[begin, :] .= 0.0
fe[end, :] .= 0.0 # speed distribution is almost 0 
fe[:, begin] .= 0.0
fe[:, end] .= 0.0 # non-emmiting wall


for n = 1:Nt # loop over time

    compute_charge!(ρi, fi, dv)
    compute_charge!(ρe, fe, dv)
    compute_charge!(ρ, fi .- fe, dv)

    J_l, J_r = compute_current(fi, fe, vv, dv)

    EE_minus, EE_plus = compute_e!(EE, ρ, λ, J_l, J_r, dx, dt)

    advection!(fi, fe, vv_plus, vv_minus, EE_plus, EE_minus, ν, μ, dx, dv, dt)

end

```