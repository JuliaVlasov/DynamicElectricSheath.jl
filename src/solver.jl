export compute_charge!


"""
$(SIGNATURES)

compute charge density
"""
function compute_charge!(rho, f, dv)
    rho .= vec(sum(f, dims = 2)) * dv
end

export compute_current

"""
$(SIGNATURES)

compute current at boundaries : rectangle integration of 

```math
\\int_{v} v * (fi(t,\\pm 1,v) - fe(t,\\pm 1,v)) dv
```
"""
function compute_current(fi, fe, vv, dv)

    @views J_l = dv * sum(vv .* (fi[begin, :] - fe[begin, :]))
    @views J_r = dv * sum(vv .* (fi[end, :] - fe[end, :]))
    return J_l, J_r

end

export compute_e!

"""
$(SIGNATURES)

Update electric field : use formula 

```math
E(x) = 0.5*(E(1) + E(-1)) + \\nu/2 * ( \\int_{-1}^{x} rho(y)dy - \\int_{x}^{1} \\rho(y)dy )
```
with ``E(t^{n+1},1) + E(t^{n+1},-1)`` approximated by ``E(t^n,1) + E(t^n,-1) + dt/\\lambda^2 * (J_l + J_r)``
(Ampère boundary conditions)
"""
function compute_e!(EE, rho, lambda, J_l, J_r, dx, dt)
    @views EE .= 0.5 * (EE[begin] + EE[end] + dt / (lambda^2) * (J_l + J_r))
    # add forward  trapeze integral
    @views EE[begin+1:end] .+=
        0.5 / (lambda^2) * 0.5 * dx .* cumsum(rho[1:end-1] .+ rho[2:end])
    # add backward trapeze integral
    @views EE[begin:end-1] .-=
        0.5 / (lambda^2) * 0.5 * dx .*
        cumsum(rho[end-1:-1:1] .+ rho[end:-1:2])[end:-1:begin]
    EE_plus = max.(EE, 0.0)
    EE_minus = min.(EE, 0.0)
    return EE_minus, EE_plus
end

export advection!

"""
$(SIGNATURES)

```math
\\frac{d f_i}{dt} + v \\frac{d f_i}{dx} + E \\frac{d f_i}{dv} = \\nu * f_e
```

```math
\\frac{d f_e}{dt} + v \\frac{d f_e}{dx} + \\frac{1}{\\mu] E \\frac{d f_e}{dv} = 0
```
"""
function advection!(fi, fe, vv_plus, vv_minus, EE_plus, EE_minus, nu, mu, dx, dv, dt)

    Nx = size(fi, 1) - 1
    Nv = size(fi, 2) - 1

    @views fi[2:Nx, 2:Nv] .+= (
        -(dt / dx) .* (
            vv_plus[2:Nv]' .* (fi[2:Nx, 2:Nv] .- fi[1:Nx-1, 2:Nv]) .+
            vv_minus[2:Nv]' .* (fi[3:Nx+1, 2:Nv] .- fi[2:Nx, 2:Nv])
        ) .-
        (dt / dv) .* (
            EE_plus[2:Nx] .* (fi[2:Nx, 2:Nv] .- fi[2:Nx, 1:Nv-1]) .+
            EE_minus[2:Nx] .* (fi[2:Nx, 3:Nv+1] .- fi[2:Nx, 2:Nv])
        ) .+ dt * nu .* fe[2:Nx, 2:Nv]
    )
    @views fe[2:Nx, 2:Nv] .+= (
        -(dt / dx) .* (
            vv_plus[2:Nv]' .* (fe[2:Nx, 2:Nv] .- fe[1:Nx-1, 2:Nv]) .+
            vv_minus[2:Nv]' .* (fe[3:Nx+1, 2:Nv] .- fe[2:Nx, 2:Nv])
        ) .+
        (dt / (dv * mu)) .* (
            EE_minus[2:Nx] .* (fe[2:Nx, 2:Nv] .- fe[2:Nx, 1:Nv-1]) .+
            EE_plus[2:Nx] .* (fe[2:Nx, 3:Nv+1] .- fe[2:Nx, 2:Nv])
        )
    )

end
