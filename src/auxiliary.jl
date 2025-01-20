# © Ylann Rouzaire 2025 `rouzaire.ylann@gmail.com`

logspace(x1, x2, n; digits=1) = unique!(round.([10.0^y for y in range(log10(x1), log10(x2), length=n)], digits=digits))

each = eachindex # handy alias, but not optimized. Do not use inside intensive loops

function remneg(x)
    x[x.<0] .= NaN
    return x
end
remove_negative = remneg # alias

function pz(z) # to print the runtime in a readable format
    ss = "$(round(Int,z)) seconds"
    mm = "$(round(z/60,digits=2)) minutes"
    hh = "$(round(z/3600,digits=2)) hours"
    # dd = "$(round(z/86400,digits=2)) days"
    println("Runtime : \n  $ss \n  $mm \n  $hh")
    return z
end
prinz = pz # alias

function arclength(theta1, theta2)
    #= This function returns the signed arclength (in radians) from theta1 to theta2 on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    WARNING : Note that the inputs thetas need to lie within [0,2π] 
    Examples : 
    arclength(-1,-2) = -1
    arclength(1,2) = 1
    arclength(1,2pi-1) = -2
    =#
    pi_32 = Float32(pi) # otherwize throws an error when comparing dtheta_abs < pi
    dtheta = theta2 - theta1
    dtheta_abs = abs(theta2 - theta1)

    shortest_unsigned_arclength = min(2 * pi_32 - dtheta_abs, dtheta_abs) # this is why one need the inputs to lie within [0,2π]
    if dtheta_abs ≤ pi_32
        signe = sign(dtheta)
    else
        signe = -sign(dtheta)
    end
    return signe * shortest_unsigned_arclength
end