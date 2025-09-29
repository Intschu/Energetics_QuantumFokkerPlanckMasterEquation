# measurement driven heat engine

using Plots

setprecision(BigFloat, 1000) # set numerical precision

# function outputting a solved recursive scheme for a given bandwidth (gam), measurement rate (lam), and reservoir coupling strength (k0)
function get_coefficients(gam, lam, k0)

    kp = k0*nb # excitation rate (kappa_+)
    km = k0*(1.0+nb) # decay rate (kappa_-)
    sig = gam/(8.0*lam)  # parameter sigma of Hermite polynomials

    # set up a system of 3*cutoff + 3 linear equations (including probability normalization)
    # basis { p_0, a^z_0, a^x_1, p_2, a^z_2, a^x_3, ..., p_(2*cutoff), a^z_(2*cutoff), a^x_(2*cutoff+1) }
    # entries p_2j, a^z_2j, a^x_(2j+1) -> indices 1 + 3*j, 1 + 3*j +1, 1 + 3*j + 2 

    m = fill(BigFloat("0.0"), (3*cutoff+3, 3*cutoff+3));
    b = fill(BigFloat("0.0"), 3*cutoff+3);
    b[1] = BigFloat("1.0") # normalization p_0 = 1

    j = 0

    m[1+3*j, 1+3*j] = BigFloat("1.0") # normalization p_0 = 1

    m[1+3*j+1, 1+3*j+1] = g
    m[1+3*j+1, 1+3*(j+1)+1] = g*2.0*(j+1)*sig
    m[1+3*j+1, 1+3*j] = gam/(2.0*sig)
    m[1+3*j+1, 1+3*j+2] = -( gam*(2.0*j+1)/2.0 + (kp+km)/4.0)


    m[1+3*j+2, 1+3*j+2] = -g*(2*j+1)*sig
    m[1+3*j+2, 1+3*j+1] = -(lam + (kp+km)/2.0 + gam*2*j/2.0)
    m[1+3*j+2, 1+3*j] = (kp-km)/2

    for j = 1:cutoff-1
        m[1+3*j, 1+3*j] = 2.0*j*BigFloat("1.0")
        m[1+3*j, 1+3*(j-1)+2] = -1/sig

        m[1+3*j+1, 1+3*j +1] = g
        m[1+3*j+1, 1+3*(j+1) +1] = g*2.0*(j+1)*sig
        m[1+3*j+1, 1+3*j] = gam/(2.0*sig)
        m[1+3*j+1, 1+3*j+2] = -( gam*(2.0*j+1)/2.0 + (kp+km)/4.0)

        m[1+3*j+2, 1+3*j+2] = -g*(2*j+1)*sig
        m[1+3*j+2, 1+3*j+1] = -(lam + (kp+km)/2.0 + gam*2*j/2.0)
        m[1+3*j+2, 1+3*j] = (kp-km)/2.0
        m[1+3*j+2, 1+3*(j-1)+2] = -g
    end

    j = cutoff

    m[1+3*j, 1+3*j] = 2*j*BigFloat("1.0")
    m[1+3*j, 1+3*(j-1)+2] = -1/sig

    m[1+3*j+1, 1+3*j +1] = g
    #m[1+3*j+1, 1+3*(j+1) +1] = g*2.0*(j+1)*sig
    m[1+3*j+1, 1+3*j] = gam/(2.0*sig)
    m[1+3*j+1, 1+3*j+2] = -( gam*(2.0*j+1)/2.0 + (kp+km)/4.0)

    m[1+3*j+2, 1+3*j+2] = -g*(2*j+1)*sig
    m[1+3*j+2, 1+3*j+1] = -(lam + (kp+km)/2.0 + gam*2*j/2)
    m[1+3*j+2, 1+3*j] = (kp-km)/2
    m[1+3*j+2, 1+3*(j-1)+2] = -g

    return  m\b # return a list of coefficients
end


function steady_state(gam, lam, k0) # compute steady-state density matrix and energy quantities

    kp = k0*nb # excitation rate (kappa_+)
    km = k0*(1.0+nb) # decay rate (kappa_-)
    sig = gam/(8.0*lam) # parameter sigma of Hermite polynomials

    d_vec = range(-d_max, d_max, points) # vector of values of d = D
    p_vec = [] # vector of p(D)
    z_vec = [] # vector of a_z(D)
    x_vec = [] # vector of a_x(D)
    
    coef = get_coefficients(gam, lam, k0) # get a list of coefficients

    for d in d_vec # loop over D

        u = d/(sqrt(2.0*sig))
        hermite = [1, 2*u] # standard Hermite polynomials of order 0 and 1
        # recursively construct standard Hermite polynomials of orders from 2 to 2*cutoff+1
        for n = 2:(2*cutoff+2)
            push!(hermite, 2.0*u*hermite[1+n-1]-2.0*(n-1)*hermite[1+n-2])
        end
        
        p = 0.0 # p(D)
        x = 0.0 # a_x(D)
        z = 0.0 # a_z(D)

        # loop over coefficients and weight by physicist's Hermite polynomials
        for l = 0:cutoff
            p += coef[1+3*l]*hermite[1+2*l]*((sig/2)^(l))
            z += coef[1+3*l+1]*hermite[1+2*l]*((sig/2)^(l))
            x += coef[1+3*l+2]*hermite[1+2*l+1]*((sig/2)^(l+1/2))
        end
        
        # weight by Gaussian distributions 
        push!(p_vec, exp(-u^2)*p/(sqrt(2.0*pi*sig))) # append p(D)
        push!(z_vec, exp(-u^2)*z/(sqrt(2.0*pi*sig))) # append a_z(D)
        push!(x_vec, exp(-u^2)*x/(sqrt(2.0*pi*sig))) # append a_x(D)

    end

    # return a list of D values and three lists corresponding to p(D), a_x(D), a_z(D)
    # together with the list [power, measurment energy rate, heat current] (in units of energy splitting omega)
    return d_vec, p_vec, x_vec, z_vec, [-g*sig*coef[3], -lam*coef[2], (kp-km)/2.0 - coef[2]*(kp + km)/2.0]

end



points = 401 # number of plotted points in a detector of D's
d_max = 2.5 # boundary value of detector D for plotting
cutoff = 100 # cutoff in the recursive scheme

# set parameters of the system
g = BigFloat("1.0") # coherent driving
nb = BigFloat("0.58") # Bose-Einstein distribution
gam = BigFloat("1.0") # bandwidth
lam = BigFloat("1.0") # measurement rate
k0 = BigFloat("0.33") # reservoir coupling strength

global d_vec , p_vec, x_vec, z_vec, energy_terms = steady_state(gam, lam, k0) # get properties of the steady state

plot_prob = plot(d_vec, [p_vec, x_vec]) # plot p(D) and a_x(D)

lam_vec = range(0.1, 1.0, 10) # vector of measurement rate lambda
power = [] # power
meas_energy = [] # measurement energy
heat = [] # heat current

for l in lam_vec # loop over measurement rate lambda
    global d_vec , p_vec, x_vec, z_vec, energy_terms = steady_state(gam, l, k0) # get energy terms
    push!(power, energy_terms[1]) # append power
    push!(meas_energy, energy_terms[2]) # append measurement energy rate
    push!(heat, energy_terms[3]) # append heat current
end

# plot entracted power, measurement energy rate, and heat current (in units of energy splitting omega) as a function of measurement rate lambda
plot_energy = plot(lam_vec, [-power, meas_energy, heat])

plot(plot_prob, plot_energy, layout = (2, 1)) # display plots
