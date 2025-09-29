# two level system

using Plots, SpecialFunctions

setprecision(BigFloat, 1000) # set numerical precision

function factorial2(x) # double factorial x!!

    y = 1.0
    while x >= 2.0
        y *= x
        x -= 2.0
    end

    return y     

end 

function fac2(x) # double factorial x!!

    return sqrt(2.0/pi)*(2.0^(x/2.0))*gamma(x/2.0 + 1.0)  

end 

# function outputting a solved recursive scheme for a given bandwidth (gam), measurement rate (lam), and reservoir coupling strength (k0)
function get_coefficients(gam, lam, k0)

    kp = k0*nb # excitation rate (kappa_+)
    km = k0*(1.0+nb) # decay rate (kappa_-)
    sig = gam/(8.0*lam) # parameter sigma of Hermite polynomials

    # set up a system of 2*cutoff + 2 linear equations (including probability normalization)
    # basis { q^+_0, q^-_1, q^+_2, q^-_3, ..., q^+_(2*cutoff), q^-_(2*cutoff+1) }
    # entries q^+_2j, q^-_(2j+1) -> indices 1 + 2*j, 1 + 2*j +1 

    m=fill(BigFloat("0.0"),(2*cutoff + 2, 2*cutoff+2));
    b=fill(BigFloat("0.0"), 2*cutoff + 2);
    b[1] = BigFloat("1.0") # normalization q^+_0 = 1

    j = 0

    #q^+ 2*0
    m[1+2*j, 1] = BigFloat("1.0") # normalization q^+_0 = 1

    #q^- 2*0 + 1
    m[1+2*j+1, 1+2*j+1] = -gam*(2.0*j+1) - (kp+km)
    m[1+2*j+1, 1+2*j] = -(gam/sig)
    for k = 0:cutoff
        m[1+2*j+1, 1+2*k] += 2*(kp-km)*((-1)^(j+k))*((sig)^(k-j-1/2))*factorial2(2*j+1)*factorial2(2*k-1)/(gamma(1+2*j+1)*(2*j-2*k+1)*sqrt(2.0*pi))
    end

    for j = 1:cutoff
        #q^+ 2*j
        m[1+2*j, 1+2*j] = -gam*(2.0*j)
        m[1+2*j, 1+2*(j-1)+1] = -gam/sig

        #q^- 2*j + 1
        m[1+2*j+1, 1+2*j+1] = -gam*(2*j+1) - (kp+km)
        m[1+2*j+1, 1+2*j] = -(gam/sig)
        for k = 0:cutoff
            m[1+2*j+1, 1+2*k] += 2*(kp-km)*((-1)^(j+k))*((sig)^(k-j-1/2))*factorial2(2*j+1)*factorial2(2*k-1)/(gamma(1+2*j+1)*(2*j-2*k+1)*sqrt(2.0*pi))
        end
        j += 1
    end

    return m\b # return a list of coefficients
end



function steady_state(gam, lam, k0) # compute steady-state density matrix and energy quantities
    
    coef = get_coefficients(gam, lam, k0) # get a list of coefficients

    sig = gam/(8.0*lam) # parameter sigma of Hermite polynomials
    d_vec = range(-d_max, d_max, points) # vector of values of d = D
    q_plus_vec = []  # vector of q_+(D)
    q_minus_vec = [] # vector of q_-(D)

    error_prob = 0.0 # detector error probability

    for d in d_vec # loop over D

        u = d/(sqrt(2.0*sig))
        hermite = [1, 2.0*u] # standard Hermite polynomials of order 0 and 1
        # recursively construct standard Hermite polynomials of orders from 2 to 2*cutoff+1
        for n = 2:(2*cutoff+2)
            push!(hermite, 2.0*u*hermite[1+n-1]-2.0*(n-1)*hermite[1+n-2])
        end

        q_plus = 0.0 # q_+(D)
        q_minus = 0.0 # q_-(D)
        
        for l = 0:cutoff # loop over coefficients and weight by physicist's Hermite polynomials
            q_plus += coef[1+2*l]*hermite[1+2*l]*((sig/2)^(l))
            q_minus += coef[1+2*l+1]*hermite[1+2*l+1]*((sig/2)^(l+1/2))
        end

        # weight by Gaussian distributions
        push!(q_plus_vec, exp(-u^2)*q_plus/(sqrt(2.0*pi*sig))) # append q_+(D)
        push!(q_minus_vec, exp(-u^2)*q_minus/(sqrt(2.0*pi*sig))) # append q_-(D)
    
    end
 
    
    # return a list of D values and two lists corresponding to q_+(D), q_-(D)
    return d_vec, q_plus_vec, q_minus_vec

end


nb = BigFloat("0.58") # Bose-Einstein distribution
cutoff = 100 # cutoff in the recursive scheme
points = 801 # number of plotted points in a detector of D's
d_max = 5.0 # boundary value of detector D for plotting
step = 2.0*d_max/(points-1)

# set parameters of the system
gam = BigFloat("1.0") # bandwidth
lam = BigFloat("0.5") # measurement rate
k0 = BigFloat("0.1") # reservoir coupling strength

global d_vec , q_plus_vec, q_minus_vec = steady_state(gam, lam, k0) # get properties of the steady state

plot_prob = plot(d_vec, q_plus_vec) # plot q_+(D)


plot(plot_prob) # display plots