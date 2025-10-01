using Random
using LinearAlgebra
using QuantumOptics
using Plots
using Distributions

# sample initian a_0
function aIni(p::Float64)
    if rand(Float64) < p
        return 1.0
    else
        return -1.0
    end
end

# sample initial d_0
function dIni(avg::Float64, sig::Float64)
    return rand(Normal(avg, sig))
end

# update d
function dNext(a::Float64, d::Float64, gam::Float64, lam::Float64)
    #return d + gam*dt*(rand(Normal(a, sig)) - d)
    return d + gam*dt*(a-d) + sqrt(dt)*gam/(2.0*sqrt(lam))*rand(Normal(0, 1))
end

# update a
function aNext(a, dNew, kp, km)
    pow = 0 # extracted energy quanta
    if a == 1.0
        if dNew >= 0.0
            if rand(Float64) < dt*kp
                pow += 1
                aNew = -1.0
            else
                aNew = a
            end
        else
            if rand(Float64) < dt*km
                pow += -1
                aNew = -1.0
            else
                aNew = a
            end
        end
    else
        if dNew >= 0.0
            if rand(Float64) < dt*km
                pow += -1
                aNew = 1.0
            else
                aNew = a
            end
        else
            if rand(Float64) < dt*kp
                pow += 1
                aNew = 1.0
            else
                aNew = a
            end
        end
    end
    return [aNew, pow]
end


# sample trajectory
function trajectory(gam::Float64, lam::Float64, kap::Float64)
    a = aIni(0.5) # sample initial state
    d = dIni(a, sqrt(gam/(8.0*lam))) # sample initial detector
    kp = kap*be # raising rate
    km = kap*(1.0+be) # lowering rate
    a_vec = [a] # trajectory of system
    d_vec = [d] # trajectory of detector
    for j = 1:N
        d = dNext(a, d, gam, lam) 
        a = aNext(a, d, kp, km)[1]
        push!(a_vec, a)
        push!(d_vec, d)
    end

    outfile_traj_sys = "C:\\Users\\trajsys.txt"
    outfile_traj_det = "C:\\Users\\trajdet.txt"

    open(outfile_traj_sys, "w") do f
        for j = 1:N
          println(f, string(dt*j*gam)*" "*string(Float64(a_vec[j])))
        end
    end

    open(outfile_traj_det, "w") do f
        for j = 1:N
          println(f, string(dt*j*gam)*" "*string(Float64(d_vec[j])))
        end
    end
    println("trajectory done")
end


#feedback run
function feedbackRun(gam::Float64, lam::Float64, kap::Float64)
    kp = kap*be # raising rate
    km = kap*(1.0+be) # lowering rate

    ent1 = 0.0 # measurement entropy version 1
    ent2 = 0.0 # measurement entropy version 12
    pow = 0.0 # extracted energy quanta
    ft = 0.0 # \exp{-sigma-sigma_m}
    ftOld = 0.0 # \exp{-sigma}

    a = aIni(0.5) # initial system
    d = dIni(a, sqrt(gam/(8.0*lam))) # initial detector
    tempIni = exp((-4.0*lam*(d-a)^2)/gam) # initial boundry term to ft
    ent1 += (a - d)^2
    dNew = dNext(a, d, gam, lam)
    ent2 += (a + a - dNew - d)*(dNew - d)
    d = dNew

    for j = 1:N        
        temp = aNext(a, d, kp, km)
        pow += temp[2]
        a = temp[1]
        ent1 += (a - d)^2

        dNew = dNext(a, d, gam, lam)
        ent2 += (a + a - dNew - d)*(dNew - d)
        d = dNew
    end
    
    tempFin = exp((-4.0*lam*(d-a)^2)/gam) # final boundry term to ft
    
    ft += tempFin*exp(pow-4.0*lam*ent2/gam )/tempIni
    ftOld += exp(pow)

    return ft, ftOld, pow, 8.0*lam*dt*ent1, 4.0*lam*ent2/gam 
end


#=
Random.seed!(123)
N = 2000
dt = 0.05
#runs = 200
be = 0.582
gam = 1.0
kap = 0.01
lam = 0.1
#trajectory(gam, lam, kap)
=#

# ensemble averaging
function simfeedback(gam::Float64, lam::Float64, kap::Float64)
    ft = 0.0
    ftOld = 0.0
    pow = 0.0
    ent1 = 0.0
    ent2 = 0.0
    
    temp1, temp2, temp3, temp4, temp5 = 0.0, 0.0, 0.0, 0.0, 0.0
    for r = 1:runs
        temp1, temp2, temp3, temp4, temp5 = feedbackRun(gam, lam, kap)
        ft += temp1
        ftOld += temp2
        pow += temp3
        ent1 += temp4
        ent2 += temp5
    end
    return [ft/runs, ftOld/runs, pow/(runs*N*dt*kap), ent1/(runs*N*dt*kap) - gam/kap, ent2/(runs*dt*N*kap)]
end


function plotSim(gam, lamVec, kap)

    outfile_ent = "C:\\Users\\entSimfdl.txt"
    

    open(outfile_ent, "w") do f
        for l in lamVec
          println(l)
          temp = simfeedback(gam, l, kap)
          println(f, string(l/gam)*" "*string(Float64(temp[3]))*" "*string(Float64(temp[4]))*" "*string(Float64(temp[5])))
        end
    end


    println("ent is done")

end
#Random.seed!(123)

N = 1000
dt = 0.01
runs = 1000000
be = 0.582
gam = 1.0
kap = 0.1

lam = 0.3

lamVec = range(0.001, 0.05, 10)

#plotSim(gam, lamVec, kap)

println("Hello")
println(simfeedback(gam, lam, kap))




