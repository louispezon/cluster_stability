using Plots
using Statistics
using LaTeXStrings
using Printf
using ProgressBars
using NLsolve

###########
function asynch(μ,J)
    """compute the stationnary activity, in units of 1/τ_h"""
    function I(A)
        return (1-exp(-1/A)) * (μ + J*A) - 1
    end
    I!(arr)=[I(arr[1])]
    #return fixedpoint(I!,[0.1]).zero[1]
    return nlsolve(I!, [1.]).zero[1]
end

function sample_asynch(μ,N)
    """sample the stationnary voltage distribution"""
    C = 1/log(μ/(μ-1))
    U = 1/N*range(1,N)#rand((N))
    #U = sort(rand((N)))
    Finv(u) = μ*(1-exp(-u/C))
    v0 = [Finv(u) for u in U]
end

# neuron parameters
mu_0 = 1.2
mu_1 = 1.2
J = 0.4

#τ_h = 10 # ms ##### !!!! all times are expressed in units of τ_h !!!!
# kernel params:
Δ = 0.5 # τ_h
τ = 0.2 # τ_h

T = 50  # τ_h
T_step = T/2

μ(t) = mu_0*(t<T_step) + mu_1*(t>=T_step) # + mu_0*(t>=2*T_step)

# intensity function
beta = 1e5 # inverse temperature 
f(v) = beta*(v>=1)# exp(beta*(v-1))

# simulation parameters
N = 100 # number of neurons
dt = 5e-4 # τ_h

delay = Int(Δ/dt)
times = 0:dt:T


### Initialisation 

v = zeros(N,length(times))
spikes = zeros(N,length(times))
r = zeros(length(times))

εA_0 = zeros(length(times))
εB_0 = zeros(length(times))
εA = zeros(length(times))
#εB = zeros(length(times))


# initialization
## synchronous state: 
#v[:,1]=.02*rand((N))
#r[1] = 1/dt

## Asynchronous state: 
v[:,1]=sample_asynch(mu_0,N)
r[1] = asynch(mu_0,J) # in units of 1/τ_h

εA_0[1] = asynch(mu_0,J)

function run(steps,μ)

    for t in ProgressBar(steps)
    
        if t-delay<1 
            εA[t] = εA_0[t]
        else
            εA[t] = εA_0[t-delay]
        end

        for n in 1:N
            if 1-exp(-f(v[n,t])*dt) > rand() #f(v[n,t]) >0 
                v[n,t+1] = 0
                spikes[n,t+1] = 1
            else
                v[n,t+1] = v[n,t] + dt*(μ[t] - v[n,t] + J*εA[t])
            end
        end

        εB_0[t+1] = εB_0[t]+ dt/τ * (-εB_0[t] + r[t])
        εA_0[t+1] = εA_0[t]+ dt/τ * (-εA_0[t]+εB_0[t])

        r[t+1] = mean(spikes[:,t+1])/dt
    end
end

nstep = length(times)-1
steps=1:nstep
mu_t = [μ(times[t]) for t in steps]
run(1:nstep,mu_t)


#using JLD2
#name = "def_name"
#@save name*".jld2"


############################################################################
################### EXTRACT DATA  ##########################################
############################################################################

function compute_A(r,τ_A = 0.05)
    """Compute the activity (in units of 1/τ_h), with a low-pass filtering of time constant τ_A (in units of τ_h)."""
    A = zeros(length(r))
    A[1] = r[1]*dt/τ_A
    for t in 1:length(r)-1
        A[t+1] = A[t] + dt/τ_A*(-A[t] + r[t+1])
    end
    return A
end
A = compute_A(r)


########## Find matching spikes

function compute_firing_delays(;nref=1, n_del=25)
    times_nref = [0.]
    indices = 1:Int(floor(N/n_del)):N
    all_delays = Array{Array}(undef,n_del+1) # firing delays wrt nref : 
    for i in 1:n_del+1
        all_delays[i] = [0.]
    end
    for t in 1:nstep
        if spikes[nref,t]==1
            append!(all_delays[n_del+1],times[t]-times_nref[end])
            append!(times_nref,times[t])
        end
        for n in 1:n_del
            if spikes[indices[n],t]==1
                append!(all_delays[n],times[t]-times_nref[end])
            end
        end
    end

    ## delete irrelevant initial values
    for a in all_delays
        deleteat!(a,1)
    end
    deleteat!(times_nref,1)
    
    return all_delays,times_nref
end

####################################### Raster plot
function get_raster_data(n_pl;n_raster=25)
    indices = 1:Int(N/n_raster):N

    raster_inds = []
    raster_times= []
    for ind in 1:length(indices)
        for t in n_pl-shortstep:n_pl
            if spikes[indices[ind],t]==1
                append!(raster_inds,ind)
                append!(raster_times,times[t])
            end
        end
    end
    return raster_times,raster_inds
end



all_delays, times_nref = compute_firing_delays(nref=1,n_del = 25)
############################################################################
################### PLOT ###################################################
############################################################################
#gr() ; 
plot_font = "Serif"
pyplot()

default(fontfamily=plot_font,framestyle=:box, label=nothing,legendfontsize=12, 
guidefontsize=14, grid=false, tickfontsize=10, titlefontsize=17, reuse=false)#, usetex=true)

title = @sprintf("\$ \\mu = %g \\rightarrow %g ,\\  J =  %g \$", mu_0,mu_1,J)
#title=LaTeXString(title)
shortstep = Int((5)/dt)

T_pl = T/2
n_pl = Int(T_pl/dt)
p2 = plot(times[n_pl-shortstep:n_pl],A[n_pl-shortstep:n_pl], lw=1.5, ylabel=L"A \ [\tau_h^{-1}]", color=:red)
raster_times,raster_inds = get_raster_data(n_pl)
p3 = scatter(raster_times,raster_inds, xlabel=L"t", ylabel="Neuron index", color=:red)

T_pl = T
n_pl = Int(T_pl/dt)
p4 = plot(times[n_pl-shortstep:n_pl],A[n_pl-shortstep:n_pl], lw=1.5, color=:magenta)
raster_times,raster_inds = get_raster_data(n_pl)
p5 = scatter(raster_times,raster_inds, xlabel=L"t \ [\tau_h]", color=:magenta)

#display(p1)
display(plot(p2,p4,p3,p5,layout=(2,2),link=:x))



# #Early activity:
# (times[1:shortstep],A[1:shortstep], lw=1.5, ylabel=L"$A$ [Hz]", color=:red)
####################################### Delays

function plot_del(pl)
    n_del = length(all_delays)-1
    for n in 2:n_del
        l = min(length(times_nref),length(all_delays[n]))
        plot!(pl,times_nref[1:l],all_delays[n][1:l], markershape=:circle, markersize = 1.5, lw=0, ls=:dot, markerstrokewidth=.15)
    end
    plot!(pl,times_nref,all_delays[n_del+1],lw=2,lc=:black, label=L"t_1^{f+1} - t_1^f", legend=:topleft)
end

p1=plot([])
plot_del(p1)
title!(p1,title )
xlabel!(L"t_1^f")
ylabel!(L"t_n^f - t_1^f")
xlims!((0,T_pl))

plot!([T_step,T_step],[0,1.2],
    #label=L"t=t_\mathrm{step}", 
    lc=:red, 
    ls=:dash)


display(p1)

