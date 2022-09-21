using Plots
using Statistics
using LaTeXStrings
using Printf

#using JLD2
using FileIO

#name = "Cluster_trans.jld2"
r,N,spikes,times, dt=load("Cluster_trans.jld2","r","N","spikes","times", "dt","T")

unit=1/100 # time is in seconds
times*=unit

mu_0 = 1.2
mu_1 = 1.1
J = 0.4
nstep = length(times)-1
steps=1:nstep
T = times[end]
T_step = T/2
dt = T/nstep

T_short = 5*unit
τ_A = 2e-2*unit
############################################################################
################### EXTRACT DATA  ##########################################
############################################################################


function compute_A(r,τ_A; τmult=1)
    A = zeros(length(r))
    A[1] = r[1]*dt/τ_A
    for t in 1:length(r)-1
        A[t+1] = A[t] + dt/τ_A*(-A[t] + r[t+1] * τmult)
    end
    return A
end
A = compute_A(r,τ_A,τmult=1)


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

all_delays, times_nref = compute_firing_delays(nref=1,n_del = 25)
all_delays*= 1/unit*10

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



############################################################################
################### PLOT ###################################################
############################################################################
#gr() ; 
plot_font = "Serif"
pyplot()

default(fontfamily=plot_font,framestyle=:box, label=nothing,legendfontsize=12, 
guidefontsize=14, grid=false, tickfontsize=10, titlefontsize=17, reuse=false)#, usetex=true)

#title = @sprintf("\$ \\mu = %g \\rightarrow %g ,\\  J =  %g \$", mu_0,mu_1,J)
title = @sprintf("\$ \\mu = %g \\rightarrow %g \$", mu_0,mu_1)
#title=LaTeXString(title)
shortstep = Int(round(T_short/dt))

T_pl = T/2
n_pl = Int(T_pl/dt)
p2 = plot(times[n_pl-shortstep:n_pl],A[n_pl-shortstep:n_pl], lw=1.5, ylabel=L"$A$ [Hz]", color=:red, xticks=nothing)
raster_times,raster_inds = get_raster_data(n_pl)
p3 = scatter(raster_times,raster_inds, xlabel=L"$t$ [s]", ylabel="Neuron index \$n\$", color=:red)

T_pl = T
n_pl = Int(T_pl/dt)
p4 = plot(times[n_pl-shortstep:n_pl],A[n_pl-shortstep:n_pl], lw=1.5, color=:magenta, xticks=nothing)
raster_times,raster_inds = get_raster_data(n_pl)
p5 = scatter(raster_times,raster_inds, xlabel=L"$t$ [s]", color=:magenta)

#display(p1)
display(plot(p2,p4,p3,p5,layout=(2,2),link=:both))


####################################### Delays

function plot_del(pl)
    n_del = length(all_delays)-1
    for n in 2:n_del
        l = min(length(times_nref),length(all_delays[n]))
        plot!(pl,times_nref[1:l],all_delays[n][1:l], markershape=:circle, markersize = 1.5, lw=0., ls=:dot, markerstrokewidth=.15)
    end
    plot!(pl,times_nref,all_delays[n_del+1],lw=2,lc=:black, label="time to next spike\nof neuron 1", legend=:topleft)
end

p1=plot([])
plot_del(p1)
#title!(p1,title )
#title!(p1,"Spike timing" )
#xlabel!(L"t_1^f")
xlabel!("firing time of neuron 1 [s]")
ylabel!(L"spike delay $\Delta t_n$ [ms]")
xlims!((0,T))

plot!([T_step,T_step],[0,13],
    #label=L"t=t_\mathrm{step}", 
    lc=:black, 
    ls=:dash)


display(p1)

