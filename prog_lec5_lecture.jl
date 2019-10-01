#####  Loading and plotting the Ca++ data

using JLD; using PyPlot

# Courtesy of Sue Ann Koay,
# 2-photon Ca++ imaging recordings from four V1 neurons of a mouse
# performing the "Accumulating Towers Task"
# (Pinto, Koay, et al. Frontiers in Behavioral Neuroscience 2018)
data = load("sue_ann_koay_data.jld")

N1 = data["neuralF"][:,1]  # Fluorescence of neuron 1 as a function of frame number
# Frames are at 15 frames/sec
Lcue = data["LcueFrame"];  # Frame numbers at which a Left tower appeared
# Lcue = convert(Array{Int64,1}, Lcue[:])
Rcue = data["RcueFrame"];  # Frame numbers at which a Right tower appeared
# Rcue = convert(Array{Int64,1}, Rcue[:])


figure(1); clf();
plot(1:length(N1), N1)
vlines(Rcue, -2, -1, color="green", linewidth=4)
vlines(Lcue, -1,  0, color="red", linewidth=4)
xlim(1900, 2050)

#####

"""
find_mu_sd(V)

Given a vector V, returns two values, its mean and its standard deviation

= PARAMETERS:

- V    vector of numbers. Strings not allowed

= RETURNS:

-  mu    mean of the elements of the vector
-  sd    standard deviation of the elements of the vector

"""
function find_mu_sd(V)
    # Always good to have clear error messages that will help people reading them
    if typeof(V)<:String
        error("find_mu_sd needs numeric vectors, can't take strings")
    end

    mumu = 0; # running tally of summed elements of V
    ss = 0;   # running tally of summed of squares of the elements of V
    # Now loop over each element of V, from first to last:
    for zuzu=1:length(V)
        mumu = mumu + V[zuzu]
        ss   = ss   + V[zuzu]^2
    end
    mumu = mumu/length(V)   # average value
    ss   = ss/length(V)     # average squared value

    # Formula for unbiased standard deviation:
    sd = sqrt((ss - mumu^2)*(length(V)/(length(V)-1)))

    # return both mean and standard deviation
    return mumu, sd
end

# EXAMPLE CALLS:
mu, sd = find_mu_sd([1, 3, 55.6, 2])

sd = find_mu_sd([1, 3, -90.8, 9])[2]



######

"""
find_big_peaks(V; peakwidth=1)

finds peak values in the vector V that are also higher than 2*standard deviation of V

= PARAMETERS:

- V    a vector of numbers

= OPTIONAL PARAMETERS:

- peakwidth=1   An integer, indicates how wide a peak has to be. For the ith
                slot of V to be considered a peak, it must be larger than
                all of the peakwidth slots before it, and also larger than
                all of the peakwidth slots after it.

= RETURNS:

- peak_indices  a vector containing the slot numbers of all the peaks found. One slot per peak, at the peak.
                If no peaks above threshold are found, this will be an empty vector
                i.e., a vector with zero slots.

"""
function find_big_peaks(V; peakwidth=1)
    sd = find_mu_sd(V)[2]
    thresh = 2*sd

    list_of_frames = falses(size(V))
    for i=peakwidth+1:length(V)-peakwidth
        if (V[i] > thresh) &
            all(V[i].>V[i-peakwidth:i-1]) &
            all(V[i] .> V[i+1:i+peakwidth])
            list_of_frames[i] = true
        end
    end

    return findall(list_of_frames)
end


######

using Statistics

"""
remove_chunks!(V, chunk_indices; remove_width=4)

Mutating function that repleaces chunks of slots in V with NaNs.

= PARAMETERS:

- V     a vector of Float64 that will be mutated.
        Use copy() to copy it first if you don't want to lose the original

- chunk_indices    a list of the slot numbers in V,
                   around which a chunk of slots will be replaced with NaNs

= OPTIONAL PARAMETERS:

- remove_width     number of slots before and number of slots after each chunk index
                   to replace with NaNs

- RETURNS:

- the mutated V

= EXAMPLE CALL:

V = [1.3 2 3 4 5]
remove_chunks!(V, 3, remove_width=1)
println(V)

1×5 Array{Float64,2}:
1.3  NaN  NaN  NaN  5.0
"""
function remove_chunks!(V, chunk_indices; remove_width=4)
    if !(typeof(V)<:Array{Float64})
        error("V must be a vector with elements that are each Float64")
    end
    for i=1:length(chunk_indices)
        start = maximum([1, chunk_indices[i]-remove_width])
        stop  = minimum([chunk_indices[i]+remove_width, length(V)])
        V[start:stop] .= NaN
    end

    return V
end

#####

# Find the slot numbers (indices) of the top of the peaks in N1:
pi = find_big_peaks(N1, peakwidth=5)

NN = copy(N1)
remove_chunks!(NN, pi, remove_width=5)

figure(1); clf();
plot(1:length(N1), N1)
vlines(Rcue, -2, -1, color="green", linewidth=4)
vlines(Lcue, -1,  0, color="red", linewidth=4)

plot(1:length(NN), NN, color="grey", linewidth=8)
xlim(1900, 2050)



#####


function stacker(V, indices; n_before=10, n_after=15)
    stack = zeros(length(indices), n_before+n_after+1)

    for i=1:length(indices)
        my_slot = indices[i]
        stack[i,:] = V[my_slot-n_before : my_slot+n_after]
    end

    return stack
end


###  Find and plot the peak-triggered average

using Statistics

n_before = 10;
n_after  = 15;

pi = find_big_peaks(N1, peakwidth=5)

stack = stacker(N1, pi, n_before=n_before, n_after=n_after)

mu = mean(stack, dims=1)[:]
sd = std(stack, dims=1)[:]

figure(2); clf();
plot(-n_before:n_after, mu)
plot(-n_before:n_after, mu[:] - sd[:]/sqrt(length(pi)), "--", color="grey")
plot(-n_before:n_after, mu[:] + sd[:]/sqrt(length(pi)), "--", color="grey")
title("peak-triggered average")
xlabel("frame number")
ylabel("fluorescence")


### find and plot the left-cue triggered average


Lstack = stacker(N1, Lcue)
figure(1); clf();
imshow(Lstack, extent=[-n_before-0.5, n_after+0.5, 0.5, size(Lstack)[1]+0.5]); axis("auto")

xlabel("frame number")
ylabel("cue number")
title("Left cue-triggered fluorescence")

mu = mean(Lstack, dims=1)[:]
sd = std(Lstack, dims=1)[:]

figure(2); clf();
plot(-n_before:n_after, mu)
plot(-n_before:n_after, mu[:] - sd[:]/sqrt(length(pi)), "--", color="grey")
plot(-n_before:n_after, mu[:] + sd[:]/sqrt(length(pi)), "--", color="grey")
title("Rcue-triggered average")
xlabel("frame number")
ylabel("fluorescence")


### find and plot the right-cue triggered average


Rstack = stacker(N1, Rcue)
figure(1); clf();
imshow(Rstack, extent=[-n_before-0.5, n_after+0.5, 0.5, size(Rstack)[1]+0.5]); axis("auto")

xlabel("frame number")
ylabel("cue number")
title("Right cue-triggered fluorescence")

mu = mean(Rstack, dims=1)[:]
sd = std(Rstack, dims=1)[:]

figure(2); clf();
plot(-n_before:n_after, mu)
plot(-n_before:n_after, mu[:] - sd[:]/sqrt(length(pi)), "--", color="grey")
plot(-n_before:n_after, mu[:] + sd[:]/sqrt(length(pi)), "--", color="grey")
title("Rcue-triggered average")
xlabel("frame number")
ylabel("fluorescence")
