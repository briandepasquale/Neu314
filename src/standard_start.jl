using PyPlot; pygui(true);   # We'll use the plotting package PyPlot and will plot outside the plots pane
using JLD;   # for saving/loading
using Printf # for string formatting

include("general_utils.jl")   # some convenient functions

figure(1); # start a figure
clf();     # clear it
set_current_fig_position(492, 27, 308, 296); # these parameters are
    # x, y, width, and height. set them to whatever is convenient for you
