"""
    module neu314

A julia package for the best course in PNI!
"""
module neu314

using PyCall, Printf, PyPlot, JLD

include("standard_start.jl")
include("general_utils.jl")

export tbin, append_to_file, print_vector
export print_vector_g, safe_axes, axisWidthChange
export axisHeightChange, axisMove, remove_xtick_labels
export remove_ytick_labels, get_current_fig_position
export set_current_fig_position, capture_current_figure_configuration

end
