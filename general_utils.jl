
using Printf


"""
    tbin(tvector, t)

Returns indmin(abs.(tvector-t)) -- just a shorthand way
for finding, in a vector of time bins tvector, the bin
that corresponds the closest to time t.
"""
function tbin(tvector, t)
    return argmin(abs.(tvector.-t))
end




"""
    append_to_file(filename, str)

Opens filename, appends str to it, and closes filename.

If filename is not a string but us type Base.PipeEndpoint (an IO stream)
then simply prints to it, without trying to open or close
"""
function append_to_file(filename, str)
    if typeof(filename)<:IO
        @printf(filename, "%s", str)
    else
        fstr = open(filename, "a")
        @printf(fstr, "%s", str)
        close(fstr)
    end
end





"""
function print_vector(vec)

Takes a vector and uses @printf to put it on the screen with [%.3f, %.3f] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector(vec)
    print_vector(STDOUT, vec)
end

"""
function print_vector(fname::String, vec)

Takes a vector and uses @printf to append it to file fname with [%.3f, %.3f] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector(fname::String, vec)
    ostream = open(fname, "a")
    print_vector(ostream, vec)
    close(ostream)
end


"""
function print_vector(stream::IO, vec)

Takes a vector and uses @sprintf to put it on stream IO with [%.3f, %.3f] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector(stream::IO, vec)

    if typeof(vec)==Symbol
        mystr = string(vec)
        @printf(stream, "%s = ", mystr);
        print_vector(stream, eval(vec))
        @printf(stream, "\n");
        return
    end

    @printf stream "["
    for p in [1:length(vec);]
        @printf(stream, "%.3f", vec[p])
        if p < length(vec) @printf(stream, ", "); end
    end
    @printf(stream, "]")
end


"""
function print_vector_g(vec)

Takes a vector and uses @printf to put it on the screen with [%g, %g] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector_g(vec)
    print_vector_g(STDOUT, vec)
end


"""
function print_vector_g(fname::String, vec)

Takes a vector and uses @printf to append it to file fname with [%g, %g] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector_g(fname::String, vec)
    ostream = open(fname, "a")
    print_vector_g(ostream, vec)
    close(ostream)
end


"""
function print_vector_g(stream::IO, vec)

Takes a vector and uses @printf to put it on stream with [%g, %g] format.

If passed a symbol (which must evaluate to a vector), then prints the string for that symbol,
an equals sign, the vector, and ends by adding a carriage return \n.
"""
function print_vector_g(stream::IO, vec)

    if typeof(vec)==Symbol
        mystr = string(vec)
        @printf(stream, "%s = ", mystr);
        print_vector_g(stream, eval(vec))
        @printf(stream, "\n");
        return
    end

    @printf stream "["
    for p in [1:length(vec);]
        @printf(stream, "%g", vec[p])
        if p < length(vec) @printf(stream, ", "); end
    end
    @printf(stream, "]")
end



"""
safe_axes(axh; further_params...)

If you're going to make axh the current axes, this function
first makes axh's figure the current figure. Some Julias
break without that.

Any optional keyword-value arguments are passed on to axes()
"""

function safe_axes(axh; further_params...)
    figure(axh.figure.number)
    PyPlot.axes(axh; Dict(further_params)...)
end


"""
ax = axisWidthChange(factor; lock="c", ax=nothing)

Changes the width of the current axes by a scalar factor.

= PARAMETERS:
 - factor      The scalar value by which to change the width, for example
               0.8 (to make them thinner) or 1.5 (to make them fatter)

= OPTIONAL PARAMETERS:
 - lock="c"    Which part of the axis to keep fixed. "c", the default does
               the changes around the middle; "l" means keep the left edge fixed
               "r" means keep the right edge fixed

 - ax = nothing   If left as the default (nothing), works on the current axis;
               otherwise should be an axis object to be modified.
"""
function axisWidthChange(factor; lock="c", ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax.get_position().bounds

    if lock=="l";
    elseif lock=="c" || lock=="m"; x = x + w*(1-factor)/2;
    elseif lock=="r"; x = x + w*(1-factor);
    else error("I don't know lock type ", lock)
    end

    w = w*factor;
    ax.set_position([x, y, w, h])

    return ax
end


"""
ax = axisHeightChange(factor; lock="c", ax=nothing)

Changes the height of the current axes by a scalar factor.

= PARAMETERS:
 - factor      The scalar value by which to change the height, for example
               0.8 (to make them shorter) or 1.5 (to make them taller)

= OPTIONAL PARAMETERS:
 - lock="c"    Which part of the axis to keep fixed. "c", the default does
               the changes around the middle; "b" means keep the bottom edge fixed
               "t" means keep the top edge fixed

 - ax = nothing   If left as the default (nothing), works on the current axis;
               otherwise should be an axis object to be modified.
"""
function axisHeightChange(factor; lock="c", ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax.get_position().bounds

    if lock=="b";
    elseif lock=="c" || lock=="m"; y = y + h*(1-factor)/2;
    elseif lock=="t"; y = y + h*(1-factor);
    else error("I don't know lock type ", lock)
    end

    h = h*factor;
    ax.set_position([x, y, w, h])

    return ax
end


"""
   ax = axisMove(xd, yd; ax=nothing)

Move an axis within a figure.

= PARAMETERS:
- xd      How much to move horizontally. Units are scaled figure units, from
           0 to 1 (with 1 meaning the full width of the figure)

- yd      How much to move vertically. Units are scaled figure units, from
            0 to 1 (with 1 meaning the full height of the figure)

= OPTIONAL PARAMETERS:
 - ax = nothing   If left as the default (nothing), works on the current axis;
               otherwise should be an axis object to be modified.

"""
function axisMove(xd, yd; ax=nothing)
    if ax==nothing; ax=gca(); end
    x, y, w, h = ax.get_position().bounds

    x += xd
    y += yd

    ax.set_position([x, y, w, h])
    return ax
end


"""
[] = remove_xtick_labels(ax=NaN)

Given an axis object, or an array of axes objects, replaces each xtick label
string with the empty string "".

If no axis is passed, uses gca() to work with the current axis.


"""
function remove_xtick_labels(ax=nothing)

    if ax==nothing
        ax = gca()
    end

    if typeof(ax) <: Array
        for i=1:length(ax)
            remove_xtick_labels(ax[i])
        end
        return
    end

    nlabels = length(ax.xaxis.get_ticklabels())

    newlabels = Array{String,1}(nlabels)
    for i=1:length(newlabels);
        newlabels[i] = ""
    end

    ax.xaxis.set_ticklabels(newlabels)
    return
end



"""
[] = remove_ytick_labels(ax=NaN)

Given an axis object, or an array of axes objects, replaces each ytick label
string with the empty string "".

If no axis is passed, uses gca() to work with the current axis.


"""
function remove_ytick_labels(ax=nothing)

    if ax==nothing
        ax = gca()
    end

    if typeof(ax) <: Array
        for i=1:length(ax)
            remove_ytick_labels(ax[i])
        end
        return
    end

    nlabels = length(ax.yaxis.get_ticklabels())

    newlabels = Array{String,1}(nlabels)
    for i=1:length(newlabels);
        newlabels[i] = ""
    end

    ax.yaxis.set_ticklabels(newlabels)
    return
end





using PyCall


"""

(x, y, w, h) = get_current_fig_position()

Returns the current figure's x, y, width and height position on the screen.

Works only when pygui(true) and when the back end is Tk or QT.
Has been tested only with PyPlot.
"""
function get_current_fig_position()
    # if !contains(pystring(plt[:get_current_fig_manager]()), "FigureManagerQT")
    try
        if occursin("Tk", PyCall.pystring(PyPlot.get_current_fig_manager()))
            g = split(plt.get_current_fig_manager().window.geometry(), ['x', '+'])
            w = parse(Int64, g[1])
            h = parse(Int64, g[2])
            x = parse(Int64, g[3])
            y = parse(Int64, g[4])
        elseif occursin("QT", pystring(plt.get_current_fig_manager()))
            x = PyPlot.get_current_fig_manager().window.pos().x()
            y = PyPlot.get_current_fig_manager().window.pos().y()
            w = PyPlot.get_current_fig_manager().window.width()
            h = PyPlot.get_current_fig_manager().window.height()
        else
            error("Only know how to work with matplotlib graphics backends that are either Tk or QT")
        end

        return (x, y, w, h)
    catch
        error("Failed to get current figure position. Is pygui(false) or are you using a back end other than QT or Tk?")
    end
end

"""

set_current_fig_position(x, y, w, h)

Sets the current figure's x, y, width and height position on the screen.

Works only when pygui(true) and when the back end is Tk or QT.
Has been tested only with PyPlot.
"""
function set_current_fig_position(x, y, w, h)
    # if !contains(pystring(plt[:get_current_fig_manager]()), "FigureManagerQT")
    try
        if occursin("Tk", PyCall.pystring(PyPlot.get_current_fig_manager()))
            PyPlot.get_current_fig_manager().window.geometry(@sprintf("%dx%d+%d+%d", w, h, x, y))
        elseif occursin("QT", pystring(plt.get_current_fig_manager()))
            PyPlot.get_current_fig_manager().window.setGeometry(x, y, w, h)
        else
            error("Only know how to work with matplotlib graphics backends that are either Tk or QT")
        end
    catch
        error("Failed to set current figure position. Is pygui(false) or are you using a back end other than QT?")
    end
end


"""
    C = capture_current_figure_configuration()

Collects the positions of all current figures and
prints out to the screen code, that can be copy-pasted,
that would reproduce that positioning configuration.

# PARAMETERS:

None

# RETURNS:

- C    A matrix that is nfigures-by-5 in size. You probably
    don't want this, you probably want the text printed to
    the screen, but here just in case.  Each row will have,
    in order: figure number, x, y, width, height
"""
function capture_current_figure_configuration()
    @printf("The following code will reproduce your current figure placement:\n\n")
    C = zeros(Int64, 0,5)
    for f in sort(PyPlot.get_fignums())
        figure(f)
        x, y, w, h = get_current_fig_position()
        @printf("figure(%d); set_current_fig_position(%d, %d, %d, %d)   # x, y, width, height\n",
            f, x, y, w, h)
        C = [C ; [f x y w h]]
    end
    return C
end
