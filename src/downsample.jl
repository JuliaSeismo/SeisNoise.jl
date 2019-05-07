export downsample

"""

  downsample(C::SeisChannel,fs::Real)

Downsample SeisChannel sampling rate to frequency `fs`.

For best results, lowpass filter data to `fs` before downsampling.
Implements the weighted average slopes interpolation scheme proposed in
[Wiggins1976] for evenly sampled data from
obspy.signal.interpolation.weighted_average_slopes.
"""
function downsample(C::SeisChannel,fs::Float64)
    dt = float(fs)

    if dt <= 0.
        throw(DomainError(dt, "Sampling rate must be nonnegative"))
    end

    dt = 1. / fs
    if fs == C.fs
        return C
    end

    old_start = C.t[1,2] * 1e-6
    old_dt = 1. / C.fs
    starttime, endtime = SeisIO.t_win(C.t,C.fs)
    starttime, endtime = starttime *1e-6, endtime * 1e-6
    npts = Int(floor((endtime-starttime)/dt)) + 1
    C.x = weighted_average_slopes(C.x,old_start,old_dt,starttime,dt,npts)
    C.fs = fs
    C.t[2,1] = npts
    return C
end

"""
    downsample(S::SeisData, fs::Float64)

Downsample SeisData sampling rate to frequency `fs`.
"""
function downsample(S::SeisData,fs::Float64)
    @inbounds for i = 1:S.n
        S[i]= downsample(S[i],fs)
    end
    return S
end

"""

  weighted_average_slopes(A::AbstractArray,old_start::Real,old_dt::Real,
                                 new_start::Real,new_dt::Real,new_npts::Int)

"""
function weighted_average_slopes(A::AbstractArray,
                                 old_start::Float64,
                                 old_dt::Float64,
                                 new_start::Float64,
                                 new_dt::Float64,
                                 new_npts::Int)
    old_end, new_end = validate_parameters(A, old_start, old_dt, new_start,
                                           new_dt, new_npts)
    new_time_array = Array(range(new_start, stop=new_end, length=new_npts))
    m = diff(A)
    m ./= old_dt
    w = abs.(m)
    w .= 1. ./ clamp.(w, eps(Float64), maximum(w))

    slope = zeros(eltype(A),length(A))
    slope[1] = m[1]
    for ii = 2:length(m)
        slope[ii] = (w[ii-1] * m[ii-1] + w[ii] * m[ii]) / (w[ii-1] + w[ii])
    end
    slope[end] = m[end]

    # If m_i and m_{i+1} have opposite signs then set the slope to zero.
    # This forces the curve to have extrema at the sample points and not
    # in-between.
    sign_change = find_nonzero(diff(sign.(m)))
    for ii = 1:length(sign_change)
        slope[sign_change[ii] + 1] = 0.
    end

    # Create interpolated value using hermite interpolation. In this case
    # it is directly applicable as the first derivatives are known.
    return hermite_interpolation(A,slope,new_time_array,old_dt,old_start)
end

"""

  hermite_interpolation(y_in, slope, x_out, h, x_start)

Hermite interpolation when zeroth and first derivatives are given for each time step.

# Arguements
- y_in: The data values to be interpolated.
- slope: The desired slope at each data point.
- x_out: The point at which to interpolate.
- h: The sample interval for y_in.
- x_start: Time of the first sample in y_in.
"""
function hermite_interpolation(y_in::AbstractArray,
                               slope::AbstractArray,
                               x_out::AbstractArray,
                               h::Float64,x_start::Float64)
    len_out = length(x_out)
    len_in = length(y_in)
    y_out = zeros(len_out)
    y_out[1] = y_in[1]

    for idx = 2:len_out
        i = (x_out[idx] - x_start) / h
        i_0 = convert(Int64,round(i))
        i_1= i_0 + 1

        if i == float(i_0)
            y_out[idx] = y_in[i_0]
            continue
        end

        t = i - float(i_0)
        a_0 = y_in[i_0]
        a_1 = y_in[i_1]
        b_minus_1 = h * slope[i_0]
        b_plus_1 = h * slope[i_1]
        b_0 = a_1 - a_0
        c_0 = b_0 - b_minus_1
        c_1 = b_plus_1 - b_0
        d_0 = c_1 - c_0

        y_out[idx] = a_0 + (b_0 + (c_0 + d_0 * t) * (t - 1.)) * t
    end
    return y_out
end

"""

  validate_parameters(data, old_start, old_dt, new_start, new_dt, new_npts)

Validates the parameters for various interpolation functions.

Returns the old and the new end.

"""
function validate_parameters(A::AbstractArray,
                             old_start::Float64,
                             old_dt::Float64,
                             new_start::Float64,
                             new_dt::Float64,
                             new_npts::Int)
    if new_dt <= 0.
        throw(DomainError(new_dt, "New sampling rate must be nonnegative"))
    end

    if ndims(A) != 1 || size(A) == (0,)
        throw(ArgumentError("Not a 1d array."))
    end

    old_end = old_start + old_dt * (length(A) - 1)
    new_end = new_start + new_dt * (new_npts - 1)

    if old_start > new_start || old_end < new_end
        throw(ArgumentError("The new array must be fully contained in the old array. No extrapolation can be performed."))
    end

    return old_end, new_end
end


function find_nonzero(c::Array{<:Union{Float32,Float64},1})
    a = similar(c, Int)
    count = 1
    @inbounds for i in eachindex(c)
        a[count] = i
        count += (c[i] != zero(eltype(c)))
    end
    return resize!(a, count-1)
end
