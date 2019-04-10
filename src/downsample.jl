export downsample

"""

  downsample(C::SeisChannel,fs::Real)

Downsample SeisChannel sampling rate to frequency `fs`.

For best results, lowpass filter data to `fs` before downsampling.
Implements the weighted average slopes interpolation scheme proposed in
[Wiggins1976] for evenly sampled data from
obspy.signal.interpolation.weighted_average_slopes.
"""
function downsample(C::SeisChannel,fs::Real)
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
    starttime, endtime = t_win(C.t,C.fs)
    starttime, endtime = starttime *1e-6, endtime * 1e-6
    npts = Int(floor((endtime-starttime)/dt)) + 1
    C.x = weighted_average_slopes(C.x,old_start,old_dt,starttime,dt,npts)
    C.fs = fs
    C.t[2,1] = npts
    return C
end

"""
    downsample(S::SeisData)

Downsample SeisData sampling rate to frequency `fs`.
"""
function downsample(S::SeisData,fs::Real)
    @inbounds for i = 1:S.n
        S[i]=downsample(S[i],fs)
    end
    return S
end

"""

  weighted_average_slopes(A::AbstractArray,old_start::Real,old_dt::Real,
                                 new_start::Real,new_dt::Real,new_npts::Int)

"""
function weighted_average_slopes(A::AbstractArray,old_start::Real,old_dt::Real,
                                 new_start::Real,new_dt::Real,new_npts::Int)
    old_end, new_end = validate_parameters(A, old_start, old_dt, new_start,
                                           new_dt, new_npts)
    new_time_array = Array(range(new_start, stop=new_end, length=new_npts))
    m = diff(A) / old_dt
    w = abs.(m)
    w .= 1. ./ clamp.(w, eps(Float64), maximum(w))

    slope = zeros(length(A))
    slope[1] = m[1]
    slope[2:end-1] = (w[1:end-1] .* m[1:end-1] .+ w[2:end] .* m[2:end]) ./ (w[1:end-1] .+ w[2:end])
    slope[end] = m[end]

    # If m_i and m_{i+1} have opposite signs then set the slope to zero.
    # This forces the curve to have extrema at the sample points and not
    # in-between.
    sign_change = findall(x -> x != 0, diff(sign.(m)))
    slope[2:end-1][sign_change] .= 0.

    # Create interpolated value using hermite interpolation. In this case
    # it is directly applicable as the first derivatives are known.
    A = hermite_interpolation(A,slope,new_time_array,old_dt,old_start)
    return A
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
function hermite_interpolation(y_in::AbstractArray,slope::AbstractArray,
                               x_out::AbstractArray,h::Real,x_start::Real)
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
function validate_parameters(data::AbstractArray, old_start::Real, old_dt::Real,
                             new_start::Real, new_dt::Real, new_npts::Int)
    if new_dt <= 0.
        throw(DomainError(new_dt, "New sampling rate must be nonnegative"))
    end

    if ndims(data) != 1 || size(data) == (0,)
        throw(ArgumentError("Not a 1d array."))
    end

    old_end = old_start + old_dt * (length(data) - 1)
    new_end = new_start + new_dt * (new_npts - 1)

    if old_start > new_start || old_end < new_end
        throw(ArgumentError("The new array must be fully contained in the old array. No extrapolation can be performed."))
    end

    return old_end, new_end
end
