using SeisIO
using DataFrames
using Dates
using JLD2

"""
    process_raw(st,fs)

Pre-process month-long stream of data.

Checks:
- sample rate is fs
- downsamples data
- checks for gaps in data
- Trims data to first and last day of month
- phase-shifts data to begin at 00:00:00.0
- chunks data into 86,400 second traces
- removes instrument response (pole-zero)

# Arguments
- `st::SeisData`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
"""
function process_raw(st, fs)

end


function downsample(st,fs)
    dt = float(fs)

    if dt <= 0.
        throw(DomainError(dt, "Sampling rate must be nonnegative"))
    end

    dt = 1. / fs

    old_start = st.t[1,2] * 1e-6
    old_dt = 1. / s.fs
    starttime, endtime = t_win(st.t,st.fs)
    starttime, endtime = starttime *1e-6, endtime * 1e-6
    npts = Int(floor((endtime-starttime)/dt)) + 1
    st.x = weighted_average_slopes(st.x,old_start,old_dt,starttime,dt,npts)
    st.fs = fs
    st.t[2,1] = npts

end


function weighted_average_slopes(data,old_start,old_dt,new_start,new_dt,new_npts)
    old_end, new_end = validate_parameters(data, old_start, old_dt, new_start,
                                           new_dt, new_pts)
    new_time_array = Array(range(new_start, stop=new_end, length=new_npts))
    m = diff(data) / old_dt
    w = abs.(m)
    w .= 1. ./ clamp.(w, eps(Float64), maximum(w))

    slope = zeros(length(data))
    slope[1] = m[1]
    slope[2:end-1] = (w[1:end-1] .* m[1:end-1] .+ w[2:end] .* m[2:end]) ./
                     (w[1:end-1] .+ w[2:end])
    slope[end] = m[end]

    # If m_i and m_{i+1} have opposite signs then set the slope to zero.
    # This forces the curve to have extrema at the sample points and not
    # in-between.
    sign_change = diff(sign.(m))
    sign_change = Array([ind for ind,val in enumerate(sign_change) if val != 0])
    slope[sign_change .+ 1] .= 0.

    # Create interpolated value using hermite interpolation. In this case
    # it is directly applicable as the first derivatives are known.
    new_data = hermite_interpolation(data,slope,new_time_array,old_dt,old_start)
    return new_data
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
function hermite_interpolation(y_in,slope,x_out,h,x_start)
    len_out = length(x_out)
    len_in = length(data)
    y_out = empty(x_out)

    for idx in 1:len_out
        i = (x_out[idx] - x_start) / h
        i_0 = Int(i)
        i_1= i_0 + 1

        if i == float(i_0)
            y_out[idx] = y_in[i_0]
            continue
        end

        t = i - float(i_0)
        a_0 = y_in[i_0]
        a1 = y_in[i_1]
        b_minus_1 = h * slope[i_0]
        b_plus_1 = h * slope[i_1]
        b_0 = a_1 - a_0
        c_0 = b_0 - b_minus_1
        c_1 = b_plus_1 - b_0
        d_0 = c_1 - c_0

        y_out[idx] = a_0 + (b_0 + (c_0 + d_0 * t) * (t - 1.)) * t

    end
end

"""

  validate_parameters(data, old_start, old_dt, new_start, new_dt, new_npts)

Validates the parameters for various interpolation functions.

Returns the old and the new end.

"""
function validate_parameters(data, old_start, old_dt, new_start, new_dt, new_npts)
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
