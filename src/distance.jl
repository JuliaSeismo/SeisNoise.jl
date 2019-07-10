export get_dist, get_azi, get_baz, get_loc

"""
  get_dist(loc1,loc2)

Get distance bewteen source `loc1` and receiver `loc2` in km.
If either location is empty, return zero.
# Arguments
- `loc1::GeoLoc`: Location of source.
- `loc2::GeoLoc`: Location of receiver.
"""
function get_dist(loc1::GeoLoc, loc2::GeoLoc)
    if isempty(loc1) | isempty(loc2)
        dist =  zero(Float64)
    else
        dist =  Geodesics.surface_distance(loc1.lon,loc1.lat,loc2.lon,loc2.lat,
                               Geodesics.EARTH_R_MAJOR_WGS84/1000)
    end
    return dist
end

"""
  get_dist(FFT1,FFT2)

Get distance bewteen source `FFT1` and receiver `FFT2` in km.
If either location is empty, return zero.
# Arguments
- `FFT1::FFTData`: FFTData with source location.
- `FFT2::FFTData`: FFTData with receiver location.
"""
function get_dist(FFT1::FFTData, FFT2::FFTData)
    return get_dist(FFT1.loc,FFT2.loc)
end

"""
  get_azi(loc1,loc2)

Get azimuth from source `loc1` to receiver `loc2` in degrees.
If either location is empty, return zero.
# Arguments
- `loc1::GeoLoc`: Location of source.
- `loc2::GeoLoc`: Location of receiver.
"""
function get_azi(loc1::GeoLoc, loc2::GeoLoc)
    if isempty(loc1) | isempty(loc2)
        azi =  zero(Float64)
    else
        azi = Geodesics.azimuth(loc1.lon,loc1.lat,loc2.lon,loc2.lat)
    end
    return azi
end

"""
  get_azi(FFT1,FFT2)

Get azimuth from source `FFT1` to receiver `FFT2` in degrees.
If either location is empty, return zero.
# Arguments
- `FFT1::FFTData`: FFTData with source location.
- `FFT2::FFTData`: FFTData with receiver location.
"""
function get_azi(FFT1::FFTData, FFT2::FFTData)
    return get_azi(FFT1.loc,FFT2.loc)
end

"""
  get_baz(loc1,loc2)

Get back azimuth from source `loc1` to receiver `loc2` in degrees.
If either location is empty, return zero.
# Arguments
- `loc1::GeoLoc`: Location of source.
- `loc2::GeoLoc`: Location of receiver.
"""
function get_baz(loc1::GeoLoc, loc2::GeoLoc)
    if isempty(loc1) | isempty(loc2)
        baz =  zero(Float64)
    else
        baz = Geodesics.inverse(deg2rad.((loc1.lon,loc1.lat,loc2.lon,loc2.lat))...,
                          Geodesics.EARTH_R_MAJOR_WGS84, Geodesics.F_WGS84)[3]
        baz = rad2deg(baz)
    end
    return baz
end

"""
  get_baz(FFT1,FFT2)

Get back azimuth from source `FFT1` to receiver `FFT2` in degrees.
If either location is empty, return zero.
# Arguments
- `FFT1::FFTData`: FFTData with source location.
- `FFT2::FFTData`: FFTData with receiver location.
"""
function get_baz(FFT1::FFTData, FFT2::FFTData)
    return get_baz(FFT1.loc,FFT2.loc)
end

"""
  get_loc(loc1,azi, dist)

Get location of receiver given location of source and dist + azi to receiver.
If location is empty, return empty location.
# Arguments
- `loc1::GeoLoc`: Location of source.
- `azi::Float64`: Azimuth from source `loc1` to receiver.
- `dist::Float64`: Distance from source `loc1` to receiver in km.
"""
function get_loc(loc1::GeoLoc, azi::Float64, dist::Float64)
    if isempty(loc1)
        loc2 = GeoLoc()
    else
        lon2,lat2,baz2 = Geodesics.forward(deg2rad.((loc1.lon, loc1.lat,azi))...,
                         dist*1000, Geodesics.EARTH_R_MAJOR_WGS84, Geodesics.F_WGS84)
        lon2, lat2, baz2 = rad2deg.((lon2, lat2, baz2))
        loc2 = GeoLoc("",lat2,lon2, zero(Float64), zero(Float64), zero(Float64),
                      loc1.inc)
    end
    return loc2
end

"""
  get_loc(C)

Get location of receiver and source from CorrData `C`.
If either location is empty, return empty locations.
# Arguments
- `C::CorrData`: CorrData with source location and distance + azimuth to receiver.
"""
function get_loc(C::CorrData)
    if isempty(C.loc)
        loc2 = GeoLoc()
    else
        loc2 = get_loc(C.loc, C.azi, C.dist)
    end

    return C.loc, loc2
end
