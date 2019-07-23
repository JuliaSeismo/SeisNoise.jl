export angular_distance, azimuth, angular_step, surface_distance, forward
export inverse, get_dist, get_azi, get_baz, get_loc

"""
angular_distance, azimuth, angular_step, surface_distance, forward and inverse
functions were written by Andy Nowacki in his (unregistered) Geodesics.jl
package https://github.com/anowacki/Geodesics.jl
"""

"Earth ellipsoid semi-axes in WGS84"
const EARTH_R_MAJOR_WGS84 = 6378137.0000
const EARTH_R_MINOR_WGS84 = 6356752.3142
"Flattening of the Earth in WGS84"
const F_WGS84 = (EARTH_R_MAJOR_WGS84 - EARTH_R_MINOR_WGS84)/EARTH_R_MAJOR_WGS84

"""
    angular_distance(lon0, lat0, lon1, lat1, degrees=true; f=F_WGS84) -> Δ
Return the angular distance between points (`lon0`,`lat0`) and (`lon1`,`lat1`)
on a flattened sphere.  The default flattening is $F_WGS84.  By default,
input and output are in degrees, but specify `degrees` as `false` to use radians.
Use `f=0` to perform computations on a sphere.
"""
function angular_distance(lon0, lat0, lon1, lat1, degrees::Bool=true; f=F_WGS84)
    if degrees
        lon0, lat0, lon1, lat1 = deg2rad(lon0), deg2rad(lat0), deg2rad(lon1), deg2rad(lat1)
    end
    gcarc, az, baz = inverse(lon0, lat0, lon1, lat1, 1.0, f)
    degrees ? rad2deg(gcarc) : gcarc
end

"""
    azimuth(lon0, lat0, lon1, lat1, degrees=true, f=F_WGS84)
Return the azimuth from point (`lon0`,`lat0`) to point (`lon1`,`lat1`) on a
flattened sphere.  The default flattening is $F_WGS84.  By default,
input and output are in degrees, but specify `degrees` as `false` to use radians.
Use `f=0` to perform computations on a sphere.
"""
function azimuth(lon0::T1, lat0::T2, lon1::T3, lat1::T4, degrees::Bool=true; f=F_WGS84) where {T1,T2,T3,T4}
    # Promote to Float64 at least to avoid error when one of the points is on the pole
    T = promote_type(Float64, T1, T2, T3, T4)
    lon0, lat0, lon1, lat1 = T(lon0), T(lat0), T(lon1), T(lat1)
    if degrees
        lon0, lat0, lon1, lat1 = deg2rad(lon0), deg2rad(lat0), deg2rad(lon1), deg2rad(lat1)
    end
    gcarc, az, baz = inverse(lon0, lat0, lon1, lat1, 1.0, f)
    degrees ? rad2deg(az) : az
end

"""
    step(lon, lat, azimuth, distance, degrees=false; f=F_WGS84) -> lon′, lat′, backazimuth
Return the longitude `lon′`, latitude `lat′` and backazimuth `baz` reached by
travelling an angular `distance` along `azimuth` from the starting point at
(`lon`,`lat`).  The default flattening is $F_WGS84.  By default,
input and output are in degrees, but specify `degrees` as `false` to use radians.
Use `f=0` to perform computations on a sphere.
"""
function angular_step(lon, lat, azimuth, distance, degrees::Bool=true; f=F_WGS84)
    if degrees
        lon, lat, azimuth, distance = deg2rad(lon), deg2rad(lat), deg2rad(azimuth), deg2rad(distance)
    end
    lon′, lat′, baz = forward(lon, lat, azimuth, distance, 1.0, f)
    if degrees
        rad2deg(lon′), rad2deg(lat′), rad2deg(baz)
    else
        lon′, lat′, baz
    end
end

"""
    surface_distance(lon0, lat0, lon1, lat1, a, degrees::Bool=false; f=F_WGS84)
Return the physical distance between points (`lon0`,`lat0`) and (`lon1`,`lat1`) on
the flattened sphere with flattening `f` and semimajor radius `a`.  Distance is given
in the same units as `a`.  By default, input angles are in degrees, but specify `degrees`
as `false` to use radians.  Use `f=0` to perform computations on a sphere.
"""
function surface_distance(lon0, lat0, lon1, lat1, a, degrees::Bool=true; f=F_WGS84)
    if degrees
        lon0, lat0, lon1, lat1 = deg2rad(lon0), deg2rad(lat0), deg2rad(lon1), deg2rad(lat1)
    end
    distance, az, baz = inverse(lon0, lat0, lon1, lat1, a, f)
    distance
end

"""
    forward(lon, lat, azimuth, distance, a, f) -> lon′, lat′, backazimuth
Return the longitude `lon′` and latitude `lat′` and `backazimuth` of a projected
point, reached by travelling along an `azimuth` for `distance` from an original
point at (`lon`, `lat`).  Specify the spheroid with the semimajor radius `a` and
flattening `f`.
Coordinates and azimuth are in radians.
Calculations use Vincenty's forward formula [1].
#### References
1. Vincenty, T. (1975). "Direct and Inverse Solutions of Geodesics on the Ellipsoid
   with application of nested equations" (PDF). Survey Review. XXIII (176): 88–93.
   doi:10.1179/sre.1975.23.176.88
"""
function forward(lon, lat, azimuth, distance, a, f)::Tuple{Float64,Float64,Float64}
    abs(lat <= π/2) || throw(ArgumentError("Latitude ($lat) must be in range [-π/2, π/2]"))
    a > 0 || throw(ArgumentError("Semimajor axis ($a) must be positive"))
    abs(f) < 1 || throw(ArgumentError("Magnitude of flattening ($f) must be less than 1"))
    # Calculations are done with Float64s internally as the tolerances are hard-wired
    lambda1, phi1, alpha12, s = Float64(lon), Float64(lat), Float64(azimuth), Float64(distance)
    a, f = Float64(a), Float64(f)
    alpha12 = mod(alpha12, 2π)
    b = a*(1 - f)

    TanU1 = (1 - f)*tan(phi1)
    U1 = atan(TanU1)
    sigma1 = atan( TanU1, cos(alpha12) )
    Sinalpha = cos(U1)*sin(alpha12)
    cosalpha_sq = 1.0 - Sinalpha*Sinalpha

    u2 = cosalpha_sq*(a*a - b*b )/(b*b)
    A = 1.0 + (u2/16384)*(4096 + u2*(-768 + u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))

    # Starting with the approximation
    sigma = (s/(b*A))

    # Not moving anywhere. We can return the location that was passed in.
    if sigma == 0
        return phi1, lambda1, alpha12
    end

    last_sigma = 2*sigma + 2 # something impossible

    # Iterate the following three equations
    # until there is no significant change in sigma
    # two_sigma_m , delta_sigma
    while abs((last_sigma - sigma)/sigma) > 1.0e-9
        global two_sigma_m = 2*sigma1 + sigma
        delta_sigma = B*sin(sigma)*(cos(two_sigma_m) + (B/4)*(cos(sigma)*(-1 + 2*cos(two_sigma_m)^2 - (B/6)*cos(two_sigma_m)*(-3 + 4*sin(sigma)^2)*(-3 + 4*cos(two_sigma_m)^2 ))))
        last_sigma = sigma
        sigma = (s/(b*A)) + delta_sigma
    end

    phi2 = atan((sin(U1)*cos(sigma) + cos(U1)*sin(sigma)*cos(alpha12)),
        ((1-f)*sqrt(Sinalpha^2 + (sin(U1)*sin(sigma) - cos(U1)*cos(sigma)*cos(alpha12))^2)))

    lambda = atan((sin(sigma)*sin(alpha12)),
        (cos(U1)*cos(sigma) - sin(U1)*sin(sigma)*cos(alpha12)))

    C = (f/16)*cosalpha_sq*(4 + f*(4 - 3*cosalpha_sq))

    omega = lambda - (1-C)*f*Sinalpha*(sigma + C*sin(sigma)*(
        cos(two_sigma_m) + C*cos(sigma)*(-1 + 2*cos(two_sigma_m)^2)))

    lambda2 = lambda1 + omega

    alpha21 = atan(Sinalpha, (-sin(U1)*sin(sigma) + cos(U1)*cos(sigma)*cos(alpha12)))
    alpha21 = mod(alpha21 + π, 2π)

    return lambda2, phi2, alpha21
end

"""
    inverse(lon1, lat1, lon2, lat2, a, f) -> distance, azimuth, backazimuth
Return the `distance`, `azimuth` and `backazimuth` between two points with longitudes
`lon1` and `lon2`, and latitudes `lat1` and `lat2`.  Specify the spheroid with the
semimajor radius `a` and flattening `f`.
Coordinates and angles are in radians, whilst `distance` is in the same units as `a`.
Calculations use Vincenty's inverse formula [1].
#### References
1. Vincenty, T. (1975). "Direct and Inverse Solutions of Geodesics on the Ellipsoid
   with application of nested equations" (PDF). Survey Review. XXIII (176): 88–93.
   doi:10.1179/sre.1975.23.176.88
"""
function inverse(lon1, lat1, lon2, lat2, a, f)::Tuple{Float64,Float64,Float64}
    for lat in (lat1, lat2)
        abs(lat <= π/2) || throw(ArgumentError("Latitude ($lat) must be in range [-π/2, π/2]"))
    end
    a > 0 || throw(ArgumentError("Semimajor axis ($a) must be positive"))
    abs(f) < 1 || throw(ArgumentError("Magnitude of flattening ($f) must be less than 1"))
    lambda1, phi1, lambda2, phi2 = Float64(lon1), Float64(lat1), Float64(lon2), Float64(lat2)
    a, f = Float64(a), Float64(f)
    tol = 1.0e-8
    if (abs(phi2 - phi1) < tol) && (abs(lambda2 - lambda1) < tol)
        return 0.0, 0.0, 0.0
    end

    b = a*(1 - f)

    TanU1 = (1 - f)*tan(phi1)
    TanU2 = (1 - f)*tan(phi2)

    U1 = atan(TanU1)
    U2 = atan(TanU2)

    lambda = lambda2 - lambda1
    last_lambda = -4000000.0 # an impossibe value
    omega = lambda

    # Iterate the following equations until there is no significant change in lambda
    alpha, sigma, Sin_sigma, Cos2sigma_m, Cos_sigma, sqr_sin_sigma =
        -999999., -999999., -999999., -999999., -999999., -999999.
    while ((last_lambda < -3000000.0) || (lambda != 0)) &&
            (abs((last_lambda - lambda)/lambda) > 1.0e-9)
        sqr_sin_sigma = (cos(U2)*sin(lambda))^2 +
                         ((cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda)))^2
        Sin_sigma = sqrt(sqr_sin_sigma)
        Cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda)
        sigma = atan(Sin_sigma, Cos_sigma)

        Sin_alpha = cos(U1)*cos(U2)*sin(lambda)/sin(sigma)

        if Sin_alpha >= 1
            Sin_alpha = 1.0
        elseif Sin_alpha <= -1
            Sin_alpha = -1.0
        end

        alpha = asin(Sin_alpha)
        Cos2sigma_m = cos(sigma) - 2*sin(U1)*sin(U2)/cos(alpha)^2
        C = (f/16)*cos(alpha)^2*(4 + f*(4 - 3*cos(alpha)^2))
        last_lambda = lambda
        lambda = omega + (1 - C)*f*sin(alpha)*(sigma +
            C*sin(sigma)*(Cos2sigma_m + C*cos(sigma)*(-1 + 2*Cos2sigma_m^2)))
    end

    u2 = cos(alpha)^2*(a*a - b*b)/(b*b)
    A = 1 + (u2/16384)*(4096 + u2*(-768 + u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
    delta_sigma = B*Sin_sigma*(Cos2sigma_m + (B/4)*(
        Cos_sigma*(-1 + 2*Cos2sigma_m^2) -
        (B/6)*Cos2sigma_m*(-3 + 4*sqr_sin_sigma)*(-3 + 4*Cos2sigma_m^2)))
    s = b*A*(sigma - delta_sigma)

    alpha12 = atan((cos(U2)*sin(lambda)), ( cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda)))
    alpha21 = atan((cos(U1)*sin(lambda)), (-sin(U1)*cos(U2) + cos(U1)*sin(U2)*cos(lambda)))

    alpha12 = mod(alpha12, 2π)
    alpha21 = mod(alpha21 + π, 2π)

    return s, alpha12, alpha21
end

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
