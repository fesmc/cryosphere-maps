# Ellipsoidal polar stereographic projection (Snyder 1987, pp. 160-162), in km.

struct PolarStereo
    lon0::Float64     # straight vertical longitude from pole [deg]
    lat_ts::Float64   # latitude of true scale [deg] (sign gives hemisphere)
end

const WGS84_A = 6378137.0
const WGS84_E = 0.0818191908426

const PROJ_GRL = PolarStereo(-45.0, 70.0)     # EPSG:3413
const PROJ_ANT = PolarStereo(0.0, -71.0)      # EPSG:3031

_tfun(φ, e) = tan(π/4 - φ/2) / ((1 - e*sin(φ)) / (1 + e*sin(φ)))^(e/2)
_mfun(φ, e) = cos(φ) / sqrt(1 - e^2*sin(φ)^2)

"Projected (x, y) in km of (lon, lat) in degrees."
function project(p::PolarStereo, lon, lat)
    s  = sign(p.lat_ts)                  # +1 north, -1 south
    φ  = deg2rad(s*lat)
    φc = deg2rad(s*p.lat_ts)
    λ  = deg2rad(s*(lon - p.lon0))
    e  = WGS84_E
    ρ  = WGS84_A * _mfun(φc, e) * _tfun(φ, e) / _tfun(φc, e)
    x  =  ρ*sin(λ)
    y  = -ρ*cos(λ)
    return (s*x/1e3, s*y/1e3)
end

"Geodetic latitude (deg) of projected (x, y) in km (iterative inverse)."
latitude(p::PolarStereo, x, y) = _latitude(p, hypot(x, y)*1e3)

"Latitude (deg) at distance ρ (m) from the pole."
function _latitude(p::PolarStereo, ρ)
    s  = sign(p.lat_ts); e = WGS84_E
    φc = deg2rad(s*p.lat_ts)
    t  = ρ * _tfun(φc, e) / (WGS84_A * _mfun(φc, e))
    φ  = π/2 - 2atan(t)
    for _ in 1:8
        φ = π/2 - 2atan(t * ((1 - e*sin(φ))/(1 + e*sin(φ)))^(e/2))
    end
    return s*rad2deg(φ)
end

"Projected (x, y) in km to (lon, lat) in degrees."
function unproject(p::PolarStereo, x, y)
    s = sign(p.lat_ts)
    λ = atan(s*x, -s*y)
    return (p.lon0 + s*rad2deg(λ), latitude(p, x, y))
end

"""
Linear scale factor k at projected (x, y) in km: true lengths are map lengths
divided by k, true areas map areas divided by k².
"""
function scale_factor(p::PolarStereo, x, y)
    e  = WGS84_E; s = sign(p.lat_ts)
    φc = deg2rad(s*p.lat_ts)
    ρ  = max(hypot(x, y)*1e3, 1.0)                    # avoid 0/0 at the pole
    φ  = deg2rad(s*_latitude(p, ρ))
    return ρ / (WGS84_A * _mfun(φ, e))
end
