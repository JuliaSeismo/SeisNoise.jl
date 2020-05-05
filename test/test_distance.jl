# test distance functions

@testset "Forward" begin
    @test_throws ArgumentError forward(0, 2π, 0, 0, 1, 0)
    @test_throws ArgumentError forward(0, 2π, 0, 0, 0, 0)
    @test_throws ArgumentError inverse(0, 0, 1, 1, 1, 1)

    # Compare to Haversine computation on the unit sphere
    let lon0 = 0, lat0 = 0, azimuth = 1, distance = 1, a = 1, f = 0
        lon1, lat1, baz = forward(lon0, lat0, azimuth, distance, a, f)
        @test lon1 ≈ 0.9189892552937267
        @test lat1 ≈ 0.4719777676633856
        @test baz ≈ mod(-1.9047280216911358, 2π) atol=1e-6
    end
end

@testset "Inverse" begin
    @test_throws ArgumentError inverse(0, 0, 1, 1, 0, 0)
    @test_throws ArgumentError inverse(0, 2π, 0, -3π, 1, 1)
    @test_throws ArgumentError inverse(0, 0, 1, 1, 1, 1)

    # Compare to Haversine computation on the unit sphere
    let lon0 = 0, lat0 = 0, lon1 = 1, lat1 = 1, a = 1, f = 0, atol = 1e-6
        distance, azimuth, backazimuth = inverse(lon0, lat0, lon1, lat1, a, f)
        @test distance ≈ 1.2745557823062943 atol=atol
        @test azimuth ≈ 0.49536728921867335 atol=atol
        @test backazimuth ≈ mod(-2.0661636160135703, 2π) atol=atol
    end
end

@testset "Ang dist" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, lon1 = 360*rand() - 180,
            lat1 = 180*rand() - 90, atol = 1e-6
        @test angular_distance(lon0, lat0, lon1, lat1, f=0) ≈
            rad2deg(inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, 0)[1]) atol=atol
        @test angular_distance(lon0, lat0, lon1, lat1) ≈
            rad2deg(inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, SeisNoise.F_WGS84)[1]) atol=atol
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), lon1 = 2π*(rand() - 0.5),
            lat1 = π*(rand() - 0.5), atol = 1e-6
        @test angular_distance(lon0, lat0, lon1, lat1, false, f=0) ≈
            inverse(lon0, lat0, lon1, lat1, 1, 0)[1] atol=atol
        @test angular_distance(lon0, lat0, lon1, lat1, false) ≈
            inverse(lon0, lat0, lon1, lat1, 1, SeisNoise.F_WGS84)[1] atol=atol
    end
end

@testset "Azimuth" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, lon1 = 360*rand() - 180,
            lat1 = 180*rand() - 90, atol = 1e-6
        @test azimuth(lon0, lat0, lon1, lat1, f=0) ≈
            rad2deg(inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, 0)[2]) atol=atol
        @test azimuth(lon0, lat0, lon1, lat1) ≈
            rad2deg(inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, SeisNoise.F_WGS84)[2]) atol=atol
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), lon1 = 2π*(rand() - 0.5),
            lat1 = π*(rand() - 0.5), atol = 1e-6
        @test azimuth(lon0, lat0, lon1, lat1, false, f=0) ≈
            inverse(lon0, lat0, lon1, lat1, 1, 0)[2] atol=atol
        @test azimuth(lon0, lat0, lon1, lat1, false) ≈
            inverse(lon0, lat0, lon1, lat1, 1, SeisNoise.F_WGS84)[2] atol=atol
    end

    # One point at pole
    let lon0 = 0, lat0 = -90, lon1 = 0, lat1 = 0
        for T in (Int16, Int32, Int64, Float16, Float32, Float64)
            @test azimuth(T(lon0), T(lat0), T(lon1), T(lat1), true) == 0.0
        end
    end
    let lon0 = 90, lat0 = 0, lon1 = 90, lat1 = 90
        for T in (Int16, Int32, Int64, Float16, Float32, Float64)
            @test azimuth(T(lon0), T(lat0), T(lon1), T(lat1), true) == 0.0
        end
    end
end

@testset "Ang step" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, dist = 360*rand(),
            az = 360*rand() - 90, atol = 1e-6
        @test all(isapprox.(angular_step(lon0, lat0, az, dist, f=0),
            rad2deg.(forward(deg2rad.((lon0, lat0, az, dist))..., 1, 0)), atol=atol))
        @test all(isapprox.(angular_step(lon0, lat0, az, dist),
            rad2deg.(forward(deg2rad.((lon0, lat0, az, dist))..., 1, SeisNoise.F_WGS84)), atol=atol))
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), dist = 2π*rand(),
            az = 2π*(rand() - 0.5), atol = 1e-6
        @test all(isapprox.(angular_step(lon0, lat0, az, dist, false, f=0),
            forward(lon0, lat0, az, dist, 1, 0), atol=atol))
        @test all(isapprox.(angular_step(lon0, lat0, az, dist, false),
            forward(lon0, lat0, az, dist, 1, SeisNoise.F_WGS84), atol=atol))
    end
end

@testset "Surf dist" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, lon1 = 360*rand() - 180,
            lat1 = 180*rand() - 90, a = 10_000_000*rand(), atol = 1e-6
        @test surface_distance(lon0, lat0, lon1, lat1, a, f=0) ≈
            inverse(deg2rad.((lon0, lat0, lon1, lat1))..., a, 0)[1] atol=atol
        @test surface_distance(lon0, lat0, lon1, lat1, a) ≈
            inverse(deg2rad.((lon0, lat0, lon1, lat1))..., a, SeisNoise.F_WGS84)[1] atol=atol
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), lon1 = 2π*(rand() - 0.5),
            lat1 = π*(rand() - 0.5), a = 10_000_000*rand(), atol = 1e-6
        @test surface_distance(lon0, lat0, lon1, lat1, a, false, f=0) ≈
            inverse(lon0, lat0, lon1, lat1, a, 0)[1] atol=atol
        @test surface_distance(lon0, lat0, lon1, lat1, a, false) ≈
            inverse(lon0, lat0, lon1, lat1, a, SeisNoise.F_WGS84)[1] atol=atol
    end
end

@testset "NoiseData Distance" begin
    # test get_dist
    loc1 = GeoLoc()
    loc2 = GeoLoc()
    @test get_dist(loc1,loc2) == 0

    # second location empty
    loc1.lat = rand(-90:90,1)[1]
    loc1.lon = rand(-180:180,1)[1]
    @test get_dist(loc1,loc2) == 0

    # first location empty
    loc1 = GeoLoc()
    loc2.lat = rand(-90:90,1)[1]
    loc2.lon = rand(-180:180,1)[1]
    @test get_dist(loc1,loc2) == 0

    # both locations contain coordinates
    loc1.lat = rand(-90:90,1)[1]
    loc1.lon = rand(-180:180,1)[1]
    # check dist in km is less than 1/2 Earth circumference
    @test 0 < get_dist(loc1,loc2) < SeisNoise.EARTH_R_MAJOR_WGS84/1000 *  π

    # test get_dist with FFTData
    F1 = FFTData()
    F1.loc = loc1
    F2 = FFTData()
    F2.loc = loc2
    @test 0 < get_dist(F1,F2) < SeisNoise.EARTH_R_MAJOR_WGS84/1000 *  π

    # test get_azi
    loc1 = GeoLoc()
    loc2 = GeoLoc()
    @test get_azi(loc1,loc2) == 0

    # second location empty
    loc1.lat = rand(-90:90,1)[1]
    loc1.lon = rand(-180:180,1)[1]
    @test get_azi(loc1,loc2) == 0

    # first location empty
    loc1 = GeoLoc()
    loc2.lat = rand(-90:90,1)[1]
    loc2.lon = rand(-180:180,1)[1]
    @test get_azi(loc1,loc2) == 0

    # both locations contain coordinates
    loc1.lat = rand(-90:90,1)[1]
    loc1.lon = rand(-180:180,1)[1]
    @test 0 <= get_azi(loc1,loc2) <= 360

    # test get_azi with FFTData
    F1 = FFTData()
    F1.loc = loc1
    F2 = FFTData()
    F2.loc = loc2
    @test 0 <= get_azi(F1,F2) <= 360

    # test get_baz
    loc1 = GeoLoc()
    loc2 = GeoLoc()
    @test get_baz(loc1,loc2) == 0

    # second location empty
    loc1.lat = rand(-90:90,1)[1]
    loc1.lon = rand(-180:180,1)[1]
    @test get_baz(loc1,loc2) == 0

    # first location empty
    loc1 = GeoLoc()
    loc2.lat = rand(-90:90,1)[1]
    loc2.lon = rand(-180:180,1)[1]
    @test get_baz(loc1,loc2) == 0

    # both locations contain coordinates
    loc1.lat = rand(-90:90,1)[1]
    loc1.lon = rand(-180:180,1)[1]
    @test 0 <= get_baz(loc1,loc2) <= 360

    # test get_azi with FFTData
    F1 = FFTData()
    F1.loc = loc1
    F2 = FFTData()
    F2.loc = loc2
    @test 0 <= get_baz(F1,F2) <= 360

    # test get_loc
    azi = rand(0:360.,1)[1]
    dist = rand(1:SeisNoise.EARTH_R_MAJOR_WGS84/1000 *  π,1)[1]
    @test isa(get_loc(loc1,azi,dist),GeoLoc)
    @test !isempty(get_loc(loc1,azi,dist))
    @test isempty(get_loc(GeoLoc(),azi,dist))
    loc3 = get_loc(loc1,get_azi(loc1,loc2),get_dist(loc1,loc2))
    @test isapprox(loc2.lat,loc3.lat)
    @test isapprox(loc2.lon,loc3.lon)

    # test get_loc with CorrData
    C = CorrData()

    # test empty location
    @test isa(get_loc(C)[2],GeoLoc)
    @test isempty(get_loc(C)[2])

    # test with loc1
    C.loc = loc1
    C.azi = get_azi(loc1,loc2)
    C.dist = get_dist(loc1,loc2)
    loc4 = get_loc(C)[2]
    @test isapprox(loc2.lat,loc4.lat)
    @test isapprox(loc2.lon,loc4.lon)
end
