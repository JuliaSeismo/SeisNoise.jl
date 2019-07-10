export rotate, rotate_corrs, rotate_corrs!, get_corr

"""
  rotate(corrpaths, ROTOUT)

Rotate all correlations from ENZ to RTZ in directory `CORROUT`.
# Arguments
- `corrpaths::Array{String,1}`: Array of filepaths for 9 component correlations.
- `ROTOUT::String`: Output directory for rotated correlations.
- `Overwrite::Bool`: Overwrite non-rotated correlations with rotated correlations.
"""
function rotate(corrpaths::Array{String,1}, ROTOUT::String; overwrite::Bool=false)

    # read correlations
    N = length(corrpaths)
    comps = Array{String}(undef,N)
    for ii = 1:N
        path, file = splitdir(corrpaths[ii])
        file = replace(file,".jld2"=>"")
        net1,sta1,chan1,net2,sta2,chan2 = split(file,'.')
        comps[ii] = chan1[end] * chan2[end]
    end

    # check if CORROUT exists and make new directory
    ROTOUT = expanduser(ROTOUT)
    isdir(ROTOUT) || mkpath(ROTOUT)

    # load data
    corrs = [load_corr(corrpaths[ii],comps[ii]) for ii = 1:N]
    azi,baz = corrs[1].azi, corrs[1].baz

    # do not rotate if no azi and baz in corrs
    if (azi == 0.) & (baz == 0.)
        println("No location data for $(corrs[1].name) $(corrs[1].id)")
    else
        # rotate and save correlations
        rotate_corrs!(corrs,azi,baz)
        for corr in corrs
            save_corr(corr,ROTOUT)
        end

        # delete non-rotated correlations
        if overwrite
            for ii = 1:N
                rm(corrpaths[ii])
            end
        end
    end

end


"""
    rotate_corrs(corrs,azi,baz)

Rotate cross-correlations from E,N,Z to R,T,Z coordinate system.
# Arguments
- `corrs::Array{CorrData,1}`: Array of EE, EN, EZ, NE, NN, NZ, ZE, ZN, ZZ CorrData
- `azi::Float64`: Azimuth (in radians) from station 1 to station 2
- `baz::Float64`: Back azimuth (in radians) from station 1 to station 2
"""
function rotate_corrs!(corrs::Array{CorrData,1},azi::Float64,baz::Float64)
    N = length(corrs)
    # get time intersect of all components
    ts = [c.t for c in corrs]
    inter = intersect(ts)[1]

    # subset so all corrs match in time
    for ii = 1:N
        ind = findall(x -> x ∈ inter, corrs[ii].t)
        corrs[ii].t = corrs[ii].t[ind]
        corrs[ii].corr = corrs[ii].corr[:,ind]
    end

    comp_dict = Dict()
    for ii = 1:N
        comp_dict[corrs[ii].comp] = ii
    end

    # get EE EN NN NE comps for rotation to TT RR TR RT
    # see Lin et al., 2008: Surface wave tomography of the western
    # United States from ambient seismic noise, GJI
    big_rot = [[-cos(azi)*cos(baz) cos(azi)*sin(baz) -sin(azi)*sin(baz) sin(azi)*cos(baz)];
               [-sin(azi)*sin(baz) -sin(azi)*cos(baz) -cos(azi)*cos(baz) -cos(azi)*sin(baz)];
               [-cos(azi)*sin(baz) -cos(azi)*cos(baz) sin(azi)*cos(baz) sin(azi)*sin(baz)];
               [-sin(azi)*cos(baz) sin(azi)*sin(baz) cos(azi)*sin(baz) -cos(azi)*cos(baz)]]
    to_rot = Array{eltype(corrs[1].corr)}(undef,4,size(corrs[1].corr)[1],length(inter))
    to_rot[1,:,:] = corrs[comp_dict["EE"]].corr
    to_rot[2,:,:] = corrs[comp_dict["EN"]].corr
    to_rot[3,:,:] = corrs[comp_dict["NN"]].corr
    to_rot[4,:,:] = corrs[comp_dict["NE"]].corr
    rotated = rotation_kernel(big_rot,to_rot)

    # update correlations from EE EN NN NE to TT RR TR RT
    update_corr!(corrs[comp_dict["EE"]],rotated[1,:,:],"TT")
    update_corr!(corrs[comp_dict["EN"]],rotated[2,:,:],"RR")
    update_corr!(corrs[comp_dict["NN"]],rotated[3,:,:],"TR")
    update_corr!(corrs[comp_dict["NE"]],rotated[3,:,:],"RT")

    # rotate from EZ NZ to RZ TZ
    azi_rot = [[cos(azi) -sin(azi)];
               [sin(azi) cos(azi)]]
    to_rot = Array{eltype(corrs[1].corr)}(undef,2,size(corrs[1].corr)[1],length(inter))
    to_rot[1,:,:] = corrs[comp_dict["NZ"]].corr
    to_rot[2,:,:] = corrs[comp_dict["EZ"]].corr
    rotated = rotation_kernel(azi_rot,to_rot)
    update_corr!(corrs[comp_dict["NZ"]],rotated[1,:,:],"RZ")
    update_corr!(corrs[comp_dict["EZ"]],rotated[2,:,:],"TZ")

    # rotate from ZE ZN to ZR ZT
    baz_rot = [[-cos(baz) -sin(baz)];
               [sin(baz) -cos(baz)]]
    to_rot = Array{eltype(corrs[1].corr)}(undef,2,size(corrs[1].corr)[1],length(inter))
    to_rot[1,:,:] = corrs[comp_dict["ZN"]].corr
    to_rot[2,:,:] = corrs[comp_dict["ZE"]].corr
    rotated = rotation_kernel(baz_rot,to_rot)
    update_corr!(corrs[comp_dict["ZN"]],rotated[1,:,:],"ZR")
    update_corr!(corrs[comp_dict["ZE"]],rotated[2,:,:],"ZT")

    return corrs
end
rotate_corrs(corrs::Array{CorrData,1},azi::Float64,baz::Float64) =
                      (U = deepcopy(corrs);rotate_corrs!(U,azi,baz);return U)

"""
  rotation_kernel(rot_matrix, to_rotate)

Rotation function for 2D correlation functions.

# Arguments
- `rot_matrix::Array`: Rotation matrix (either 2x2 or 4x4).
- `to_rotate::Array`: Array to be rotated column-wise.
"""
function rotation_kernel(rot_matrix::Array, to_rotate::Array)
    rotated = similar(to_rotate)
    for ii = 1:size(to_rotate,2)
        for jj = 1:size(to_rotate,3)
            rotated[:,ii,jj] .= rot_matrix * to_rotate[:,ii,jj]
        end
    end
    return rotated
end

"""
  update_corr!(C,comp)

Change component, name and corr of rotated correlation `C` to component `comp`.

# Arguments
- `C::CorrData`: CorrData object.
- `new_corr::Array`: Rotated correlation matrix.
- `comp::String`: New component after rotation, e.g. "RR"
"""
function update_corr!(C::CorrData,new_corr::Array,comp::String)
    # change comp of correlation
    C.comp = comp

    # change name of correlation
    newname = split(C.name,'.')
    newname[4] = newname[4][1:end-1] * comp[1]
    newname[8] = newname[8][1:end-1] * comp[2]
    C.name = join(newname,'.')

    # change rotated flag
    C.rotated = true

    # update correlation
    C.corr = new_corr
    return C
end

"""
  get_corr(CORROUT)

Get name of all correlations in directory `CORROUT` for rotation.

# Arguments
- `CORROUT::String`: Path to directory of correaltion files, e.g. /home/CORR/
"""
function get_corr(CORROUT::String; autocorr::Bool=false)
    files = glob("*",CORROUT)
    corrs = [splitdir(f)[2] for f in files]
    corrs = [replace(c,".jld2"=>"") for c in corrs]

    # get net1.sta1.net2.sta2, comp and if autocorrelation from name of
    # correlation jld2 file
    N = length(corrs)
    NETSTA = Array{String,1}(undef,N)
    COMP = Array{String,1}(undef,N)
    AUTO = Array{Bool,1}(undef,N)

    for ii = 1:N
        net1,sta1,chan1,net2,sta2,chan2 = split(corrs[ii],'.')
        NETSTA[ii] = join([net1,sta1,net2,sta2],'.')
        COMP[ii] = chan1[end] * chan2[end]
        AUTO[ii] = net1*sta1 == net2*sta2
    end

    # create dataframe
    df = DataFrame(NETSTA = NETSTA, COMP = COMP, AUTO=AUTO, FILE = files)

    # remove autocorrelations
    df = df[df[:AUTO] .== autocorr,:]

    return df2files(df)
end

"""
  df2files(df)

Helper function to group correlation pairs.

Returns an array of file paths for each correlation pair in dataframe.

# Arguments
- `df::DataFrame`: DataFrame with correlation pair name, filepath and component.
"""
function df2files(df::DataFrame)
    grouped = groupby(df,:NETSTA)
    N = length(grouped)
    out = Array{Array{String,1},1}(undef,N)
    for ii = 1:N
        out[ii] = collect(grouped[ii][:FILE])
    end
    return out
end
