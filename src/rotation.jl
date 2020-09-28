export rotate, rotate!

"""
    rotate!(corrs,azi,baz)

Rotate cross-correlations from E,N,Z to R,T,Z coordinate system.
# Arguments
- `corrs::Array{CorrData,1}`: Array of EE, EN, EZ, NE, NN, NZ, ZE, ZN, ZZ CorrData
          for a single station-pair.
- `azi::Real`: Azimuth (in radians) from station 1 to station 2
- `baz::Real`: Back azimuth (in radians) from station 1 to station 2
"""
function rotate!(corrs::Array{CorrData,1},azi::Real,baz::Real)
    T = eltype(corrs[1].corr)
    azi = T(azi)
    baz = T(baz)
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
    # update ZZ rotation
    corrs[comp_dict["ZZ"]].rotated = true

    return nothing
end
rotate(corrs::Array{CorrData,1},azi::Real,baz::Real) =
                      (U = deepcopy(corrs);rotate!(U,azi,baz);return U)

"""
  rotation_kernel(rot_matrix, to_rotate)

Rotation function for 2D correlation functions.

# Arguments
- `rot_matrix::Array`: Rotation matrix (either 2x2 or 4x4).
- `to_rotate::Array`: Array to be rotated column-wise.
"""
function rotation_kernel(rot_matrix::AbstractArray, to_rotate::AbstractArray)
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
function update_corr!(C::CorrData,new_corr::AbstractArray,comp::String)
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
