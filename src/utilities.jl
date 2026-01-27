# This is useful for others.
export crs3

"""
Returns the skew-symmetric "cross product matrix" for a given vector, such that
`crs3(v) * x == cross(v, x)`.
"""
function crs3(r)
    @assert length(r) == 3 "The length of the input to `crs3` was $(length(r)) but should be 3."
    return @SMatrix [0 -r[3] r[2]; r[3] 0 -r[1]; -r[2] r[1] 0]
end
