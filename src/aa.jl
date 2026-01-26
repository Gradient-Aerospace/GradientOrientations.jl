export AxisAngle, AA, AA_F64
export aax, aay, aaz

"TODO"
@kwdef struct AxisAngle{T} <: AbstractOrientation{T}
    axis::SVector{3, T}
    angle::T
end
const AA = AxisAngle
const AA_F64 = AxisAngle{Float64}

"TODO"
@kwdef struct AxisAngleDeg{T} <: AbstractOrientation{T}
    axis::SVector{3, T}
    angle::T
end
const AADeg = AxisAngleDeg
const AADeg_F64 = AxisAngleDeg{Float64}

Base.rad2deg(aa::AxisAngle) = AADeg(aa.axis, rad2deg(aa.angle))
Base.deg2rad(aa_deg::AxisAngleDeg) = AA(aa_deg.axis, deg2rad(aa_deg.angle))
Base.convert(::Type{AxisAngle}, aa_deg::AxisAngleDeg) = deg2rad(aa_deg)
Base.convert(::Type{AxisAngle{T}}, aa_deg::AxisAngleDeg{T}) where {T} = deg2rad(aa_deg)
Base.convert(::Type{AxisAngleDeg}, aa::AxisAngle) = rad2deg(aa)
Base.convert(::Type{AxisAngleDeg{T}}, aa::AxisAngle{T}) where {T} = rad2deg(aa)

################
# Constructors #
################

"TODO"
aax(θ::T) where {T} = AA{T}(SA[one(T), zero(T), zero(T)], θ)

"TODO"
aay(θ::T) where {T} = AA{T}(SA[zero(T), one(T), zero(T)], θ)

"TODO"
aaz(θ::T) where {T} = AA{T}(SA[zero(T), zero(T), one(T)], θ)

Base.zero(aa::AA) = zero(typeof(aa))
Base.zero(::Type{<:AA}) = zero(AA_F64)
Base.zero(::Type{AA{T}}) where {T} = AA{T}(SA[one(T), zero(T), zero(T)], zero(T))

"Constructs a random AxisAngle orientation following a uniform distribution on SO(3)."
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{AA{T}}) where {T}
    return erp2aa(rand(rng, ERP{T}))
end

##############
# Operations #
##############

Base.inv(aa::AA) = typeof(aa)(-aa.axis, aa.angle)
compose(a::AA, b::AA) = erp2aa(compose(aa2erp(a), aa2erp(b))) # TODO: Be more direct.
reframe(aa::AA, v) = reframe(aa2erp(aa),  v) # TODO: Be more direct.
difference(a::AA, b::AA) = erp2aa(difference(aa2erp(a), aa2erp(b))) # TODO: Be more direct.
distance(aa::AA) = abs(aa.angle)
interpolate(a::AA, b::AA, f) = erp2aa(interpolate(aa2erp(a), aa2erp(b), f)) # TODO: Be more direct.
# TODO: rate?

#############
# Iteration #
#############

# Helpful functions for iterating over AxisAngle.
Base.length(::AA) = 2
Base.eltype(::AA{T}) where {T} = Union{T, SVector{3, T}}
Base.size(::AA) = (2,)

# Allow a user to iterate over the elements of AA, e.g. for splatting.
function Base.iterate(aa::AA, state = 1)
    state == 1 && return (aa.axis, state + 1)
    state == 2 && return (aa.angle, state + 1)
    return nothing
end

# Provide linear indexing behavior.
function Base.getindex(aa::AA, k)
    k == 1 && return aa.axis
    k == 2 && return aa.angle
    throw(BoundsError(aa, k))
end
Base.firstindex(aa::AA) = 1
Base.lastindex(aa::AA) = 2

# Let's tell the user they can't do this.
function Base.setindex!(aa::AA, k, value)
    error("AxisAngle is immutable and cannot support setindex!.")
end

#################
# Miscellaneous #
#################

function Base.show(io::IO, aa::AxisAngle)
    print(io, "AxisAngle(axis = $(aa.axis), angle = $(aa.angle))")
end

########################################
# Conversions to Non-Orientation Types #
########################################

# What would be a normal type to convert to? A tuple?
