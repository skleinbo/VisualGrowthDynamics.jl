function project(P::Plane, v::Point3f)
    Point2f(dot(P.u, v), dot(P.v, v))
end

function validate_css(T, s::AbstractString, l=0)
    v = split(s, ',')
    if (length(v)==l || l==0) && all(x->!isnothing(tryparse(T, x)), v)
        return true
    end
    return false
end

xyz_string_to_point3(s::AbstractString) = Point3(tryparse.(Float64, split(s, ',')))

cs_string_to_vec(T, s::AbstractString) = tryparse.(T, split(s, ','))
