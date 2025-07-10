using NearestNeighbors
using StaticArrays
using BenchmarkTools

function inrange(
    f::F, 
    tree::KDTree,
    point::AbstractVector,
    radius::Number,
    # idx_in_ball::Union{Nothing, Vector{<:Integer}} = Int[]
    ) where {F}
    init_min = NearestNeighbors.get_min_distance_no_end(tree.metric, tree.hyper_rec, point)
    return inrange_kernel(f, tree, 1, point, NearestNeighbors.eval_pow(tree.metric, radius), tree.hyper_rec, init_min)
end

function inrange_kernel(
    f::F, 
    tree::KDTree,
    index::Int,
    point::AbstractVector,
    r::Number,
    # idx_in_ball::Union{Nothing, Vector{<:Integer}},
    hyper_rec::NearestNeighbors.HyperRectangle,
    min_dist
    ) where {F}
    # Point is outside hyper rectangle, skip the whole sub tree
    if min_dist > r
        # return 0
        return nothing
    end

    # At a leaf node. Go through all points in node and add those in range
    if NearestNeighbors.isleaf(tree.tree_data.n_internal_nodes, index)
        return apply_inrange(f, tree, index, point, r)
    end

    split_val = tree.split_vals[index]
    split_dim = tree.split_dims[index]
    lo = hyper_rec.mins[split_dim]
    hi = hyper_rec.maxes[split_dim]
    p_dim = point[split_dim]
    split_diff = p_dim - split_val
    M = tree.metric

    # count = 0

    if split_diff > 0 # Point is to the right of the split value
        close = NearestNeighbors.getright(index)
        far = NearestNeighbors.getleft(index)
        hyper_rec_far = NearestNeighbors.HyperRectangle(hyper_rec.mins, @inbounds setindex(hyper_rec.maxes, split_val, split_dim))
        hyper_rec_close = NearestNeighbors.HyperRectangle(@inbounds(setindex(hyper_rec.mins, split_val, split_dim)), hyper_rec.maxes)
        ddiff = max(zero(p_dim - hi), p_dim - hi)
    else # Point is to the left of the split value
        close = NearestNeighbors.getleft(index)
        far = NearestNeighbors.getright(index)
        hyper_rec_far = NearestNeighbors.HyperRectangle(@inbounds(setindex(hyper_rec.mins, split_val, split_dim)), hyper_rec.maxes)
        hyper_rec_close = NearestNeighbors.HyperRectangle(hyper_rec.mins, @inbounds setindex(hyper_rec.maxes, split_val, split_dim))
        ddiff = max(zero(lo - p_dim), lo - p_dim)
    end
    # Call closer sub tree
    inrange_kernel(f, tree, close, point, r, hyper_rec_close, min_dist)

    # Call further sub tree with the new min distance
    split_diff_pow = NearestNeighbors.eval_pow(M, split_diff)
    ddiff_pow = NearestNeighbors.eval_pow(M, ddiff)
    diff_tot = NearestNeighbors.eval_diff(M, split_diff_pow, ddiff_pow, split_dim)
    new_min = NearestNeighbors.eval_reduce(M, min_dist, diff_tot)
    inrange_kernel(f, tree, far, point, r, hyper_rec_far, new_min)
    # return count
    return nothing
end

function apply_inrange(f::F, tree::NNTree, index::Int, point::AbstractVector, r::Number) where {F}
    # count = 0
    for z in NearestNeighbors.get_leaf_range(tree.tree_data, index)
        idx = tree.reordered ? z : tree.indices[z]
        if NearestNeighbors.check_in_range(tree.metric, tree.data[idx], point, r)
            # count += 1
            # idx_in_ball !== nothing && push!(idx_in_ball, idx)
            f(idx)
        end
    end
    # return count
end

points = rand(SVector{2, Float32}, 10_000)
tree = KDTree(points)

function test(labels, tree, p, r)
    inrange(i -> labels[i] = 1, tree, p, r)
end

p = @SVector zeros(Float32, 2)
r = 0.1

labels = zeros(length(points))
@btime test($labels, $tree, $p, $r)

function count_inrange(
    tree::KDTree,
    point::AbstractVector,
    radius::Number,
    # idx_in_ball::Union{Nothing, Vector{<:Integer}} = Int[]
    )
    init_min = NearestNeighbors.get_min_distance_no_end(tree.metric, tree.hyper_rec, point)
    return count_inrange_kernel(tree, 1, point, NearestNeighbors.eval_pow(tree.metric, radius), tree.hyper_rec, init_min)
end

function count_inrange_kernel(
    tree::KDTree,
    index::Int,
    point::AbstractVector,
    r::Number,
    # idx_in_ball::Union{Nothing, Vector{<:Integer}},
    hyper_rec::NearestNeighbors.HyperRectangle,
    min_dist
    )
    # Point is outside hyper rectangle, skip the whole sub tree
    if min_dist > r
        return 0
        # return nothing
    end

    # At a leaf node. Go through all points in node and add those in range
    if NearestNeighbors.isleaf(tree.tree_data.n_internal_nodes, index)
        return _count_inrange(tree, index, point, r)
    end

    split_val = tree.split_vals[index]
    split_dim = tree.split_dims[index]
    lo = hyper_rec.mins[split_dim]
    hi = hyper_rec.maxes[split_dim]
    p_dim = point[split_dim]
    split_diff = p_dim - split_val
    M = tree.metric

    count = 0

    if split_diff > 0 # Point is to the right of the split value
        close = NearestNeighbors.getright(index)
        far = NearestNeighbors.getleft(index)
        hyper_rec_far = NearestNeighbors.HyperRectangle(hyper_rec.mins, @inbounds setindex(hyper_rec.maxes, split_val, split_dim))
        hyper_rec_close = NearestNeighbors.HyperRectangle(@inbounds(setindex(hyper_rec.mins, split_val, split_dim)), hyper_rec.maxes)
        ddiff = max(zero(p_dim - hi), p_dim - hi)
    else # Point is to the left of the split value
        close = NearestNeighbors.getleft(index)
        far = NearestNeighbors.getright(index)
        hyper_rec_far = NearestNeighbors.HyperRectangle(@inbounds(setindex(hyper_rec.mins, split_val, split_dim)), hyper_rec.maxes)
        hyper_rec_close = NearestNeighbors.HyperRectangle(hyper_rec.mins, @inbounds setindex(hyper_rec.maxes, split_val, split_dim))
        ddiff = max(zero(lo - p_dim), lo - p_dim)
    end
    # Call closer sub tree
    count += count_inrange_kernel(tree, close, point, r, hyper_rec_close, min_dist)

    # Call further sub tree with the new min distance
    split_diff_pow = NearestNeighbors.eval_pow(M, split_diff)
    ddiff_pow = NearestNeighbors.eval_pow(M, ddiff)
    diff_tot = NearestNeighbors.eval_diff(M, split_diff_pow, ddiff_pow, split_dim)
    new_min = NearestNeighbors.eval_reduce(M, min_dist, diff_tot)
    count += count_inrange_kernel(tree, far, point, r, hyper_rec_far, new_min)
    # return count
    return count
end

function _count_inrange(tree::NNTree, index::Int, point::AbstractVector, r::Number)
    count = 0
    for z in NearestNeighbors.get_leaf_range(tree.tree_data, index)
        idx = tree.reordered ? z : tree.indices[z]
        if NearestNeighbors.check_in_range(tree.metric, tree.data[idx], point, r)
            count += 1
            # idx_in_ball !== nothing && push!(idx_in_ball, idx)
            # f(idx)
        end
    end
    return count
end


p = @SVector zeros(Float32, 2)
r = 0.1

@btime count_inrange($tree, $p, $r)

function test3(v, tree)
    p = SVector{2}((0.0, 0.0))
    empty!(v)
    inrange!(v, tree, p, 0.1)
end

v = Int[]
@btime test3($v, $tree)