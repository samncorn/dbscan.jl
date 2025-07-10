module dbscan

using StaticArrays
using StatsBase
using NearestNeighbors
using Logging
using Dates
using Logging
using Printf

# convenience function so that we don't need linear algebra as a deopendecy
dot(x, y) = sum(x .* y)

# """ abstarcts out the query, partition, and metric functions
# """
# # function DBSCAN_kernel(points, radius, min_pts, query, partition, metric; n_threads = 1)
# function DBSCAN_kernel(points, radius, min_pts, query, partition, metric; n_threads = 1)
#     labels = zeros(length(points))

#     # for 
#     chunks = partition(points)
#     Threads.@threads for i_th in 1:n_threads
#         for i_c in i_th:n_threads:length(chunks)

#         end
#     end


# end

""" multithreaded fixed sized cell implementation

bins points into cell with width equal to the clustering radius. Then the neighborhood check is equivalent to checking 
the surrounding cells (if they exist, only cells with points are instantiated).
"""
function DBSCAN_cells(points::AbstractVector{SVector{D, T}}, radius, min_pts; n_threads = 1, max_pts = Inf) where {D, T}
    cells = Dict{SVector{D, Int32}, Vector{Int}}()
    center = mean(points)

    # assign points to cells
    for (i, p) in pairs(points)
        cell = SVector{D, Int32}(find_cell(p - center, radius))
        if haskey(cells, cell)
            push!(cells[cell], i)
        else
            cells[cell] = [i]
        end
    end

    # skipped overfilled cells
    n_deleted = 0
    for (cell, idx) in cells
        if length(idx) > max_pts
            delete!(cells, cell)
            n_deleted += 1
        end
    end
    @info "deleted $(n_deleted) cells with > $(max_pts) points"

    _neighbor_cells = map(x -> SVector{D}(x), Iterators.product(((0, 1) for _ in 1:D)...))
    # query cells
    labels = zeros(UInt32, length(points))

    # old chunking (almost random)
    # chunks = collect(Iterators.partition(keys(cells), floor(Int, length(cells) / n_threads)))

    # new chunking
    # use larger chunks to distribute the threads
    # chunk_width = chunk_scale*radius
    # chunks = Dict{SVector{D, Int32}, Vector{SVector{D, Int32}}}()
    # for cell in keys(cells)
    #     # with chunk side length, compute a broader key 
    #     chunk = find_cell(cell, chunk_width)
    #     if haskey(chunks, chunk)
    #         push!(chunks[chunk], cell)
    #     else
    #         chunks[chunk] = [cell]
    #     end
    # end
    

    # new new chunking
    # based on number of threads, choose level to partition "quad"rants
    # n_quads  = 2^D
    # n_chunks = 1
    # while n_chunks <= n_threads
    #     n_chunks *= n_quads
    # end
    # @info @sprintf "divided cells into %i chunks of width %.3e" length(chunks) chunk_width

    # new new new chunking
    depth = ceil(Int, log2(length(points)))
    leafsize = 2^depth
    tree = KDTree(points; leafsize = leafsize, reorder = false)

    chunks = [NearestNeighbors.get_leaf_range(tree.tree_data, i+tree.tree_data.n_internal_nodes) for i in 1:tree.tree_data.n_leafs]
    merges = [Tuple{Int, Int}[] for _ in chunks]
    @info "using $(2^depth) chunks"
    for chunk in chunks
        @info "chunk $(chunk)"
    end


    chunk_keys = collect(keys(chunks))
    # cell_chunk = Dict{SVector{D, Int32}, Int}()
    # for (i, chunk) in enumerate(chunk_keys)
    #     for cell in chunks[chunk]
    #         cell_chunk[cell] = i
    #     end
    # end

    Threads.@threads for i_th in 1:n_threads
        for i_c in i_th:n_threads:length(chunks)
            chunk = chunks[chunk_keys[i_c]]
            merge = merges[i_c]

            for i_tree in chunk
                i = tree.indices[i_tree]
                p_i = points[i]
                celli = find_cell(p_i, radius)
                n_core = 0
                for _n in _neighbor_cells
                    cellj = celli + _n

                    if !haskey(cells, cellj)
                        continue
                    end

                    for j in cells[cellj]
                        p_j = points[j]
                        if dot(p_i - p_j, p_i - p_j) < radius^2
                            n_core += 1
                        end
                    end
                end
                
                if n_core >= min_pts # is core point
                    for _n in _neighbor_cells
                        cellj = celli + _n

                        if !haskey(cells, cellj)
                            continue
                        end

                        for j in cells[cellj]
                            p_j = points[j]
                            if dot(p_i - p_j, p_i - p_j) < radius^2
                                if in(j, chunk)
                                    join_labels!(labels, i, j)
                                else
                                    push!(merge, (i, j))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    @info @sprintf "merge list has %i pairs to check (%.3e Bytes)" sum(length.(merges)) sum(sizeof.(merges))

    for (i, j) in Iterators.flatten(merges)
        join_labels!(labels, i, j)
    end

    empty!.(merges)
    empty!(merges)

    promote_labels!(labels)

    return labels
end

# function count_neighbors(p::SVector{D, T}, radius, cells) where {D, T}
#     cell = find_cell(p, radius)
#     # for i in 2^length(p)
#     #     n = zero(MVector{D, Int})
#     #     for j in 1:length(p)

#     #     end
#     # end
#     for _neighbor in Iterators.product(((0, 1) for _ in 1:d)...)
#         cellj = SVector{D, Int}(_neighbor) + cell
# end

# function kernel_cells!(labels, merge, p::SVector{D, T}, cells, radius) where {T}

# end

function find_cell(point::SVector{D, T}, radius) where {D, T}
    cell = @MVector zeros(Int, D)
    for d in eachindex(point)
        cell[d] = sign(point[d]) * floor(Int, abs(point[d] / radius)) 
    end
    return SVector{D}(cell)
    # return join(cell, ',')
end

function dbscan_tree(points, radius, min_pts; leaf_size = 100, metric = Euclidean(), n_threads = 1)
    # construct the tree
    tree = KDTree(points, metric; leafsize = leaf_size)

    # construct the partitioning bounds

    # loop through points building clusters
    labels = zeros(UInt, length(points))

end

# we need a custom function for rolling through points
# function tree_search(point, radius, tree)
#     # search the nodes for potential overlap
#     # climb the tree from left to right, ascending nodes the intersect the search box
#     node = 1
#     while node <= # condition for being able to ascend one of the branch nodes

#     end
# end

# struct TreeNode{V}
#     index::Int
#     range::UnitRange{Int}
#     bounds::NearestNeighbors.HyperRectangle{V}
# end

""" apply a function f to all points in a tree with radius of point
f must be a single function of the point index (any other information has to be passed by closure)

hacked together code in the main nearest nerighbors branch
"""
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

"""
single threaded test
"""
function DBSCAN_tree(points, radius, min_pts; leafsize = 25, reorder = true, n_threads = 1, metric = Euclidean(), max_pts = Inf)
    tree = KDTree(points, metric; leafsize = leafsize, reorder = reorder)
    
    labels = zeros(UInt, length(points))
    for (i, p) in enumerate(points)
        n = inrangecount(tree, p, radius)
        if n >= min_pts
            inrange(j -> join_labels!(labels, i, j), tree, p, radius)
        end
    end

    return labels
end

# function DBSCAN_tree(points, r, min_pts; leafsize = 25, reorder = true, n_threads = 1, metric = Euclidean(), max_pts = Inf)
#     t0 = now()
#     tree = KDTree(points, metric; leafsize = leafsize, reorder = reorder)
#     tf = now()
#     @debug("$(canonicalize(tf - t0)) to build kdtree")

#     N = length(points)
#     labels = zeros(Int, N)

#     # no point in chunking, override
#     n_threads = N > 100n_threads ? n_threads : 1 

#     mergers = Vector{Vector{Tuple{Int, Int}}}(undef, n_threads)
#     chunksize = ceil(Int, N / n_threads)
#     chunks = collect(Iterators.partition(eachindex(points), chunksize))
#     @debug "performing range searches"
#     Threads.@threads for i_thread in 1:n_threads
#         @debug "starting thread $(i_thread)"
#         t0 = now()
#         idxs = chunks[i_thread]
#         mergers[i_thread] = _dbscan_kernel!(labels, tree, points, idxs, r, min_pts, max_pts; reordered = reorder)
#         tf = now()
#         @debug("$(canonicalize(tf - t0)) to process thread $(i_thread)")
#     end 

#     # locks = [ReentrantLock() for _ in labels]
#     for (m_i, m) in enumerate(mergers)
#         t0 = now()
#         # Threads.@threads for i_thread in 1:n_threads
#         #     for k in i_thread:n_threads:length(m)
#         #         (i, j) = m[k]
#         #         join_labels_locking!(labels, locks, i, j)
#         #     end
#         # end
#         # non parallel version
#         for (i, j) in m
#             join_labels!(labels, i, j)
#         end
#         tf = now()
#         @debug("$(canonicalize(tf - t0)) to merge clusters $(m_i)")
#     end

#     @debug "updating labels to cluster roots"
#     t0 = now()
#     promote_labels!(labels)
#     tf = now()
#     @debug("$(canonicalize(tf - t0)) to update labels")
#     # collect labels into clusters

#     @debug "DBSCAN completed"
#     return collect_labels(labels)
# end

function _dbscan_kernel!(labels, tree, points, points_idx, r, min_pts, max_pts; reordered = false)
    # if reordered, do the operations on the trees internal ordering, for better locality
    # still have to map out the 
    mergers = Tuple{Int, Int}[] # save edges that cross a chunk boundary for post processing
    neighborhood = Int[]
    for i in points_idx
        p_i = points[i]

        empty!(neighborhood)
        inrange!(neighborhood, tree, p_i, r)
        (length(neighborhood) <= min_pts || length(neighborhood) > max_pts) && continue # need equality because the point will also be counted
        # length(neighborhood) > max_pts && throw(DomainError("dbscan exceeded maximum number of points in a neighborhood"))

        if labels[i] == 0 # unvisited, flag as core point
            labels[i] = i
        end

        # ci = find_root(i, labels)
        for j in neighborhood
            i == j && continue
            if in(j, points_idx)
                join_labels!(labels, i, j)
                # if labels[j] == 0
                #     labels[j] == ci
                # else
                #     cj = find_root(j, labels)
                #     if ci > cj
                #         labels[cj] = ci
                #         labels[j]  = ci
                #     elseif ci < cj
                #         labels[ci] = cj
                #         labels[i]  = cj 
                #     end
                # end
            else
                push!(mergers, (i, j))
            end
        end
    end
    return mergers
end

# """
# labels
#     total point -> cluster label array
# mergers
#     records whether to check points in the clean up sweep
# range_query
#     queries the region around a point
# chunk
#     points to run clustering on
# chunk map
#     takes a points label, return whether it is in the chunk

# """
# function _dbscan_kernel_2!(labels, mergers, range_query, chunk, chunk_map, r, min_pts, max_pts)
#     neighborhood = Int[]

# end


# """ New and improved dbscan

# the partition arg sends point to user chunks.
# """
# function DBSCAN2(points, r, min_pts; n_threads = 1, metric = Euclidean(), max_pts = Inf)

# end

# """ neighborhood must return an iterable of tuples of (point, bool), where the bool indicates whether the point requires follow up
# to resolve a boundary crossing. The neighborhood return iterable must also have a length method defined on it
# """
# function _dbscan_kernel_2!(labels, mergers, points, neighborhood, r, min_pts, max_pts)
#     for (i, p_i) in enumerate(points)
#         neighbors = neighborhood(p_i, r)
#         (length(neighbors) <= min_pts || length(neighbors) > max_pts) && continue # need equality because the point will also be counted
#         for (j, in_chunk) in neighbors
#             if in_chunk
#                 join_labels!(labels, i, j)
#             else
#                 push!(mergers, (i, j))
#             end
#         end
#     end
# end

"""
loop through labels, changing the label to match the root label
should add splicing w/ some locks
""" 
function promote_labels!(labels)
    Threads.@threads for i in eachindex(labels)
        # find cluster root
        if labels[i] != 0 
            j = find_root(i, labels)
            labels[i] = j
        end
    end
end

function find_root(i, labels)
    j = i
    n = 0
    max_n = length(labels)
    while j != labels[j]
        j = labels[j]
        n += 1
        j == 0 && throw("ERROR: noise point root of cluster")
        n > max_n && throw("ERROR: probably stuck in a loop")
    end
    return j
end

""" appropriatley labels points. 

Uses Rem's algorithm w/ splicing from 
M. Patwary, J. Blair, and F. Manne, “Experiments on union-find al-
gorithms for the disjoint-set data structure,” in Proceedings of the
9th International Symposium on Experimental Algorithms (SEA 2010).
Springer, LNCS 6049, 2010, pp. 411–423.
"""
function join_labels!(labels, ii, jj)
    # handle noise points
    if labels[ii] == labels[jj] == 0
        l = max(ii, jj)
        labels[ii] = l
        labels[jj] = l
    elseif labels[jj] == 0
        labels[jj] = jj
    elseif labels[ii] == 0
        labels[ii] = ii
    end

    # splicing to flatten the tree
    i = ii
    j = jj
    while labels[i] != labels[j]
        if labels[i] < labels[j]
            if labels[i] == i
                labels[i] = labels[j]
            else
                k = labels[i]
                labels[i] = labels[j]
                i = k
            end
        else
            if labels[j] == j
                labels[j] = labels[i]
            else
                k = labels[j]
                labels[j] = labels[i]
                j = k
            end
        end
    end
    # update the original points
    l = labels[i]
    labels[ii] = l
    labels[jj] = l
end

function join_labels_locking!(labels, locks, ii, jj)
    i = ii
    j = jj
    # check for noise points
    if labels[i] == labels[j] == 0
        lock(locks[i]) do 
            lock(locks[j]) do
                if labels[i] == labels[j] == 0 
                    l = max(i, j)
                    labels[i] == l
                    labels[j] == l
                end
            end
        end
    elseif labels[i] == 0
        lock(locks[i]) do
            if labels[i] == 0
                labels[i] = labels[j]
            end
        end
    elseif labels[j] == 0
        lock(locks[j]) do
            if labels[j] == 0
                labels[j] = labels[i]
            end
        end
    end

    # splicing
    while labels[i] != labels[j]
        if labels[i] < labels[j]
            if labels[i] == i
                lock(locks[i]) do 
                    if labels[i] == i # double check in case it's been changed
                        labels[i] = labels[j]
                    end
                end
            end
            i = labels[i]
        else
            if labels[j] == j
                lock(locks[j]) do
                    if labels[j] == j  # double check in case it's been changed
                        labels[j] = labels[i]
                    end
                end
            end
            j = labels[j]
        end
    end
    l = labels[i]
    lock(locks[ii]) do
        labels[ii] = l
    end
    lock(locks[jj]) do
        labels[jj] = l
    end
end

function collect_labels(labels::Vector{I}) where I
    clusters = Dict{I, Vector{Int}}()
    for (i, l) in pairs(labels)
        if l == 0
            continue
        elseif haskey(clusters, l)
            push!(clusters[l], i)
        else
            clusters[l] = [i]
        end
    end
    return collect(values(clusters)) # cluster labels not important past this
end

end
