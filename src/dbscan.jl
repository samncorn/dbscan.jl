module dbscan

using NearestNeighbors
using Logging
using Dates

function DBSCAN(points, r, min_pts; leafsize = 25, reorder = true, n_threads = 1, metric = Euclidean(), max_pts = Inf)
    t0 = now()
    tree = KDTree(points, metric; leafsize = leafsize, reorder = reorder)
    tf = now()
    @debug("$(canonicalize(tf - t0)) to build kdtree")

    N = length(points)
    labels = zeros(Int, N)

    # no point in chunking, override
    n_threads = N > 100n_threads ? n_threads : 1 

    mergers = Vector{Vector{Tuple{Int, Int}}}(undef, n_threads)
    chunksize = ceil(Int, N / n_threads)
    chunks = collect(Iterators.partition(eachindex(points), chunksize))
    @debug "performing range searches"
    Threads.@threads for i_thread in 1:n_threads
        @debug "starting thread $(i_thread)"
        t0 = now()
        idxs = chunks[i_thread]
        mergers[i_thread] = _dbscan_kernel!(labels, tree, points, idxs, r, min_pts, max_pts; reordered = reorder)
        tf = now()
        @debug("$(canonicalize(tf - t0)) to process thread $(i_thread)")
    end 

    locks = [ReentrantLock() for _ in labels]
    for (m_i, m) in enumerate(mergers)
        t0 = now()
        # Threads.@threads for i_thread in 1:n_threads
        #     for k in i_thread:n_threads:length(m)
        #         (i, j) = m[k]
        #         join_labels_locking!(labels, locks, i, j)
        #     end
        # end
        # non parallel version
        for (i, j) in m
            join_labels!(labels, i, j)
        end
        tf = now()
        @debug("$(canonicalize(tf - t0)) to merge clusters $(m_i)")
    end

    @debug "updating labels to cluster roots"
    t0 = now()
    promote_labels!(labels)
    tf = now()
    @debug("$(canonicalize(tf - t0)) to update labels")
    # collect labels into clusters

    @debug "DBSCAN completed"
    return collect_labels(labels)
end

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


# """ New and improved dbscan

# the partition arg sends point to user chunks.
# """
# function DBSCAN2(points, r, min_pts; n_threads = 1, metric = Euclidean(), max_pts = Inf)

# end

""" neighborhood must return an iterable of tuples of (point, bool), where the bool indicates whether the point requires follow up
to resolve a boundary crossing. The neighborhood return iterable must also have a length method defined on it
"""
function _dbscan_kernel_2!(labels, mergers, points, neighborhood, r, min_pts, max_pts)
    @assert length(labels) == length(points)
    for (i, p_i) in enumerate(points)
        neighbors = neighborhood(p_i, r)
        (length(neighbors) <= min_pts || length(neighbors) > max_pts) && continue # need equality because the point will also be counted
        for (j, in_chunk) in neighbors
            if in_chunk
                join_labels!(labels, i, j)
            else
                push!(mergers, (i, j))
            end
        end
    end
end

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
    return values(clusters) # cluster labels not important past this
end

end
