module dbscan

using NearestNeighbors

function DBSCAN(points, r, min_pts; leafsize = 25, reorder = true, n_chunks = 1, metric = Euclidean())
    tree = KDTree(points, metric; leafsize = leafsize, reorder = reorder)

    N = length(points)
    labels = zeros(Int, N)

    mergers = Vector{Vector{Tuple{Int, Int}}}(undef, n_chunks)
    chunksize = ceil(Int, N / n_chunks)
    chunks = collect(Iterators.partition(eachindex(points), chunksize))

    Threads.@threads for i in 1:n_chunks
        points_idx = chunks[i]
        mergers[i] = _dbscan_kernel!(labels, tree, points, points_idx, r, min_pts, Inf)
    end 
    for m in mergers
        locks = [ReentrantLock() for _ in labels]
        Threads.@threads for (i, j) in m
            join_labels_locking!(labels, locks, i, j, min_pts)
        end
    end
    promote_labels!(labels)
    # collect labels into clusters
    return collect_labels(labels)
end

function _dbscan_kernel!(labels, tree, points, points_idx, r, min_pts, max_pts)
    mergers = Tuple{Int, Int}[] # save edges that cross a chunk boundary for post processing
    neighborhood = Int[]
    for i in points_idx
        p_i = points[i]

        empty!(neighborhood)
        inrange!(neighborhood, tree, p_i, r)
        length(neighborhood) < min_pts && continue
        length(neighborhood) > max_pts && throw("dbscan exceeded maximum number of points in a neighborhood")

        if labels[i] == 0 # unvisited, flag as core point
            labels[i] = i
        end

        for j in neighborhood
            i <= j && continue
            if !in(j, points_idx)
                push!(mergers, (i, j))
                continue # only assign
            end
            join_labels!(labels, i, j)
        end
    end
    return mergers
end

"""
loop through labels, changing the label to match the root label
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
    while j != labels[j]
        j = labels[j]
        j == j 
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
    i = ii
    j = jj
    while labels[i] != labels[j]
        if labels[i] < labels[j]
            if labels[i] == i
                labels[i] = labels[j]
                continue
            end
            k = i 
            labels[i] = labels[j]
            i = labels[k]
        else
            if labels[j] == j
                labels[j] = labels[i]
                continue
            end
            k = j
            labels[k] = labels[i]
            j = labels[k]
        end
    end
end

function join_labels_locking!(labels, locks, ii, jj, min_pts)
    i = ii
    j = jj

    if labels[i] == 0 && labels[j] == 0 && min_pts <= 2
        lock(locks[i])
        lock(locks[j])
        try
            k = max(i, j)
            labels[i] = k
            labels[j] = k
        finally
            unlock(locks[i])
            unlock(locks[j])
        end
        return nothing
    elseif labels[i] == 0
        labels[i] = j
    elseif labels[j] == 0
        labels[j] = i
    end

    while labels[i] != labels[j]
        if labels[i] < labels[j]
            if labels[i] == i
                lock(locks[i])
                try
                    labels[i] = labels[j]
                finally
                    unlock(locks[i])
                end
            end
            i = labels[i]
        else
            if labels[j] == j
                lock(locks[j])
                try
                    labels[j] = labels[i]
                finally
                    unlock(locks[j])
                end
            end
            j = labels[j]
        end
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
