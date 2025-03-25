module dbscan

using NearestNeighbors
using Logging

function DBSCAN(points, r, min_pts; leafsize = 25, reorder = true, n_threads = 1, metric = Euclidean(), max_pts = Inf)
    tree = KDTree(points, metric; leafsize = leafsize, reorder = reorder)

    N = length(points)
    labels = zeros(Int, N)

    mergers = Vector{Vector{Tuple{Int, Int}}}(undef, n_threads)
    chunksize = ceil(Int, N / n_threads)
    chunks = collect(Iterators.partition(eachindex(points), chunksize))
    @debug "performing range searches"
    Threads.@threads for i_thread in 1:n_threads
        @debug "starting thread $(i_thread)"
        idxs = chunks[i_thread]
        mergers[i_thread] = _dbscan_kernel!(labels, tree, points, idxs, r, min_pts, max_pts)
        @debug "thread $(i_thread) completed"
    end 

    locks = [ReentrantLock() for _ in labels]
    for (m_i, m) in enumerate(mergers)
        Threads.@threads for (i, j) in m
            join_labels_locking!(labels, locks, i, j)
        end
        @debug "merged $(m_i) edge clusters"
    end

    # return labels

    @debug "updating labels to cluster roots"
    promote_labels!(labels)
    # collect labels into clusters

    @debug "DBSCAN completed"
    return collect_labels(labels)
end

function _dbscan_kernel!(labels, tree, points, points_idx, r, min_pts, max_pts)
    mergers = Tuple{Int, Int}[] # save edges that cross a chunk boundary for post processing
    neighborhood = Int[]
    for i in points_idx
        p_i = points[i]

        empty!(neighborhood)
        inrange!(neighborhood, tree, p_i, r)
        length(neighborhood) <= min_pts && continue # need equality because the point will also be counted
        length(neighborhood) > max_pts && throw(DomainError("dbscan exceeded maximum number of points in a neighborhood"))

        if labels[i] == 0 # unvisited, flag as core point
            labels[i] = i
        end

        ci = find_root(i, labels)
        for j in neighborhood
            if in(j, points_idx)
                # join_labels!(labels, i, j)
                if labels[j] == 0
                    labels[j] == ci
                else
                    cj = find_root(j, labels)
                    if ci > cj
                        labels[cj] = ci
                        labels[j]  = ci
                    elseif ci < cj
                        labels[ci] = cj
                        labels[i]  = cj 
                    end
                end
            else
                push!(mergers, (i, j))
            end
        end
    end
    return mergers
end

"""
loop through labels, changing the label to match the root label

not currently multihteraded, but could be using some locks
""" 
function promote_labels!(labels)
    for i in eachindex(labels)
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
    i = ii
    j = jj
    # splicing to flatten the tree
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
                labels[k] = labels[i]
                j = k
            end
        end
    end
    # update the original points
    l = max(labels[i], labels[j])
    labels[ii] = l
    labels[jj] = l
end

function join_labels_locking!(labels, locks, ii, jj)
    i = ii
    j = jj

    if labels[i] == 0
        lock(locks[i]) do
            labels[i] = labels[j]
        end
    end

    if labels[j] == 0
        lock(locks[j]) do
            labels[j] = labels[i]
        end
    end

    while labels[i] != labels[j]
        if labels[i] < labels[j]
            if labels[i] == i
                lock(locks[i]) do 
                    labels[i] = labels[j]
                end
            end
            i = labels[i]
        else
            if labels[j] == j
                lock(locks[j]) do
                    labels[j] = labels[i]
                end
            end
            j = labels[j]
        end
    end
    l = max(labels[i], labels[j])
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
