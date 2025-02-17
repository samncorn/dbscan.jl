module dbscan

using NearestNeighbors

function DBSCAN(points, r, min_pts; leafsize = 25, reorder = true, n_chunks = 1)
    tree = KDTree(points; leafsize = leafsize, reorder = reorder)

    N = length(points)
    labels = zeros(Int, N)

    mergers = Vector{Vector{Tuple{Int, Int}}}(undef, n_chunks)
    chunksize = ceil(Int, N / n_chunks)
    chunks = collect(Iterators.partition(eachindex(points), chunksize))

    Threads.@threads for i in 1:n_chunks
        points_idx = chunks[i]
        mergers[i] = _dbscan_kernel!(labels, tree, points, points_idx, r, min_pts, Inf)
    end
    resolve_cluster_merge!(labels, collect(Iterators.flatten(mergers)))
    promote_labels!(labels)
    return labels
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

function resolve_cluster_merge!(labels, mergers)
    locks = [ReentrantLock() for _ in labels]
    Threads.@threads for (i, j) in mergers
        begin 
            lock(locks[i])
            lock(locks[j])
            try
                if labels[i] == 0 && labels[j] == 0 
                    # only join two noise points if min pts is small enough
                    if min_pts == 2
                        labels[i] = i
                        labels[j] = i
                    end
                else
                    join_labels!(labels, i, j)
                end
            finally
                unlock(locks[i])
                unlock(locks[j])
            end
        end
    end
end

end
