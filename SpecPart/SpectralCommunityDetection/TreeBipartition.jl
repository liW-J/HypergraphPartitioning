mutable struct TreeCutToken
    nforced_0::Int
    nforced_1::Int
    nforced_01::Int
    total_vwts::Int
    exc0::Vector{Float64}
    exc1::Vector{Float64}
    areaPart::Vector{Int}
    cutCost0::Vector{Float64}
    cutCost1::Vector{Float64}
    ratioCost::Vector{Float64}
    cutCost::Vector{Float64}
    areaCost::Vector{Float64}
    tCost_v::Vector{Float64}
    polarity::Vector{Int}
    statusFlag::Vector{Int}
    areaUtil0::Vector{Int}
    areaUtil1::Vector{Int}
    pred::Vector{Int}
    hyperedges_flag::Vector{Int}

    # clique graph 切割成本 
    cliqueCutCost0::Vector{Float64}
    cliqueCutCost1::Vector{Float64}
    clique_weight::Float64
end

function TreeBipartition(tree::SimpleWeightedGraphs.SimpleGraph, tree_cut_token::TreeCutToken, hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vtxs::Pindex, capacities::Vector{Int}, ncuts::Int)
    cutInfo = analyzeCutsOnTree(hypergraph, incidence_struct, hypergraph.vwts, fixed_vtxs, tree, tree_cut_token.hyperedges_flag)
    clique_adj = buildCliqueGraph(hypergraph)
    
    hyperedge_wts = hypergraph.hwts
    n = hypergraph.n
    vtxCuts = cutInfo.vtxCuts
    edgeCuts = cutInfo.edgeCuts
    pred = cutInfo.pred
    pindex = cutInfo.pindex
    p1 = pindex.p1
    p2 = pindex.p2
    forced_0 = cutInfo.forced_0
    forced_1 = cutInfo.forced_1
    forced_01 = cutInfo.forced_01
    FB0 = cutInfo.FB0
    FB1 = cutInfo.FB1
    edgeCuts0 = cutInfo.edgeCuts0
    edgeCuts1 = cutInfo.edgeCuts1
    tree_cut_token.pred = cutInfo.pred
    tree_cut_token.nforced_0 = sum(hyperedge_wts[forced_0])
    tree_cut_token.nforced_1 = sum(hyperedge_wts[forced_1])
    tree_cut_token.nforced_01 = sum(hyperedge_wts[forced_01])
    total_vwts = sum(hypergraph.vwts)

    for i in 1:hypergraph.n
        tree_cut_token.exc0[i] = edgeCuts[i] + tree_cut_token.nforced_0 - FB0[i] + edgeCuts1[i] + tree_cut_token.nforced_01
        tree_cut_token.exc1[i] = edgeCuts[i] + tree_cut_token.nforced_1 - FB1[i] + edgeCuts0[i] + tree_cut_token.nforced_01
        # tree_cut_token.exc0[i] = - edgeCuts[i] 
        # tree_cut_token.exc1[i] = - edgeCuts[i]
        # println(tree_cut_token.exc0[i], " ", tree_cut_token.exc1[i])
    end

    tree_cut_token.exc0[1] = tree_cut_token.exc1[1] = hypergraph.e

    for i in 1:length(tree_cut_token.exc0)
        tree_cut_token.cutCost0[i] = tree_cut_token.exc0[i]
        tree_cut_token.cutCost1[i] = tree_cut_token.exc1[i]
        (tree_cut_token.cutCost[i], pol) = findmin([tree_cut_token.cutCost0[i], tree_cut_token.cutCost1[i]])
        tree_cut_token.polarity[i] = pol - 1

        if pol == 0
            tree_cut_token.areaUtil0[i] = vtxCuts[i]
            tree_cut_token.areaUtil1[i] = total_vwts - vtxCuts[i]

        else
            tree_cut_token.areaUtil1[i] = vtxCuts[i]
            tree_cut_token.areaUtil0[i] = total_vwts - vtxCuts[i]
        end

        if tree_cut_token.areaUtil0[i] > capacities[1] || tree_cut_token.areaUtil1[i] > capacities[2]
            tree_cut_token.areaCost[i] = 1e09
        end

        tree_cut_token.ratioCost[i] = tree_cut_token.cutCost[i] / (tree_cut_token.areaUtil0[i] * tree_cut_token.areaUtil1[i])
        tree_cut_token.tCost_v[i] = tree_cut_token.cutCost[i] + tree_cut_token.areaCost[i]
    end

    return cutInfo
end

function FindBestCutOnTree(tree::SimpleWeightedGraphs.SimpleGraph, hypergraph::Hypergraph, incidence_struct::Incidence, fixed_vtxs::Pindex, capacities::Vector{Int}, ncuts::Int)
    tree_cut_token = TreeCutToken(0, 0, 0, 0, zeros(Int, hypergraph.n), zeros(Int, hypergraph.n), zeros(Int, 2), zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n),
        zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n), zeros(Int, hypergraph.n), zeros(Int, hypergraph.n),
        zeros(Int, hypergraph.n), zeros(Int, hypergraph.n), zeros(Int, hypergraph.n), ones(Int, hypergraph.e), zeros(Float64, hypergraph.n), zeros(Float64, hypergraph.n), 0.5)
    done = false
    tree_copy = deepcopy(tree)
    part_vector = zeros(Int, hypergraph.n)
    cut_size = -1
    part_area = zeros(Int, 2)
    ci = TreeBipartition(tree_copy, tree_cut_token, hypergraph, incidence_struct, fixed_vtxs, capacities, ncuts)
    cut_size, cut_idx = findmin(tree_cut_token.tCost_v)
    components = Vector{Vector{Int}}()

    if cut_size == -1 #>= 1e09
        overflow = true
        edges_removed = Vector{Vector{Int}}()

        while overflow == true
            overflow = false
            (~, cut_idx) = findmin(tree_cut_token.ratioCost)
            rem_edge!(tree_copy, cut_idx, tree_cut_token.pred[cut_idx])
            push!(edges_removed, [cut_idx, tree_cut_token.pred[cut_idx]])
            components = connected_components(tree_copy)

            i = 0

            clusters = findLabels(components, hypergraph.n)

            for vertices in components
                total_weight = sum(hypergraph.vwts[vertices])

                i += 1

                println(i, ": ", total_weight, "::", capacities[1])

                if total_weight > capacities[1]
                    overflow = true
                    break
                end
            end

            #println("Final overflow: ", overflow)

            if overflow == false
                for edge in edges_removed
                    add_edge!(tree_copy, edge[1], edges[2])
                end
                break
            end

            (cut_hyperedges_mrk, ~, ~, ~, ~) = cutProfile(hypergraph, incidence_struct, clusters)
            cut_hyperedges = findall(!iszero, cut_hyperedges_mrk)
            tree_cut_token.hyperedges_flag[cut_hyperedges] .= 0
            ci = TreeBipartition(tree, tree_cut_token, hypergraph, incidence_struct, fixed_vtxs, capacities, ncuts)

            return ci
        end

        vwts_cc = ContractVtxWts(hypergraph.vwts, clusters)
        n_cc, e_cc, hedges_cc, eptr_cc, hwts_cc = ContractHyperGraph(hypergraph, clusters)
        hypergraph_cc = Hypergraph(n_cc, e_cc, hedges_cc, eptr_cc, vwts_cc, hwts_cc)
        fixed_part_cc = -ones(Int, hypergraph_cc.n)
        part_vector = PartitionILP(H_cc, fixed_part_cc, capacities)
        cut_size, part_area = FindCutSize(partition_vector, hypergraph_cc)
    else
        done = true
        SimpleWeightedGraphs.rem_edge!(tree, cut_idx, tree_cut_token.pred[cut_idx])
        components = SimpleWeightedGraphs.connected_components(tree)
        part_vector = findLabels(components, hypergraph.n)
        SimpleWeightedGraphs.add_edge!(tree, cut_idx, tree_cut_token.pred[cut_idx])
        part_area[1] = tree_cut_token.areaUtil0[cut_idx]
        part_area[2] = tree_cut_token.areaUtil1[cut_idx]
    end

    @info "[TREE CUT] BEST CUT RECORDED ON TREE: $cut_size WITH AREA SPLIT: $(part_area[1]) and $(part_area[2])"

    return part_vector .- 1, cut_size
end

function buildCliqueGraph(hypergraph::Hypergraph)
    n = hypergraph.n
    e = hypergraph.e
    hedges = hypergraph.hedges
    eptr = hypergraph.eptr
    hwts = hypergraph.hwts

    # 使用稀疏矩阵构建 clique graph
    clique_adj = spzeros(Int, n, n)

    for i in 1:e
        start_idx = eptr[i]
        end_idx = eptr[i+1] - 1
        vertices = hedges[start_idx:end_idx]
        weight = hwts[i]

        # 在 clique 中的每对顶点之间添加边
        for j in 1:length(vertices)
            for k in j+1:length(vertices)
                u = vertices[j]
                v = vertices[k]
                clique_adj[u, v] += weight
                clique_adj[v, u] += weight
            end
        end
    end

    return clique_adj
end

function computeCliqueCutCost(clique_adj::SparseMatrixCSC{Int,Int}, cut_vertex::Int,
    vtxCuts::Vector{Int}, total_vwts::Int, polarity::Int)
    n = size(clique_adj, 1)
    clique_cut_cost = 0.0

    # 根据极性确定分割点
    if polarity == 0
        split_point = vtxCuts[cut_vertex]
    else
        split_point = total_vwts - vtxCuts[cut_vertex]
    end

    # 计算每个顶点属于哪个分区
    vertex_partition = zeros(Int, n)
    current_sum = 0

    # 这里需要根据树的遍历顺序来确定顶点分区
    # 简化版本：假设顶点按某种顺序排列
    for i in 1:n
        current_sum += 1  # 假设每个顶点权重为1
        if current_sum <= split_point
            vertex_partition[i] = 0
        else
            vertex_partition[i] = 1
        end
    end

    # 计算跨越分区的边权重
    for i in 1:n
        for j in i+1:n
            if clique_adj[i, j] > 0 && vertex_partition[i] != vertex_partition[j]
                clique_cut_cost += clique_adj[i, j]
            end
        end
    end

    return clique_cut_cost
end