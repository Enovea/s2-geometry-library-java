/**
 * This project is a kotlin port of the Google s2 geometry library (Copyright 2005 Google Inc. All Rights Reserved.):
 *                                 https://github.com/google/s2geometry.git
 *
 * Copyright © 2020 Dilivia (contact@dilivia.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package dilivia.s2.builder.graph

import com.google.common.collect.ComparisonChain
import dilivia.s2.Assertions.assertEQ
import dilivia.s2.Assertions.assertTrue
import dilivia.s2.S2Error
import dilivia.s2.builder.Edge
import dilivia.s2.builder.EdgeId
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.InputEdgeId
import dilivia.s2.builder.InputEdgeIdSetId
import dilivia.s2.collections.assign
import dilivia.s2.collections.reverse
import kotlin.math.max


class EdgeProcessor(
        val options: GraphOptions,
        val edges: ArrayList<Edge>,
        val input_ids: ArrayList<InputEdgeIdSetId>,
        val id_set_lexicon: IdSetLexicon
) {

    private val out_edges_ = mutableListOf<EdgeId>()
    private val in_edges_ = mutableListOf<EdgeId>()

    private val new_edges_ = ArrayList<Edge>()
    private val new_input_ids_ = ArrayList<InputEdgeIdSetId>()

    private val tmp_ids_ = mutableListOf<InputEdgeId>()

    init {
        // Sort the outgoing and incoming edges in lexigraphic order.  We use a
        // stable sort to ensure that each undirected edge becomes a sibling pair,
        // even if there are multiple identical input edges.
        out_edges_.assign(edges.size, 0)
        out_edges_.sortWith { a, b -> if (Graph.stableLessThan(edges[a], edges[b], a, b)) -1 else 1 }
        in_edges_.assign(edges.size, 0)
        in_edges_.sortWith { a, b -> if (Graph.stableLessThan(edges[a].reverse(), edges[b].reverse(), a, b)) -1 else 1 }

        new_edges_.ensureCapacity(edges.size)
        new_input_ids_.ensureCapacity(edges.size)
    }

    fun run(error: S2Error) {
        val num_edges = edges.size
        if (num_edges == 0) return

        // Walk through the two sorted arrays performing a merge join.  For each
        // edge, gather all the duplicate copies of the edge in both directions
        // (outgoing and incoming).  Then decide what to do based on "options_" and
        // how many copies of the edge there are in each direction.
        var outIdx = 0
        var inIdx = 0
        var out_edge = edges[out_edges_[outIdx]]
        var in_edge = edges[in_edges_[inIdx]]
        val sentinel = Edge(VertexId.MAX_VALUE, VertexId.MAX_VALUE)
        while (true) {
            val edge = minOf(out_edge, in_edge.reverse()) { e1, e2 -> ComparisonChain.start().compare(e1.first, e1.second).compare(e1.second, e2.second).result() }
            if (edge == sentinel) break

            val out_begin = outIdx
            val in_begin = inIdx
            while (out_edge == edge) {
                out_edge = if (++outIdx == num_edges) sentinel else edges[out_edges_[outIdx]]
            }
            while (in_edge.reverse() == edge) {
                in_edge = if (++inIdx == num_edges) sentinel else edges[in_edges_[inIdx]]
            }
            val n_out = outIdx - out_begin
            val n_in = inIdx - in_begin
            if (edge.first == edge.second) {
                assertEQ(n_out, n_in);
                if (options.degenerate_edges == DegenerateEdges.DISCARD) {
                    continue
                }
                if (options.degenerate_edges == DegenerateEdges.DISCARD_EXCESS &&
                        ((out_begin > 0 && edges[out_edges_[out_begin - 1]].first == edge.first) ||
                                (outIdx < num_edges && edges[out_edges_[outIdx]].first == edge.first) ||
                                (in_begin > 0 && edges[in_edges_[in_begin - 1]].second == edge.first) ||
                                (inIdx < num_edges && edges[in_edges_[inIdx]].second == edge.first))) {
                    continue  // There were non-degenerate incident edges, so discard.
                }
                if (options.edge_type == EdgeType.UNDIRECTED &&
                        (options.sibling_pairs == SiblingPairs.REQUIRE || options.sibling_pairs == SiblingPairs.CREATE)) {
                    // When we have undirected edges and are guaranteed to have siblings,
                    // we cut the number of edges in half (see s2builder.h).
                    assertEQ(0, n_out and 1)  // Number of edges is always even.
                    addEdges(if (options.duplicate_edges == DuplicateEdges.MERGE) 1 else (n_out / 2), edge, mergeInputIds(out_begin, outIdx))
                } else if (options.duplicate_edges == DuplicateEdges.MERGE) {
                    addEdges(if (options.edge_type == EdgeType.UNDIRECTED) 2 else 1, edge, mergeInputIds(out_begin, outIdx))
                } else if (options.sibling_pairs == SiblingPairs.DISCARD || options.sibling_pairs == SiblingPairs.DISCARD_EXCESS) {
                    // Any SiblingPair option that discards edges causes the labels of all
                    // duplicate edges to be merged together (see s2builder.h).
                    addEdges(n_out, edge, mergeInputIds(out_begin, outIdx))
                } else {
                    copyEdges(out_begin, outIdx)
                }
            } else if (options.sibling_pairs == SiblingPairs.KEEP) {
                if (n_out > 1 && options.duplicate_edges == DuplicateEdges.MERGE) {
                    addEdge(edge, mergeInputIds(out_begin, outIdx))
                } else {
                    copyEdges(out_begin, outIdx)
                }
            } else if (options.sibling_pairs == SiblingPairs.DISCARD) {
                if (options.edge_type == EdgeType.DIRECTED) {
                    // If n_out == n_in: balanced sibling pairs
                    // If n_out < n_in:  unbalanced siblings, in the form AB, BA, BA
                    // If n_out > n_in:  unbalanced siblings, in the form AB, AB, BA
                    if (n_out <= n_in) continue
                    // Any option that discards edges causes the labels of all duplicate
                    // edges to be merged together (see s2builder.h).
                    addEdges(if (options.duplicate_edges == DuplicateEdges.MERGE) 1 else (n_out - n_in), edge, mergeInputIds(out_begin, outIdx))
                } else {
                    if ((n_out and 1) == 0) continue
                    addEdge(edge, mergeInputIds(out_begin, outIdx))
                }
            } else if (options.sibling_pairs == SiblingPairs.DISCARD_EXCESS) {
                if (options.edge_type == EdgeType.DIRECTED) {
                    // See comments above.  The only difference is that if there are
                    // balanced sibling pairs, we want to keep one such pair.
                    if (n_out < n_in) continue
                    addEdges(if (options.duplicate_edges == DuplicateEdges.MERGE) 1 else max(1, n_out - n_in), edge, mergeInputIds(out_begin, outIdx))
                } else {
                    addEdges(if ((n_out and 1) != 0) 1 else 2, edge, mergeInputIds(out_begin, outIdx))
                }
            } else {
                assertTrue(options.sibling_pairs == SiblingPairs.REQUIRE || options.sibling_pairs == SiblingPairs.CREATE)
                if (error.isOk() && options.sibling_pairs == SiblingPairs.REQUIRE &&
                        if (options.edge_type == EdgeType.DIRECTED) (n_out != n_in) else ((n_out and 1) != 0)
                ) {
                    error.init(S2Error.BUILDER_MISSING_EXPECTED_SIBLING_EDGES,
                            "Expected all input edges to have siblings,  but some were missing")
                }
                if (options.duplicate_edges == DuplicateEdges.MERGE) {
                    addEdge(edge, mergeInputIds(out_begin, outIdx))
                } else if (options.edge_type == EdgeType.UNDIRECTED) {
                    // Convert graph to use directed edges instead (see documentation of
                    // REQUIRE/CREATE for undirected edges).
                    addEdges((n_out + 1) / 2, edge, mergeInputIds(out_begin, outIdx))
                } else {
                    copyEdges(out_begin, outIdx)
                    if (n_in > n_out) {
                        // Automatically created edges have no input edge ids or labels.
                        addEdges(n_in - n_out, edge, IdSetLexicon.emptySetId())
                    }
                }
            }
        }

        edges.clear()
        edges.addAll(new_edges_)
        edges.trimToSize()
        input_ids.clear()
        input_ids.addAll(new_input_ids_)
        input_ids.trimToSize()
    }


    private fun addEdge(edge: Edge, inputEdgeIdSetId: InputEdgeIdSetId) {
        new_edges_.add(edge)
        new_input_ids_.add(inputEdgeIdSetId)
    }

    private fun addEdges(numEdges: Int, edge: Edge, inputEdgeIdSetId: InputEdgeIdSetId) {
        for (i in 0 until numEdges) {
            addEdge(edge, inputEdgeIdSetId)
        }
    }

    private fun copyEdges(out_begin: Int, out_end: Int) {
        for (i in out_begin until out_end) {
            addEdge(edges[out_edges_[i]], input_ids[out_edges_[i]])
        }
    }

    private fun mergeInputIds(out_begin: Int, out_end: Int): InputEdgeIdSetId {
        if (out_end - out_begin == 1) {
            return input_ids[out_edges_[out_begin]]
        }
        tmp_ids_.clear()
        for (i in out_begin until out_end) {
            for (id in id_set_lexicon.idSet(input_ids[out_edges_[i]])) {
                tmp_ids_.add(id)
            }
        }
        return id_set_lexicon.add(tmp_ids_)
    }

}
