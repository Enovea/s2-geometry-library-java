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
package dilivia.s2.index

import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.builder.DegenerateEdges
import dilivia.s2.builder.DuplicateEdges
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.Graph
import dilivia.s2.builder.GraphOptions
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.Label
import dilivia.s2.builder.LabelSetId
import dilivia.s2.builder.layers.Layer
import dilivia.s2.builder.SiblingPairs
import dilivia.s2.index.S2PointVectorLayer.Options
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.shape.S2PointVectorShape


typealias LabelSetIds = MutableList<LabelSetId>

// A layer type that collects degenerate edges as points.
// This layer expects all edges to be degenerate. In case of finding
// non-degenerate edges it sets S2Error but it still generates the
// output with degenerate edges.
class S2PointVectorLayer(
        private val points: MutableList<S2Point>,
        private val labelSetIds: LabelSetIds? = null,
        private val labelSetLexicon: IdSetLexicon? = null,
        private val options: Options = Options()
) : Layer() {

  data class Options(var duplicateEdges: DuplicateEdges = DuplicateEdges.MERGE)

  // Layer interface:

    override fun graphOptions(): GraphOptions = GraphOptions(
            edge_type = EdgeType.DIRECTED,
            degenerate_edges = DegenerateEdges.KEEP,
            duplicate_edges = options.duplicateEdges,
            sibling_pairs = SiblingPairs.KEEP
    )

    override fun build(g: Graph, error: S2Error) {
        val fetcher = Graph.LabelFetcher(g, EdgeType.DIRECTED)

        val labels = mutableListOf<Label>()  // Temporary storage for labels.
        for (edge_id in g.edges.indices) {
            val edge = g.edge(edge_id)
            if (edge.first != edge.second) {
                error.init(S2Error.INVALID_ARGUMENT, "Found non-degenerate edges");
                continue;
            }
            points.add(g.vertex(edge.first));
            if (labelSetIds != null && labelSetLexicon != null) {
                fetcher.fetch(edge_id, labels)
                val set_id = labelSetLexicon.add(labels)
                labelSetIds.add(set_id)
            }
        }
    }

}

// Like S2PointVectorLayer, but adds the points to a MutableS2ShapeIndex (if
// the point vector is non-empty).
class IndexedS2PointVectorLayer(
        private val index: MutableS2ShapeIndex,
        private val options: Options = Options()
) : Layer() {

    private val points: MutableList<S2Point> = mutableListOf()
    private val layer: S2PointVectorLayer = S2PointVectorLayer(points, options = options)

    override fun graphOptions(): GraphOptions = layer.graphOptions()

    override fun build(g: Graph, error: S2Error) {
        layer.build(g, error)
        if (error.isOk() && points.isNotEmpty()) {
            index.add(S2PointVectorShape(points = points))
        }
    }

}
