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
package dilivia.s2.builder.layers

import dilivia.s2.Assertions
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.builder.graph.DegenerateEdges
import dilivia.s2.builder.graph.DuplicateEdges
import dilivia.s2.builder.EdgeType
import dilivia.s2.builder.graph.Graph
import dilivia.s2.builder.graph.GraphOptions
import dilivia.s2.builder.IdSetLexicon
import dilivia.s2.builder.Label
import dilivia.s2.builder.graph.SiblingPairs
import dilivia.s2.index.LabelSetIds
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.region.S2Polyline


// A layer type that assembles edges (directed or undirected) into an
// S2Polyline.  Returns an error if the edges cannot be assembled into a
// single unbroken polyline.
//
// Duplicate edges are handled correctly (e.g., if a polyline backtracks on
// itself, or loops around and retraces some of its previous edges.)  The
// implementation attempts to preserve the order of directed input edges
// whenever possible, so that if the input is a polyline and it is not
// modified by S2Builder, then the output will be the same polyline (even if
// the polyline backtracks on itself or forms a loop).  With undirected edges,
// there are no such guarantees; for example, even if the input consists of a
// single undirected edge, then either directed edge may be returned.
//
// S2PolylineLayer does not support options such as discarding sibling pairs
// or merging duplicate edges because these options can split the polyline
// into several pieces.  Use S2PolylineVectorLayer if you need these features.
class S2PolylineLayer(
        private val polyline: S2Polyline,
        private val labelSetIds: LabelSetIds? = null,
        private val labelSetLexicon: IdSetLexicon? = null,
        private val options: Options = Options()
) : Layer() {

    init {
        Assertions.assertEQ(labelSetIds == null, labelSetLexicon == null)
    }

  data class Options(

          // Indicates whether the input edges provided to S2Builder are directed or
          // undirected.  Directed edges should be used whenever possible to avoid
          // ambiguity.
          //
          // DEFAULT: S2Builder::EdgeType::DIRECTED
          var edgeType: EdgeType = EdgeType.DIRECTED,

        // If true, calls FindValidationError() on the output polyline.  If any
        // error is found, it will be returned by S2Builder::Build().
        //
        // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
        // order to turn off the default error checking in debug builds.
        //
        // DEFAULT: false
        var validate: Boolean = false
  )


  // Layer interface:

    override fun graphOptions(): GraphOptions {
        // Remove edges that collapse to a single vertex, but keep duplicate and
        // sibling edges, since merging duplicates or discarding siblings can make
        // it impossible to assemble the edges into a single polyline.
        return GraphOptions(options.edgeType, DegenerateEdges.DISCARD, DuplicateEdges.KEEP, SiblingPairs.KEEP)
    }

    override fun build(g: Graph, error: S2Error) {
        if (g.numEdges() == 0) {
            polyline.init(mutableListOf())
            return
        }
        val edge_polylines = g.getPolylines(Graph.PolylineType.WALK)
        if (edge_polylines.size != 1) {
            error.init(S2Error.BUILDER_EDGES_DO_NOT_FORM_POLYLINE, "Input edges cannot be assembled into polyline")
            return
        }
        val edge_polyline = edge_polylines[0]
        val vertices = ArrayList<S2Point>(edge_polyline.size)  // Temporary storage for vertices.
        vertices.add(g.vertex(g.edge(edge_polyline[0]).first))
        for (e in edge_polyline) {
            vertices.add(g.vertex(g.edge(e).second))
        }
        if (labelSetIds != null && labelSetLexicon != null) {
            val fetcher = Graph.LabelFetcher(g, options.edgeType)
            val labels = mutableListOf<Label>();  // Temporary storage for labels.
            if (labelSetIds is ArrayList) {
                labelSetIds.ensureCapacity(edge_polyline.size)
            }
            for (e in edge_polyline) {
                fetcher.fetch(e, labels)
                labelSetIds.add(labelSetLexicon.add(labels))
            }
        }
        polyline.init(vertices, check = !options.validate)
        if (options.validate) {
            polyline.findValidationError(error)
        }
    }

}


// Like S2PolylineLayer, but adds the polyline to a MutableS2ShapeIndex (if the
// polyline is non-empty).
class IndexedS2PolylineLayer(
        private val index: MutableS2ShapeIndex,
        options: S2PolylineLayer.Options = S2PolylineLayer.Options()
) : Layer() {

    private val polyline: S2Polyline = S2Polyline()

    private val layer: S2PolylineLayer = S2PolylineLayer(polyline, options = options)

    override fun graphOptions(): GraphOptions = layer.graphOptions()

    override fun build(g: Graph, error: S2Error) {
        layer.build(g, error)
        if (error.isOk() && polyline.numVertices() > 0) {
            index.add(S2Polyline.Shape(polyline = polyline))
        }
    }

}


