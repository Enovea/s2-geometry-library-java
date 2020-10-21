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
package dilivia.s2.shape

import dilivia.s2.Assertions.assertEQ
import dilivia.s2.Assertions.assertLT
import dilivia.s2.S2Point
import dilivia.s2.region.S2Polyline
import dilivia.s2.shape.TypeTags.kLaxPolylineTypeTag
import mu.KotlinLogging
import kotlin.math.max
import kotlin.math.min

// S2LaxPolylineShape represents a polyline.  It is similar to
// S2Polyline::Shape except that duplicate vertices are allowed, and the
// representation is slightly more compact.
//
// Polylines may have any number of vertices, but note that polylines with
// fewer than 2 vertices do not define any edges.  (To create a polyline
// consisting of a single degenerate edge, either repeat the same vertex twice
// or use S2LaxClosedPolylineShape defined in s2_lax_loop_shape.h.)
// Constructs an S2LaxPolylineShape with the given vertices.
class S2LaxPolylineShape(val vertices: List<S2Point>) : S2Shape() {

    // Constructs an empty polyline.
    constructor() : this(emptyList()) {}

    // Constructs an S2LaxPolylineShape from the given S2Polyline, by copying
    // its data.
    constructor(polyline: S2Polyline) : this(polyline.vertices())

    // Initializes an S2LaxPolylineShape with the given vertices.
    init {
        if (vertices.size == 1) {
            logger.warn { "S2LaxPolylineShape with one vertex has no edges" }
        }
    }

    fun numVertices(): Int = vertices.size
    fun vertex(i: Int): S2Point = vertices[i]

    // S2Shape interface:
    override val numEdges: Int = max(0, numVertices() - 1)
    override fun edge(edgeId: Int): Edge {
        assertLT(edgeId, numEdges)
        return Edge(vertex(edgeId), vertex(edgeId + 1))
    }

    override val dimension: Int = 1
    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)
    override val numChains: Int = min(1, numEdges)
    override fun chain(chain_id: Int): Chain = Chain(0, numEdges)

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        assertEQ(chainId, 0)
        assertLT(offset, numEdges)
        return Edge(vertex(offset), vertex(offset + 1))
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

    override val typeTag: TypeTag = kLaxPolylineTypeTag

    companion object {
        val logger = KotlinLogging.logger(S2LaxPolylineShape::class.java.name)
    }
}

// TODO
/*
// Exactly like S2LaxPolylineShape, except that the vertices are kept in an
// encoded form and are decoded only as they are accessed.  This allows for
// very fast initialization and no additional memory use beyond the encoded
// data.  The encoded data is not owned by this class; typically it points
// into a large contiguous buffer that contains other encoded data as well.
class EncodedS2LaxPolylineShape : public S2Shape {
 public:
  // Constructs an uninitialized object; requires Init() to be called.
  EncodedS2LaxPolylineShape() {}

  // Initializes an EncodedS2LaxPolylineShape.
  //
  // REQUIRES: The Decoder data buffer must outlive this object.
  bool Init(Decoder* decoder);

  int num_vertices() const { return vertices_.size(); }
  S2Point vertex(int i) const { return vertices_[i]; }

  // S2Shape interface:
  int num_edges() const final { return std::max(0, num_vertices() - 1); }
  Edge edge(int e) const final;
  int dimension() const final { return 1; }
  ReferencePoint GetReferencePoint() const final {
    return ReferencePoint::Contained(false);
  }
  int num_chains() const final;
  Chain chain(int i) const final;
  Edge chain_edge(int i, int j) const final;
  ChainPosition chain_position(int e) const final;

 private:
  s2coding::EncodedS2PointVector vertices_;


}
*/


