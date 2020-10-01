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

import dilivia.s2.Assertions
import dilivia.s2.S2Point
import dilivia.s2.region.S2ContainsVertexQuery

// A 32-bit tag that can be used to identify the type of an encoded S2Shape.
// All encodable types have a non-zero type tag.  The tag associated with a
// given shape type can be accessed as Shape::kTypeTag, while the tag
// associated with a given object can be accessed as shape.type_tag().
//
// Type tags in the range 0..8191 are reserved for use by the S2 library.
typealias TypeTag = UInt

object TypeTags {
    val kNoTypeTag: TypeTag = 0U
    val kPolylineTypeTag: TypeTag = 2U
    val kPointVectorTypeTag: TypeTag = 3U
}

// An edge, consisting of two vertices "v0" and "v1".  Zero-length edges are
// allowed, and can be used to represent points.
data class Edge(val v0: S2Point = S2Point(), val v1: S2Point = S2Point()): Comparable<Edge> {

    override fun compareTo(other: Edge): Int {
        val v0Comparison = v0.compareTo(other.v0)
        return if (v0Comparison == 0) {
            v1.compareTo(other.v1)
        } else v0Comparison
    }

};
// The purpose of S2Shape is to represent polygonal geometry in a flexible
// way.  It is organized as a collection of edges that optionally defines an
// interior.  All geometry represented by an S2Shape must have the same
// dimension, which means that an S2Shape can represent either a set of
// points, a set of polylines, or a set of polygons.
//
// S2Shape is defined as an abstract base class in order to give clients
// control over the underlying data representation.  Sometimes an S2Shape does
// not have any data of its own, but instead "wraps" some other class.  There
// are various useful subtypes defined in *_shape.h, and some S2 classes also
// have a nested "Shape" class (e.g., S2Polygon::Shape).  It is easy for
// clients to implement their own subtypes, since the interface is minimal.
//
// S2Shape operations are typically defined on S2ShapeIndex objects rather
// than individual shapes.  An S2ShapeIndex is simply a collection of
// S2Shapes, possibly of different dimensions (e.g. 10 points and 3 polygons),
// organized into a data structure for efficient edge access.
//
// The edges of an S2Shape are identified by a contiguous range of "edge ids"
// starting at 0.  The edges are further subdivided into "chains", where each
// chain consists of a sequence of edges connected end-to-end (a polyline).
// For example, an S2Shape representing two polylines AB and CDE would have
// three edges (AB, CD, DE) grouped into two chains: (AB) and (CD, DE).
// Similarly, an S2Shape representing 5 points would have 5 chains consisting
// of one edge each.
//
// S2Shape has methods that allow edges to be accessed either using the global
// numbering (edge id) or within a particular chain.  The global numbering is
// sufficient for most purposes, but the chain representation is useful for
// certain algorithms such as intersection (see S2BooleanOperation).
abstract class S2Shape(var id: Int = -1) {

    // A range of edge ids corresponding to a chain of zero or more connected
    // edges, specified as a (start, length) pair.  The chain is defined to
    // consist of edge ids {start, start + 1, ..., start + length - 1}.
    data class Chain(val start: Int, val length: Int)

    // The position of an edge within a given edge chain, specified as a
    // (chain_id, offset) pair.  Chains are numbered sequentially starting from
    // zero, and offsets are measured from the start of each chain.
    data class ChainPosition(val chainId: Int, val offset: Int)

    // A ReferencePoint consists of a point P and a boolean indicating whether P
    // is contained by a particular shape.
    data class ReferencePoint(val point: S2Point = S2Point.origin(), val contained: Boolean)

    // Returns the number of edges in this shape.  Edges have ids ranging from 0
    // to num_edges() - 1.
    abstract val numEdges: Int

    // Returns the endpoints of the given edge id.
    //
    // REQUIRES: 0 <= id < num_edges()
    abstract fun edge(edgeId: Int): Edge

    // Returns the dimension of the geometry represented by this shape.
    //
    //  0 - Point geometry.  Each point is represented as a degenerate edge.
    //
    //  1 - Polyline geometry.  Polyline edges may be degenerate.  A shape may
    //      represent any number of polylines.  Polylines edges may intersect.
    //
    //  2 - Polygon geometry.  Edges should be oriented such that the polygon
    //      interior is always on the left.  In theory the edges may be returned
    //      in any order, but typically the edges are organized as a collection
    //      of edge chains where each chain represents one polygon loop.
    //      Polygons may have degeneracies (e.g., degenerate edges or sibling
    //      pairs consisting of an edge and its corresponding reversed edge).
    //      A polygon loop may also be full (containing all points on the
    //      sphere); by convention this is represented as a chain with no edges.
    //      (See S2LaxPolygonShape for details.)
    //
    // Note that this method allows degenerate geometry of different dimensions
    // to be distinguished, e.g. it allows a point to be distinguished from a
    // polyline or polygon that has been simplified to a single point.
    abstract val dimension: Int

    // Returns true if the shape contains no points.  (Note that the full
    // polygon is represented as a chain with zero edges.)
    fun isEmpty(): Boolean {
        return numEdges == 0 && (dimension < 2 || numChains == 0);
    }
    // Returns true if the shape contains all points on the sphere.
    fun isFull(): Boolean {
        return numEdges == 0 && dimension == 2 && numChains > 0;
    }

    // Returns an arbitrary point P along with a boolean indicating whether P is
    // contained by the shape.  (The boolean value must be false for shapes that
    // do not have an interior.)
    //
    // This ReferencePoint may then be used to compute the containment of other
    // points by counting edge crossings.
    abstract fun getReferencePoint(): ReferencePoint

    // Returns the number of contiguous edge chains in the shape.  For example,
    // a shape whose edges are [AB, BC, CD, AE, EF] would consist of two chains
    // (AB,BC,CD and AE,EF).  Every chain is assigned a "chain id" numbered
    // sequentially starting from zero.
    //
    // Note that it is always acceptable to implement this method by returning
    // num_edges() (i.e. every chain consists of a single edge), but this may
    // reduce the efficiency of some algorithms.
    abstract val numChains: Int

    // Returns the range of edge ids corresponding to the given edge chain.  The
    // edge chains must form contiguous, non-overlapping ranges that cover the
    // entire range of edge ids.  This is spelled out more formally below:
    //
    // REQUIRES: 0 <= i < num_chains()
    // REQUIRES: chain(i).length >= 0, for all i
    // REQUIRES: chain(0).start == 0
    // REQUIRES: chain(i).start + chain(i).length == chain(i+1).start,
    //           for i < num_chains() - 1
    // REQUIRES: chain(i).start + chain(i).length == num_edges(),
    //           for i == num_chains() - 1
    abstract fun chain(chain_id: Int): Chain

    // Returns the edge at offset "offset" within edge chain "chain_id".
    // Equivalent to "shape.edge(shape.chain(chain_id).start + offset)"
    // but may be more efficient.
    abstract fun chainEdge(chainId: Int, offset: Int): Edge

    // Finds the chain containing the given edge, and returns the position of
    // that edge as a (chain_id, offset) pair.
    //
    // REQUIRES: shape.chain(pos.chain_id).start + pos.offset == edge_id
    // REQUIRES: shape.chain(pos.chain_id + 1).start > edge_id
    //
    // where     pos == shape.chain_position(edge_id).
    abstract fun chainPosition(edgeId: Int): ChainPosition

    // Returns an integer that can be used to identify the type of an encoded
    // S2Shape (see TypeTag above).
    abstract val typeTag: TypeTag

    // Virtual methods that return pointers of your choice.  These methods are
    // intended to help with the problem of attaching additional data to S2Shape
    // objects.  For example, you could return a pointer to a source object, or
    // a pointer to a bundle of additional data allocated with the S2Shape.
    // Because this method exists in all S2Shapes, you can override it in each
    // type of shape you have and call it without knowing the concrete subtype.
    // For example, if you have polyline and polygon shapes, you can do this:
    //
    //   class MyPolyline : public S2Polyline::Shape {
    //    public:
    //     virtual void* mutable_user_data() { return &my_data_; }
    //    private:
    //     MyData my_data_;
    //   };
    //   class MyPolygon : public S2Polygon::Shape {
    //    public:
    //     virtual void* mutable_user_data() { return &my_data_; }
    //    private:
    //     MyData my_data_;
    //   };
    //   ...
    //   S2Shape* shape = index.shape(id);
    //   const MyData* data = static_cast<const MyData*>(shape->user_data());
    //
    // This is not the only way to map from an S2Shape back to your source
    // data.  Other reasonable techniques include:
    //
    //  - Every shape has an id() assigned by S2ShapeIndex.  Ids are assigned
    //    sequentially starting from 0 in the order the shapes are added to the
    //    index.  You can use this id to look up arbitrary data stored in your
    //    own vector.
    //
    //  - If all of your shapes are the same type, then you can create your own
    //    subclass of some existing S2Shape type (such as S2Polyline::Shape) and
    //    add your own methods and fields.  You can access this data by
    //    downcasting the S2Shape pointers returned by S2ShapeIndex methods.
    open fun userData(): Any? = null

    companion object {

        // Indicates that a given S2Shape type cannot be encoded.
        val kNoTypeTag: TypeTag = 0U

        // The minimum allowable tag for user-defined S2Shape types.
        val kMinUserTypeTag: TypeTag = 8192.toUInt()
        // This is a helper function for implementing S2Shape::GetReferencePoint().
        //
        // Given a shape consisting of closed polygonal loops, the interior of the
        // shape is defined as the region to the left of all edges (which must be
        // oriented consistently).  This function then chooses an arbitrary point and
        // returns true if that point is contained by the shape.
        //
        // Unlike S2Loop and S2Polygon, this method allows duplicate vertices and
        // edges, which requires some extra care with definitions.  The rule that we
        // apply is that an edge and its reverse edge "cancel" each other: the result
        // is the same as if that edge pair were not present.  Therefore shapes that
        // consist only of degenerate loop(s) are either empty or full; by convention,
        // the shape is considered full if and only if it contains an empty loop (see
        // S2LaxPolygonShape for details).
        //
        // Determining whether a loop on the sphere contains a point is harder than
        // the corresponding problem in 2D plane geometry.  It cannot be implemented
        // just by counting edge crossings because there is no such thing as a "point
        // at infinity" that is guaranteed to be outside the loop.
        fun getReferencePoint(shape: S2Shape): ReferencePoint {
            Assertions.assertEQ(shape.dimension, 2)
            if (shape.numEdges == 0) {
                // A shape with no edges is defined to be full if and only if it
                // contains at least one chain.
                return ReferencePoint(contained = shape.numChains > 0)
            }
            // Define a "matched" edge as one that can be paired with a corresponding
            // reversed edge.  Define a vertex as "balanced" if all of its edges are
            // matched. In order to determine containment, we must find an unbalanced
            // vertex.  Often every vertex is unbalanced, so we start by trying an
            // arbitrary vertex.
            val edge = shape.edge(0)
            var result = getReferencePointAtVertex(shape, edge.v0)
            if (result != null) {
                return result
            }
            // That didn't work, so now we do some extra work to find an unbalanced
            // vertex (if any).  Essentially we gather a list of edges and a list of
            // reversed edges, and then sort them.  The first edge that appears in one
            // list but not the other is guaranteed to be unmatched.
            val n = shape.numEdges
            val edges = ArrayList<Edge>(n)
            val rev_edges = ArrayList<Edge>(n)
            for (i in 0 until n) {
                val e = shape.edge(i)
                edges.add(e)
                rev_edges.add(Edge(e.v1, e.v0))
            }
            edges.sort()
            rev_edges.sort()
            for (i in 0 until n) {
                if (edges[i] < rev_edges[i]) {  // edges[i] is unmatched
                    result = getReferencePointAtVertex(shape, edges[i].v0)
                    check(result != null)
                    return result;
                }
                if (rev_edges[i] < edges[i]) {  // rev_edges[i] is unmatched
                    result = getReferencePointAtVertex(shape, rev_edges[i].v0)
                    check(result != null)
                    return result;
                }
            }
            // All vertices are balanced, so this polygon is either empty or full except
            // for degeneracies.  By convention it is defined to be full if it contains
            // any chain with no edges.
            for (i in 0 until shape.numChains) {
                if (shape.chain(i).length == 0) return ReferencePoint(contained = true)
            }
            return ReferencePoint(contained = false)
        }

        // This is a helper function for GetReferencePoint() above.
        //
        // If the given vertex "vtest" is unbalanced (see definition below), sets
        // "result" to a ReferencePoint indicating whther "vtest" is contained and
        // returns true.  Otherwise returns false.
        private fun getReferencePointAtVertex(shape: S2Shape, vtest: S2Point): ReferencePoint? {
            // Let P be an unbalanced vertex.  Vertex P is defined to be inside the
            // region if the region contains a particular direction vector starting from
            // P, namely the direction S2::Ortho(P).  This can be calculated using
            // S2ContainsVertexQuery.
            val contains_query = S2ContainsVertexQuery(vtest)
            val n = shape.numEdges
            for (e in 0 until n) {
                val edge = shape.edge(e)
                if (edge.v0 == vtest) contains_query.addEdge(edge.v1, 1)
                if (edge.v1 == vtest) contains_query.addEdge(edge.v0, -1)
            }
            val contains_sign = contains_query.containsSign()
            if (contains_sign == 0) {
                return null;  // There are no unmatched edges incident to this vertex.
            }
            return ReferencePoint(point = vtest, contained = contains_sign > 0)
        }
    }

}
