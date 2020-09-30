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

import dilivia.s2.Assertions
import dilivia.s2.S2EdgeCrosser
import dilivia.s2.S2Error
import dilivia.s2.S2WedgeRelations
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.S2ShapeIndex
import dilivia.s2.shape.S2ShapeIndexCell
import dilivia.s2.shape.ShapeEdge
import mu.KotlinLogging

object S2CrossingEdgePairsScanner {

    private val logger = KotlinLogging.logger {  }

    // Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
    // early if the given EdgePairVisitor function returns false (in which case
    // VisitCrossings returns false as well).  "type" indicates whether all
    // crossings should be visited, or only interior crossings.
    //
    // CAVEAT: Crossings may be visited more than once.
    fun visitCrossingEdgePairs(index: S2ShapeIndex, type: CrossingType, visitor: EdgePairVisitor): Boolean = TODO()

    // Like the above, but visits all pairs of crossing edges where one edge comes
    // from each S2ShapeIndex.
    //
    // CAVEAT: Crossings may be visited more than once.
    fun visitCrossingEdgePairs(a_index: S2ShapeIndex, b_index: S2ShapeIndex, type: CrossingType, visitor: EdgePairVisitor): Boolean = TODO()

    // Given an S2ShapeIndex containing a single polygonal shape (e.g., an
    // S2Polygon or S2Loop), return true if any loop has a self-intersection
    // (including duplicate vertices) or crosses any other loop (including vertex
    // crossings and duplicate edges) and set "error" to a human-readable error
    // message.  Otherwise return false and leave "error" unchanged.
    //
    // This method is used to implement the FindValidationError methods of S2Loop
    // and S2Polygon.
    //
    // TODO(ericv): Add an option to support S2LaxPolygonShape rules (i.e.,
    // duplicate vertices and edges are allowed, but loop crossings are not).
    fun findSelfIntersection(index: S2ShapeIndex): S2Error {
        if (index.numShapeIds() == 0) return S2Error(code = S2Error.OK);
        Assertions.assertEQ(1, index.numShapeIds())
        val shape = index.shape(0)!!

        // Visit all crossing pairs except possibly for ones of the form (AB, BC),
        // since such pairs are very common and FindCrossingError() only needs pairs
        // of the form (AB, AC).
        var error = S2Error(S2Error.OK)
        visitCrossings(
                index, CrossingType.ALL, false /*need_adjacent*/,
                object : EdgePairVisitor {
                    override fun visit(a: ShapeEdge, b: ShapeEdge, is_interior: Boolean): Boolean {
                        error = findCrossingError(shape, a, b, is_interior)
                        return error.code == S2Error.OK
                    }

                });
        return error
    }

    // Appends all edges in the given S2ShapeIndexCell to the given vector.
    private fun appendShapeEdges(index: S2ShapeIndex, cell: S2ShapeIndexCell, shape_edges: ShapeEdgeVector) {
        for (s in 0 until cell.numClipped()) {
            val clipped = cell.clipped(s)
            val shape = index.shape(clipped.shapeId) ?: continue
            val num_edges = clipped.numEdges()
            for (i in 0 until num_edges) {
                shape_edges.add(ShapeEdge(shape, clipped.edge(i)));
            }
        }
    }

    // Returns a vector containing all edges in the given S2ShapeIndexCell.
    // (The result is returned as an output parameter so that the same storage can
    // be reused, rather than allocating a new temporary vector each time.)
    private fun getShapeEdges(index: S2ShapeIndex, cell: S2ShapeIndexCell, shape_edges: ShapeEdgeVector) {
        shape_edges.clear()
        appendShapeEdges(index, cell, shape_edges)
    }

    // Returns a vector containing all edges in the given S2ShapeIndexCell vector.
    // (The result is returned as an output parameter so that the same storage can
    // be reused, rather than allocating a new temporary vector each time.)
    fun getShapeEdges(index: S2ShapeIndex, cells: List<S2ShapeIndexCell>, shape_edges: ShapeEdgeVector) {
        shape_edges.clear()
        for (cell in cells) {
            appendShapeEdges(index, cell, shape_edges)
        }
    }

    // Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
    // early if the given EdgePairVisitor function returns false (in which case
    // VisitCrossings returns false as well).  "type" indicates whether all
    // crossings should be visited, or only interior crossings.
    //
    // If "need_adjacent" is false, then edge pairs of the form (AB, BC) may
    // optionally be ignored (even if the two edges belong to different edge
    // chains).  This option exists for the benefit of FindSelfIntersection(),
    // which does not need such edge pairs (see below).
    fun visitCrossings(index: S2ShapeIndex, type: CrossingType, need_adjacent: Boolean, visitor: EdgePairVisitor): Boolean {
        logger.trace { "Visit crossings(type = $type, needAdjacent = $need_adjacent)" }
        // TODO(ericv): Use brute force if the total number of edges is small enough
        // (using a larger threshold if the S2ShapeIndex is not constructed yet).
        val shape_edges = mutableListOf<ShapeEdge>()
        val iter = index.iterator(S2ShapeIndex.InitialPosition.BEGIN)
        while (!iter.done()) {
            logger.trace { "Iterate on cell: ${iter.cell()}" }
            getShapeEdges(index, iter.cell()!!, shape_edges)
            if (!visitCrossings(shape_edges, type, need_adjacent, visitor)) {
                return false
            }
            iter.next()
        }
        return true
    }

    // Given a vector of edges within an S2ShapeIndexCell, visit all pairs of
    // crossing edges (of the given CrossingType).
    private fun visitCrossings(shape_edges: ShapeEdgeVector, type: CrossingType, need_adjacent: Boolean, visitor: EdgePairVisitor): Boolean {
        val min_crossing_sign = if (type == CrossingType.INTERIOR) 1 else 0
        val num_edges = shape_edges.size
        var i = 0
        while (i + 1 < num_edges) {
            val a = shape_edges[i]
            var j = i + 1
            // A common situation is that an edge AB is followed by an edge BC.  We
            // only need to visit such crossings if "need_adjacent" is true (even if
            // AB and BC belong to different edge chains).
            if (!need_adjacent && a.v1 == shape_edges[j].v0) {
                if (++j >= num_edges) break
            }
            val crosser = S2EdgeCrosser(a.v0, a.v1);
            while (j < num_edges) {
                val b = shape_edges[j]
                if (crosser.c == null || crosser.c != b.v0) {
                    crosser.restartAt(b.v0)
                }
                val sign = crosser.crossingSign(b.v1)
                if (sign >= min_crossing_sign) {
                    if (!visitor.visit(a, b, sign == 1)) return false
                }
                ++j
            }
            ++i
        }
        return true
    }

    // Helper function that formats a loop error message.  If the loop belongs to
    // a multi-loop polygon, adds a prefix indicating which loop is affected.
    private fun loopError(code: Int, format: String, ap: S2Shape.ChainPosition, bp: S2Shape.ChainPosition, is_polygon: Boolean): S2Error {
        var message = format.format(ap.offset, bp.offset)
        if (is_polygon) {
            message = "Loop ${ap.chainId}: $message"
        }
        return S2Error(code = code, text = message)
    }

    // Given two loop edges that cross (including at a shared vertex), return true
    // if there is a crossing error and set "error" to a human-readable message.
    private fun findCrossingError(shape: S2Shape, a: ShapeEdge, b: ShapeEdge, is_interior: Boolean): S2Error {
        val error: S2Error
        val is_polygon = shape.numChains > 1
        val ap = shape.chainPosition(a.id.edgeId)
        val bp = shape.chainPosition(b.id.edgeId)
        if (is_interior) {
            return if (ap.chainId != bp.chainId) {
                S2Error(code = S2Error.POLYGON_LOOPS_CROSS, "Loop %d edge %d crosses loop %d edge %d".format(
                        ap.chainId, ap.offset, bp.chainId, bp.offset))
            } else {
                loopError(S2Error.LOOP_SELF_INTERSECTION, "Edge %d crosses edge %d", ap, bp, is_polygon)
            }
        }
        // Loops are not allowed to have duplicate vertices, and separate loops
        // are not allowed to share edges or cross at vertices.  We only need to
        // check a given vertex once, so we also require that the two edges have
        // the same end vertex.
        if (a.v1 != b.v1) return S2Error(S2Error.OK)
        if (ap.chainId == bp.chainId) {
            return loopError(S2Error.DUPLICATE_VERTICES, "Edge %d has duplicate vertex with edge %d", ap, bp, is_polygon)
        }
        val a_len = shape.chain(ap.chainId).length;
        val b_len = shape.chain(bp.chainId).length;
        val a_next = if(ap.offset + 1 == a_len) 0 else (ap.offset + 1)
        val b_next = if(bp.offset + 1 == b_len) 0 else (bp.offset + 1)
        val a2 = shape.chainEdge(ap.chainId, a_next).v1
        val b2 = shape.chainEdge(bp.chainId, b_next).v1
        if (a.v0 == b.v0 || a.v0 == b2) {
            // The second edge index is sometimes off by one, hence "near".
            return S2Error(code = S2Error.POLYGON_LOOPS_SHARE_EDGE, text = "Loop %d edge %d has duplicate near loop %d edge %d".format(
                ap.chainId, ap.offset, bp.chainId, bp.offset))
        }
        // Since S2ShapeIndex loops are oriented such that the polygon interior is
        // always on the left, we need to handle the case where one wedge contains
        // the complement of the other wedge.  This is not specifically detected by
        // GetWedgeRelation, so there are two cases to check for.
        //
        // Note that we don't need to maintain any state regarding loop crossings
        // because duplicate edges are detected and rejected above.
        if (S2WedgeRelations.getWedgeRelation(a.v0, a.v1, a2, b.v0, b2) == S2WedgeRelations.WedgeRelation.WEDGE_PROPERLY_OVERLAPS &&
                S2WedgeRelations.getWedgeRelation(a.v0, a.v1, a2, b2, b.v0) == S2WedgeRelations.WedgeRelation.WEDGE_PROPERLY_OVERLAPS) {
            return S2Error(code = S2Error.POLYGON_LOOPS_CROSS, "Loop %d edge %d crosses loop %d edge %d".format(
                    ap.chainId, ap.offset, bp.chainId, bp.offset))
        }
        return S2Error(S2Error.OK)
    }
}

typealias ShapeEdgeVector = MutableList<ShapeEdge>

// A function that is called with pairs of crossing edges.  The function may
// return false in order to request that the algorithm should be terminated,
// i.e. no further crossings are needed.
//
// "is_interior" indicates that the crossing is at a point interior to both
// edges (i.e., not at a vertex).  (The calling function already has this
// information and it is moderately expensive to recompute.)
@FunctionalInterface
interface EdgePairVisitor {

    fun visit(a: ShapeEdge, b: ShapeEdge, is_interior: Boolean): Boolean

}