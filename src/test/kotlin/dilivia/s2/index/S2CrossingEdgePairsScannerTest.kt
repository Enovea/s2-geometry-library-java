package dilivia.s2.index

import dilivia.s2.S2EdgeCrossings
import dilivia.s2.S2Error
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2LatLng
import dilivia.s2.shape.S2EdgeVectorShape
import dilivia.s2.shape.ShapeEdge
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging

// A set of edge pairs within an S2ShapeIndex.
typealias EdgePairVector = Set<Pair<ShapeEdgeId, ShapeEdgeId>>

class S2CrossingEdgePairsScannerTest : S2GeometryTestCase() {

    fun getCrossings(index: S2ShapeIndex, type: CrossingType): EdgePairVector {
        val edge_pairs = HashSet<Pair<ShapeEdgeId, ShapeEdgeId>>()
        S2CrossingEdgePairsScanner.visitCrossingEdgePairs(
                index, type, object : EdgePairVisitor {
            override fun visit(a: ShapeEdge, b: ShapeEdge, is_interior: Boolean): Boolean {
                edge_pairs.add(Pair(a.id, b.id))
                return true;  // Continue visiting.
            }
        })
        return edge_pairs
    }

    fun getCrossingEdgePairsBruteForce(index: S2ShapeIndex, type: CrossingType): EdgePairVector {
        val result = HashSet<Pair<ShapeEdgeId, ShapeEdgeId>>()
        val min_sign = if (type == CrossingType.ALL) 0 else 1
        val a_iter = EdgeIterator(index)
        while (!a_iter.done()) {
            val a = a_iter.edge()
            val b_iter = EdgeIterator(a_iter)
            b_iter.next()
            while (!b_iter.done()) {
                val b = b_iter.edge()
                if (S2EdgeCrossings.crossingSign(a.v0, a.v1, b.v0, b.v1) >= min_sign) {
                    result.add(Pair(a_iter.shape_edge_id(), b_iter.shape_edge_id()))
                }
                b_iter.next()
            }
            a_iter.next()
        }
        return result
    }

    fun testGetCrossingEdgePairs(index: S2ShapeIndex, type: CrossingType) {
        val expected = getCrossingEdgePairsBruteForce(index, type)
        val actual = getCrossings(index, type)
        if (actual != expected) {
            var message = """
            |Unexpected edge pairs; see details below.
            |Expected number of edge pairs: ${expected.size}
            |Actual number of edge pairs: ${actual.size}
            |
          """.trimMargin()
            expected.filter { edgePair -> !actual.contains(edgePair) }.forEach { edgePair -> message += "Missing value: $edgePair\n" }
            actual.filter { edgePair -> !expected.contains(edgePair) }.forEach { edgePair -> message += "Extra value: $edgePair\n" }
            fail(message)
        }
    }

    fun testGetCrossingEdgePairsNoIntersections() {
        val index = MutableS2ShapeIndex()
        testGetCrossingEdgePairs(index, CrossingType.ALL)
        testGetCrossingEdgePairs(index, CrossingType.INTERIOR)
    }

    fun testGetCrossingEdgePairsEdgeGrid() {
        val kGridSize = 10;  // (kGridSize + 1) * (kGridSize + 1) crossings
        val index = MutableS2ShapeIndex()
        val shape = S2EdgeVectorShape()
        for (i in 0..kGridSize) {
            shape.add(S2LatLng.fromDegrees(0, i).toPoint(), S2LatLng.fromDegrees(kGridSize, i).toPoint())
            shape.add(S2LatLng.fromDegrees(i, 0).toPoint(), S2LatLng.fromDegrees(i, kGridSize).toPoint())
        }
        index.add(shape)
        testGetCrossingEdgePairs(index, CrossingType.ALL)
        testGetCrossingEdgePairs(index, CrossingType.INTERIOR)
    }
/*
    // This function recursively verifies that HasCrossing returns the given
    // result for all possible cyclic permutations of the loop vertices for the
    // given set of loops.
    fun testHasCrossingPermutations(loops: MutableList<S2Loop>, i: Int, has_crossing: Boolean) {
        if (i == loops.size) {
            val index = MutableS2ShapeIndex()
            val polygon = S2Polygon(loops, check = false)
            index.add(polygon)
            assertEquals(has_crossing, hasSelfIntersection(index))
        } else {
            val orig_loop = loops[i]
            for (j in 0 until orig_loop.numVertices()) {
                val vertices = mutableListOf<S2Point>()
                for (k in 0 until orig_loop.numVertices()) {
                    vertices.add(orig_loop.vertex(j + k))
                }
                loops[i] = S2Loop(vertices, check = false)
                testHasCrossingPermutations(loops, i + 1, has_crossing)
            }
            loops[i] = orig_loop
        }
    }

    // Given a string reprsenting a polygon, and a boolean indicating whether this
    // polygon has any self-intersections or loop crossings, verify that all
    // HasSelfIntersection returns the expected result for all possible cyclic
    // permutations of the loop vertices.
    fun testHasCrossing(polygon_str: String, has_crossing: Boolean) {
        // Set S2Debug::DISABLE to allow invalid polygons.
        unique_ptr<S2Polygon> polygon = s2textformat::MakePolygonOrDie(polygon_str, S2Debug::DISABLE);
        val loops = polygon->Release();
        testHasCrossingPermutations(loops, 0, has_crossing);
    }

    fun testFindSelfIntersectionBasic() {
        // Coordinates are (lat,lng), which can be visualized as (y,x).
        testHasCrossing("0:0, 0:1, 0:2, 1:2, 1:1, 1:0", false);
        testHasCrossing("0:0, 0:1, 0:2, 1:2, 0:1, 1:0", true);  // duplicate vertex
        testHasCrossing("0:0, 0:1, 1:0, 1:1", true);  // edge crossing
        testHasCrossing("0:0, 1:1, 0:1; 0:0, 1:1, 1:0", true);  // duplicate edge
        testHasCrossing("0:0, 1:1, 0:1; 1:1, 0:0, 1:0", true);  // reversed edge
        testHasCrossing("0:0, 0:2, 2:2, 2:0; 1:1, 0:2, 3:1, 2:0", true);  // vertex crossing
    }

 */

    companion object {

        private val logger = KotlinLogging.logger { }

        // Return true if any loop crosses any other loop (including vertex crossings
        // and duplicate edges), or any loop has a self-intersection (including
        // duplicate vertices).
        fun hasSelfIntersection(index: MutableS2ShapeIndex): Boolean {
            val error = S2CrossingEdgePairsScanner.findSelfIntersection(index)
            if (error.code != S2Error.OK) {
                logger.error { error }
                return true
            }
            return false
        }

    }
}
