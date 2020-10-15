package dilivia.s2.shape

import dilivia.s2.S1Angle
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2LatLng
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.index.MutableS2ShapeIndex
import dilivia.s2.index.S2ContainsPointQuery.Companion.makeS2ContainsPointQuery
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon

class S2LaxPolygonShapeTest : S2GeometryTestCase() {

    fun testEmptyPolygon() {
        val shape = S2LaxPolygonShape((S2Polygon()))
        assertEquals(0, shape.numLoops())
        assertEquals(0, shape.numVertices())
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(2, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testFullPolygon() {
        val shape = S2LaxPolygonShape(S2Polygon(S2Loop(listOf(S2Point(0, 0, -1)))))
        assertEquals(1, shape.numLoops())
        assertEquals(0, shape.numVertices())
        assertEquals(0, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertEquals(2, shape.dimension)
        assertFalse(shape.isEmpty())
        assertTrue(shape.isFull())
        assertTrue(shape.getReferencePoint().contained)
    }

    fun testSingleVertexPolygon() {
        // S2Polygon doesn't support single-vertex loops, so we need to construct
        // the S2LaxPolygonShape directly.
        val loops = mutableListOf<List<S2Point>>()
        loops.add(S2TextParser.parsePoints("0:0"))
        val shape = S2LaxPolygonShape(loops)
        assertEquals(1, shape.numLoops())
        assertEquals(1, shape.numVertices())
        assertEquals(1, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertEquals(0, shape.chain(0).start)
        assertEquals(1, shape.chain(0).length)
        val edge = shape.edge(0)
        assertEquals(loops[0][0], edge.v0)
        assertEquals(loops[0][0], edge.v1)
        assertTrue(edge == shape.chainEdge(0, 0))
        assertEquals(2, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testSingleLoopPolygon() {
        // Test S2Polygon constructor.
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1, 1:0")
        val shape = S2LaxPolygonShape(S2Polygon(S2Loop(vertices)))
        assertEquals(1, shape.numLoops())
        assertEquals(vertices.size, shape.numVertices())
        assertEquals(vertices.size, shape.numLoopVertices(0))
        assertEquals(vertices.size, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertEquals(0, shape.chain(0).start)
        assertEquals(vertices.size, shape.chain(0).length)
        for (i in 0 until vertices.size) {
            assertEquals(vertices[i], shape.loopVertex(0, i))
            val edge = shape.edge(i)
            assertEquals(vertices[i], edge.v0)
            assertEquals(vertices[(i + 1) % vertices.size], edge.v1)
            assertEquals(edge.v0, shape.chainEdge(0, i).v0)
            assertEquals(edge.v1, shape.chainEdge(0, i).v1)
        }
        assertEquals(2, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(S2Shape.containsBruteForce(shape, S2Point.origin()))
    }

    fun testMultiLoopPolygon() {
        // Test vector<vector<S2Point>> constructor.  Make sure that the loops are
        // oriented so that the interior of the polygon is always on the left.
        val loops = listOf(
                S2TextParser.parsePoints("0:0, 0:3, 3:3"),  // CCW
                S2TextParser.parsePoints("1:1, 2:2, 1:2")   // CW
        )
        val shape = S2LaxPolygonShape(loops)

        assertEquals(loops.size, shape.numLoops())
        var num_vertices = 0
        assertEquals(loops.size, shape.numChains)
        for (i in 0 until loops.size) {
            assertEquals(loops[i].size, shape.numLoopVertices(i))
            assertEquals(num_vertices, shape.chain(i).start)
            assertEquals(loops[i].size, shape.chain(i).length)
            for (j in loops[i].indices) {
                assertEquals(loops[i][j], shape.loopVertex(i, j))
                val edge = shape.edge(num_vertices + j)
                assertEquals(loops[i][j], edge.v0)
                assertEquals(loops[i][(j + 1) % loops[i].size], edge.v1)
            }
            num_vertices += loops[i].size
        }
        assertEquals(num_vertices, shape.numVertices())
        assertEquals(num_vertices, shape.numEdges)
        assertEquals(2, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(S2Shape.containsBruteForce(shape, S2Point.origin()))
    }

    fun testMultiLoopS2Polygon() {
      fail("TODO")
        // Verify that the orientation of loops representing holes is reversed when
        // converting from an S2Polygon to an S2LaxPolygonShape.
//        val polygon = MakePolygonOrDie ("0:0, 0:3, 3:3; 1:1, 1:2, 2:2")
//        S2LaxPolygonShape shape ( * polygon)
//        for (int i = 0; i < polygon->num_loops(); ++i) {
//            S2Loop * loop = polygon->loop(i)
//            for (int j = 0; j < loop->num_vertices(); ++j) {
//            assertEquals(loop->oriented_vertex(j),
//            shape.loop_vertex(i, j))
//        }
//        }
    }

    fun testManyLoopPolygon() {
        // Test a polygon with enough loops so that cumulative_vertices_ is used.
        val loops = mutableListOf<List<S2Point>>()
                repeat(100) { i ->
            val center = S2LatLng.fromDegrees(0, i).toPoint()
            loops.add(makeRegularPoints(center, S1Angle.degrees(0.1), S2Random.randomInt(3)))
        }
        val shape = S2LaxPolygonShape(loops)

        assertEquals(loops.size, shape.numLoops())
        var num_vertices = 0
        assertEquals(loops.size, shape.numChains)
        for (i in 0 until loops.size) {
            assertEquals(loops[i].size, shape.numLoopVertices(i))
            assertEquals(num_vertices, shape.chain(i).start)
            assertEquals(loops[i].size, shape.chain(i).length)
            for (j in loops[i].indices) {
            assertEquals(loops[i][j], shape.loopVertex(i, j))
            val edge = shape . edge (num_vertices + j)
            assertEquals(loops[i][j], edge.v0)
            assertEquals(loops[i][(j + 1) % loops[i].size], edge.v1)
        }
            num_vertices += loops[i].size
        }
        assertEquals(num_vertices, shape.numVertices())
        assertEquals(num_vertices, shape.numEdges)
    }

    fun testDegenerateLoops() {
        val loops = mutableListOf(
            S2TextParser.parsePoints("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
                S2TextParser.parsePoints("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
                S2TextParser.parsePoints("5:5, 6:6")
        )
        val shape = S2LaxPolygonShape(loops)
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testInvertedLoops() {
        val loops = mutableListOf(
                S2TextParser.parsePoints("1:2, 1:1, 2:2"),
                S2TextParser.parsePoints("3:4, 3:3, 4:4")
        )
        val shape = S2LaxPolygonShape(loops)
        assertTrue(S2Shape.containsBruteForce(shape, S2Point.origin()))
    }

    fun compareS2LoopToShape(loop: S2Loop, shape: S2Shape) {
        val index = MutableS2ShapeIndex()
                index.add(shape)
        val cap = loop .capBound
        val query = makeS2ContainsPointQuery(index)
        repeat(100) {
        val point = S2Random.samplePoint(cap)
        assertEquals(loop.contains(point), query.shapeContains(index.shape(0)!!, point))
    }
    }

    fun testCompareToS2Loop() {
        fail("TODO")
//        repeat(100) {
//            S2Testing::Fractal fractal fractal.set_max_level(S2Testing::rnd.Uniform(5))
//            fractal.set_fractal_dimension(1 + S2Testing::rnd.RandDouble())
//            S2Point center = S2Testing ::RandomPoint()
//            unique_ptr<S2Loop> loop (fractal.MakeLoop(
//                    S2Testing::GetRandomFrameAt(center), S1Angle::Degrees(5)))
//
//            // Compare S2Loop to S2LaxLoopShape.
//            CompareS2LoopToShape(*loop, make_unique<S2LaxLoopShape>(*loop))
//
//            // Compare S2Loop to S2LaxPolygonShape.
//            vector < S2LaxPolygonShape::Loop > loops(
//                    1, vector<S2Point>(& loop->vertex(0),
//            &loop->vertex(0)+loop->num_vertices()))
//            CompareS2LoopToShape(*loop, make_unique<S2LaxPolygonShape>(loops))
//        }
    }
}
