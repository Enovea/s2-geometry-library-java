package dilivia.s2.shape

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2Random

class S2PointVectorShapeTest : S2GeometryTestCase() {

    fun testEmpty() {
        val shape = S2PointVectorShape()
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(0, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testConstructionAndAccess() {
        val points = mutableListOf<S2Point>()
        S2Random.reset(0);
        val kNumPoints = 100;
        repeat(kNumPoints) {
            points.add(S2Random.randomPoint())
        }
        val shape = S2PointVectorShape(points = points)

        assertEquals(kNumPoints, shape.numEdges)
        assertEquals(kNumPoints, shape.numChains)
        assertEquals(0, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        S2Random.reset(0)
        repeat(kNumPoints) { i ->
            assertEquals(i, shape.chain(i).start)
            assertEquals(1, shape.chain(i).length)
            val edge = shape.edge(i)
            val pt = S2Random.randomPoint()
            assertEquals(pt, edge.v0)
            assertEquals(pt, edge.v1)
        }
    }

}
