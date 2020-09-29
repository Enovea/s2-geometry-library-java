package dilivia.s2.shape

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser

class S2LaxPolylineShapeTest : S2GeometryTestCase() {

    fun testNoVertices() {
        val shape = S2LaxPolylineShape()
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(1, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testOneVertex() {
        val vertices = listOf(S2Point(1, 0, 0))
        val shape = S2LaxPolylineShape(vertices)
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(1, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
    }

    fun testEdgeAccess() {
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1")
        val shape = S2LaxPolylineShape(vertices)
        assertEquals(2, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertEquals(0, shape.chain(0).start)
        assertEquals(2, shape.chain(0).length)
        assertEquals(1, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        val edge0 = shape.edge(0)
        assertEquals(vertices[0], edge0.v0)
        assertEquals(vertices[1], edge0.v1)
        val edge1 = shape.edge(1)
        assertEquals(vertices[1], edge1.v0)
        assertEquals(vertices[2], edge1.v1)
    }

}