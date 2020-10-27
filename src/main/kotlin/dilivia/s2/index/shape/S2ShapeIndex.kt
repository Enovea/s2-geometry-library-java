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
package dilivia.s2.index.shape

import dilivia.s2.shape.S2Shape
import kotlin.collections.Iterator

/**
 * S2ShapeIndex is an abstract base class for indexing polygonal geometry. The objects in the index are known as
 * "shapes", and may consist of points, polylines, and/or polygons, possibly overlapping.  The index makes it very fast
 * to answer queries such as finding nearby shapes, measuring distances, testing for intersection and containment, etc.
 *
 * Each object in the index implements the S2Shape interface. An S2Shape is a collection of edges that optionally
 * defines an interior. The edges do not need to be connected, so for example an S2Shape can represent a polygon
 * with multiple shells and/or holes, or a set of polylines, or a set of points.
 *
 * All geometry within a single S2Shape must have the same dimension, so for example if you want to create an
 * S2ShapeIndex containing a polyline and 10 points, then you will need at least two different S2Shape objects.
 *
 * The most important type of S2ShapeIndex is MutableS2ShapeIndex, which allows you to build an index incrementally by
 * adding or removing shapes.
 *
 * There are a number of built-in classes that work with S2ShapeIndex objects. Generally these classes accept any
 * collection of geometry that can be represented by an S2ShapeIndex, i.e. any combination of points, polylines,
 * and polygons.  Such classes include:
 *
 * - S2ContainsPointQuery: returns the shape(s) that contain a given point.
 * - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge, S2CellId, or S2ShapeIndex.
 * - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
 * - S2BooleanOperation: computes boolean operations such as union, and boolean predicates such as containment.
 * - S2ShapeIndexRegion: computes approximations for a collection of geometry.
 * - S2ShapeIndexBufferedRegion: computes approximations that have been expanded by a given radius.
 *
 * Here is an example showing how to index a set of polygons and then
 * determine which polygon(s) contain each of a set of query points:
 *
 * <pre>
 *   fun testContainment(points: List<S2Point>, polygons: List<S2Polygon>) {
 *     val index = MutableS2ShapeIndex()
 *     for (polygon in polygons) {
 *       index.add(polygon)
 *     }
 *     val query = makeS2ContainsPointQuery(index)
 *     for (point in points) {
 *       for (shape : query.getContainingShapes(point)) {
 *         val polygon = polygons[shape.id()]
 *         ... do something with (point, polygon) ...
 *       }
 *     }
 *   }
 * </pre>
 *
 * This example uses S2Polygon.Shape, which is one example of an S2Shape object. S2Polyline and S2Loop also have
 * nested Shape classes, and there are additional S2Shape types defined in *Shape.kt
 *
 * Internally, an S2ShapeIndex is essentially a map from S2CellIds to the set of shapes that intersect each S2CellId.
 * It is adaptively refined to ensure that no cell contains more than a small number of edges.
 *
 * In addition to implementing a set of abstract methods, all S2ShapeIndex subtypes define an Iterator type with the
 * same API. This makes it easy to convert code that uses a particular S2ShapeIndex subtype to instead use the abstract
 * base class (or vice versa).
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
abstract class S2ShapeIndex: Iterable<S2Shape?> {

    /**
     * Gets the number of distinct shape ids in the index. This is the same as the number of shapes provided that no
     * shapes have ever been removed. (Shape ids are never reused.)
     *
     * @return The num of used shape ids.
     */
    abstract fun numShapeIds(): Int

    /**
     * Gets the shape with the given id, or null if the shape has been removed from the index.
     *
     * @param id A shape id.
     * @return The shape with the given id or null if it has been removed.
     */
    abstract fun shape(id: Int): S2Shape?

    fun begin(): ShapeIterator = ShapeIterator(this, 0)

    fun end(): ShapeIterator = ShapeIterator(this, numShapeIds())

    /**
     * Gets an iterator over the indexed shape. Example usage:
     * <pre>
     *     for(shape in index.shapeIterator()) {
     *          ...
     *     }
     * </pre>
     */
    override fun iterator(): Iterator<S2Shape?> = object : Iterator<S2Shape?> {
        private var shapeId: Int = 0
        override fun hasNext(): Boolean = shapeId < numShapeIds()
        override fun next(): S2Shape? = shape(shapeId++)
    }

    fun cellIterator(pos: InitialPosition = InitialPosition.UNPOSITIONED): S2ShapeIndexCellIterator = S2ShapeIndexCellIterator(index = this, pos = pos)

    /**
     * Gets a new iterator positioned as specified.
     *
     * @param pos The initial position of the iterator.
     * @return An iterator instance.
     */
    internal abstract fun newIterator(pos: InitialPosition): IteratorBase


    fun toDebugString(): String {
        var cellMap = ""
        val iterator = newIterator(pos = InitialPosition.BEGIN)
        while (!iterator.done()) {
            val cell = iterator.cell()
            val id = iterator.id()
            cellMap += "\n$id => ${cell.toDebugString(separator = "", prefix = "\n - ")}"
            iterator.next()
        }
        return """
            |${this.javaClass.simpleName} content
            |---------------------------------------------
            |shapes: ${this.iterator().asSequence().joinToString(separator = ",\n", prefix = "[\n", postfix = "\n]")}
            |CellMap:$cellMap
            |---------------------------------------------
        """.trimMargin()
    }


}
