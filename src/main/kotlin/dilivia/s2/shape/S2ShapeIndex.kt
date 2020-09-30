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

//
// S2ShapeIndex is an abstract base class for indexing polygonal geometry in
// memory.  The main documentation is with the class definition below.
// (Some helper classes are defined first.)


// S2ShapeIndex is an abstract base class for indexing polygonal geometry in
// memory.  The objects in the index are known as "shapes", and may consist of
// points, polylines, and/or polygons, possibly overlapping.  The index makes
// it very fast to answer queries such as finding nearby shapes, measuring
// distances, testing for intersection and containment, etc.
//
// Each object in the index implements the S2Shape interface.  An S2Shape is a
// collection of edges that optionally defines an interior.  The edges do not
// need to be connected, so for example an S2Shape can represent a polygon
// with multiple shells and/or holes, or a set of polylines, or a set of
// points.  All geometry within a single S2Shape must have the same dimension,
// so for example if you want to create an S2ShapeIndex containing a polyline
// and 10 points, then you will need at least two different S2Shape objects.
//
// The most important type of S2ShapeIndex is MutableS2ShapeIndex, which
// allows you to build an index incrementally by adding or removing shapes.
// Soon there will also be an EncodedS2ShapeIndex type that makes it possible
// to keep the index data in encoded form.  Code that only needs read-only
// ("const") access to an index should use the S2ShapeIndex base class as the
// parameter type, so that it will work with any S2ShapeIndex subtype.  For
// example:
//
//   void DoSomething(const S2ShapeIndex& index) {
//     ... works with MutableS2ShapeIndex or EncodedS2ShapeIndex ...
//   }
//
// There are a number of built-in classes that work with S2ShapeIndex objects.
// Generally these classes accept any collection of geometry that can be
// represented by an S2ShapeIndex, i.e. any combination of points, polylines,
// and polygons.  Such classes include:
//
// - S2ContainsPointQuery: returns the shape(s) that contain a given point.
//
// - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge,
//                       S2CellId, or S2ShapeIndex.
//
// - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
//
// - S2BooleanOperation: computes boolean operations such as union,
//                       and boolean predicates such as containment.
//
// - S2ShapeIndexRegion: computes approximations for a collection of geometry.
//
// - S2ShapeIndexBufferedRegion: computes approximations that have been
//                               expanded by a given radius.
//
// Here is an example showing how to index a set of polygons and then
// determine which polygon(s) contain each of a set of query points:
//
//   void TestContainment(const vector<S2Point>& points,
//                        const vector<S2Polygon*>& polygons) {
//     MutableS2ShapeIndex index
//     for (auto polygon : polygons) {
//       index.Add(absl::make_unique<S2Polygon::Shape>(polygon))
//     }
//     auto query = MakeS2ContainsPointQuery(&index)
//     for (const auto& point : points) {
//       for (S2Shape* shape : query.GetContainingShapes(point)) {
//         S2Polygon* polygon = polygons[shape->id()]
//         ... do something with (point, polygon) ...
//       }
//     }
//   }
//
// This example uses S2Polygon::Shape, which is one example of an S2Shape
// object.  S2Polyline and S2Loop also have nested Shape classes, and there are
// additional S2Shape types defined in *_shape.h.
//
// Internally, an S2ShapeIndex is essentially a map from S2CellIds to the set
// of shapes that intersect each S2CellId.  It is adaptively refined to ensure
// that no cell contains more than a small number of edges.
//
// In addition to implementing a shared set of virtual methods, all
// S2ShapeIndex subtypes define an Iterator type with the same API.  This
// makes it easy to convert code that uses a particular S2ShapeIndex subtype
// to instead use the abstract base class (or vice versa).  You can also
// choose to avoid the overhead of virtual method calls by making the
// S2ShapeIndex type a template argument, like this:
//
//   template <class IndexType>
//   void DoSomething(const IndexType& index) {
//     for (typename IndexType::Iterator it(&index, S2ShapeIndex::BEGIN)
//          !it.done(); it.Next()) {
//       ...
//     }
//   }
//
// Subtypes provided by the S2 library have the same thread-safety properties
// as std::vector.  That is, const methods may be called concurrently from
// multiple threads, and non-const methods require exclusive access to the
// S2ShapeIndex.
abstract class S2ShapeIndex {

    // Returns the number of distinct shape ids in the index.  This is the same
    // as the number of shapes provided that no shapes have ever been removed.
    // (Shape ids are never reused.)
    abstract fun numShapeIds(): Int

    // Returns a pointer to the shape with the given id, or nullptr if the shape
    // has been removed from the index.
    abstract fun shape(id: Int): S2Shape?

    // The possible relationships between a "target" cell and the cells of the
    // S2ShapeIndex.  If the target is an index cell or is contained by an index
    // cell, it is "INDEXED".  If the target is subdivided into one or more
    // index cells, it is "SUBDIVIDED".  Otherwise it is "DISJOINT".
    enum class CellRelation {
        INDEXED,       // Target is contained by an index cell
        SUBDIVIDED,    // Target is subdivided into one or more index cells
        DISJOINT       // Target does not intersect any index cells
    }

    // When passed to an Iterator constructor, specifies whether the iterator
    // should be positioned at the beginning of the index (BEGIN), the end of
    // the index (END), or arbitrarily (UNPOSITIONED).  By default iterators are
    // unpositioned, since this avoids an extra seek in this situation where one
    // of the seek methods (such as Locate) is immediately called.
    enum class InitialPosition { BEGIN, END, UNPOSITIONED }
/*

    // ShapeFactory is an interface for decoding vectors of S2Shapes.  It allows
    // random access to the shapes in order to support lazy decoding.  See
    // s2shapeutil_coding.h for useful subtypes.
    class ShapeFactory {
        public :
        virtual ~ShapeFactory()
        {}

        // Returns the number of S2Shapes in the vector.
        virtual int size() const = 0

        // Returns the S2Shape object corresponding to the given "shape_id".
        // Returns nullptr if a shape cannot be decoded or a shape is missing
        // (e.g., because MutableS2ShapeIndex::Release() was called).
        virtual std::unique_ptr<S2Shape> operator [](int shape_id) const = 0

        // Returns a deep copy of this ShapeFactory.
        virtual std::unique_ptr<ShapeFactory> Clone() const = 0
    }
*/

    // Returns a new iterator positioned as specified.
    abstract fun iterator(pos: InitialPosition = InitialPosition.UNPOSITIONED): S2ShapeIndexIteratorBase

    abstract fun begin(): ListIterator<S2Shape?>
    abstract fun end(): ListIterator<S2Shape?>
}