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

    /*
    // Allows iterating over the indexed shapes using range-based for loops:
    //
    //   for (S2Shape* shape : index) { ... }
    //
    // CAVEAT: Returns nullptr for shapes that have been removed from the index.
    class ShapeIterator
        : public std::iterator<std::forward_iterator_tag, S2Shape*>
    {
        public:
        ShapeIterator() = default
        S2Shape * operator * () const
        ShapeIterator& operator++()
        ShapeIterator operator ++(int)

        // REQUIRES: "it" and *this must reference the same S2ShapeIndex.
        bool operator ==(ShapeIterator it) const

        // REQUIRES: "it" and *this must reference the same S2ShapeIndex.
        bool operator !=(ShapeIterator it) const

        private:
        friend class S2ShapeIndex
        ShapeIterator(const S2ShapeIndex * index, int shape_id)
        : index_(index), shape_id_(shape_id) {}

        const S2ShapeIndex * index_ = nullptr
        int shape_id_ = 0
    }
    ShapeIterator begin() const 
    ShapeIterator end() const 
     */

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
    // A random access iterator that provides low-level access to the cells of
    // the index.  Cells are sorted in increasing order of S2CellId.
    class Iterator {
        public :
        // Default constructor; must be followed by a call to Init().
        Iterator() : iter_(nullptr)
        {}

        // Constructs an iterator positioned as specified.  By default iterators
        // are unpositioned, since this avoids an extra seek in this situation
        // where one of the seek methods (such as Locate) is immediately called.
        //
        // If you want to position the iterator at the beginning, e.g. in order to
        // loop through the entire index, do this instead:
        //
        //   for (S2ShapeIndex::Iterator it(&index, S2ShapeIndex::BEGIN)
        //        !it.done(); it.Next()) { ... }
        explicit Iterator(const S2ShapeIndex* index,
        InitialPosition pos = UNPOSITIONED)
        : iter_(index->NewIterator(pos))
        {}

        // Initializes an iterator for the given S2ShapeIndex.  This method may
        // also be called in order to restore an iterator to a valid state after
        // the underlying index has been updated (although it is usually easier
        // just to declare a new iterator whenever required, since iterator
        // construction is cheap).
        void Init(const S2ShapeIndex* index,
        InitialPosition pos = UNPOSITIONED)
        {
            iter_ = index->NewIterator(pos)
        }

        // Iterators are copyable and movable.
        Iterator(const Iterator&)
        Iterator& operator =(const Iterator&)
        Iterator(Iterator&&)
        Iterator& operator =(Iterator&&)

        // Returns the S2CellId of the current index cell.  If done() is true,
        // returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
        S2CellId id()
        const { return iter_->id(); }

        // Returns the center point of the cell.
        // REQUIRES: !done()
        S2Point center()
        const { return id().ToPoint(); }

        // Returns a reference to the contents of the current index cell.
        // REQUIRES: !done()
        const S2ShapeIndexCell& cell()
        const { return iter_->cell(); }

        // Returns true if the iterator is positioned past the last index cell.
        bool done()
        const { return iter_->done(); }

        // Positions the iterator at the first index cell (if any).
        void Begin()
        { iter_->Begin(); }

        // Positions the iterator past the last index cell.
        void Finish()
        { iter_->Finish(); }

        // Positions the iterator at the next index cell.
        // REQUIRES: !done()
        void Next()
        { iter_->Next(); }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        bool Prev()
        { return iter_->Prev(); }

        // Positions the iterator at the first cell with id() >= target, or at the
        // end of the index if no such cell exists.
        void Seek(S2CellId target)
        { iter_->Seek(target); }

        // Positions the iterator at the cell containing "target".  If no such cell
        // exists, returns false and leaves the iterator positioned arbitrarily.
        // The returned index cell is guaranteed to contain all edges that might
        // intersect the line segment between "target" and the cell center.
        bool Locate(const S2Point& target)
        {
            return IteratorBase::LocateImpl(target, this)
        }

        // Let T be the target S2CellId.  If T is contained by some index cell I
        // (including equality), this method positions the iterator at I and
        // returns INDEXED.  Otherwise if T contains one or more (smaller) index
        // cells, it positions the iterator at the first such cell I and returns
        // SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
        // positioned arbitrarily.
        CellRelation Locate(S2CellId target)
        {
            return IteratorBase::LocateImpl(target, this)
        }

        private :
        // Although S2ShapeIndex::Iterator can be used to iterate over any
        // index subtype, it is more efficient to use the subtype's iterator when
        // the subtype is known at compile time.  For example, MutableS2ShapeIndex
        // should use a MutableS2ShapeIndex::Iterator.
        //
        // The following declarations prevent accidental use of
        // S2ShapeIndex::Iterator when the actual subtype is known.  (If you
        // really want to do this, you can down_cast the index argument to
        // S2ShapeIndex.)
        template <
        class T>
        explicit Iterator(const T* index, InitialPosition pos = UNPOSITIONED)
        {}

        template <
        class T>
        void Init(const T* index, InitialPosition pos = UNPOSITIONED)
        {}

        std::unique_ptr<IteratorBase> iter_
    }

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

    protected :
    // Each subtype of S2ShapeIndex should define an Iterator type derived
    // from the following base class.
    class IteratorBase {
        public :
        virtual ~IteratorBase()
        {}

        IteratorBase(const IteratorBase&)
        IteratorBase& operator =(const IteratorBase&)

        // Returns the S2CellId of the current index cell.  If done() is true,
        // returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
        S2CellId id() const 

        // Returns the center point of the cell.
        // REQUIRES: !done()
        S2Point center() const 

        // Returns a reference to the contents of the current index cell.
        // REQUIRES: !done()
        const S2ShapeIndexCell& cell() const 

        // Returns true if the iterator is positioned past the last index cell.
        bool done() const 

        // Positions the iterator at the first index cell (if any).
        virtual void Begin() = 0

        // Positions the iterator past the last index cell.
        virtual void Finish() = 0

        // Positions the iterator at the next index cell.
        // REQUIRES: !done()
        virtual void Next() = 0

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        virtual bool Prev() = 0

        // Positions the iterator at the first cell with id() >= target, or at the
        // end of the index if no such cell exists.
        virtual void Seek(S2CellId target) = 0

        // Positions the iterator at the cell containing "target".  If no such cell
        // exists, returns false and leaves the iterator positioned arbitrarily.
        // The returned index cell is guaranteed to contain all edges that might
        // intersect the line segment between "target" and the cell center.
        virtual bool Locate(const S2Point& target) = 0

        // Let T be the target S2CellId.  If T is contained by some index cell I
        // (including equality), this method positions the iterator at I and
        // returns INDEXED.  Otherwise if T contains one or more (smaller) index
        // cells, it positions the iterator at the first such cell I and returns
        // SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
        // positioned arbitrarily.
        virtual CellRelation Locate(S2CellId target) = 0

        protected :
        IteratorBase() : id_(S2CellId::Sentinel()), cell_(nullptr)
        {}

        // Sets the iterator state.  "cell" typically points to the cell contents,
        // but may also be given as "nullptr" in order to implement decoding on
        // demand.  In that situation, the first that the client attempts to
        // access the cell contents, the GetCell() method is called and "cell_" is
        // updated in a thread-safe way.
        void set_state(S2CellId id, const S2ShapeIndexCell* cell)

        // Sets the iterator state so that done() is true.
        void set_finished()

        // Returns the current contents of the "cell_" field, which may be nullptr
        // if the cell contents have not been decoded yet.
        const S2ShapeIndexCell* raw_cell() const 

        // This method is called to decode the contents of the current cell, if
        // set_state() was previously called with a nullptr "cell" argument.  This
        // allows decoding on demand for subtypes that keep the cell contents in
        // an encoded state.  It does not need to be implemented at all if
        // set_state() is always called with (cell != nullptr).
        //
        // REQUIRES: This method is thread-safe.
        // REQUIRES: Multiple calls to this method return the same value.
        virtual const S2ShapeIndexCell* GetCell() const = 0

        // Returns an exact copy of this iterator.
        virtual std::unique_ptr<IteratorBase> Clone() const = 0

        // Makes a copy of the given source iterator.
        // REQUIRES: "other" has the same concrete type as "this".
        virtual void Copy(const IteratorBase& other) = 0

        // The default implementation of Locate(S2Point).  It is instantiated by
        // each subtype in order to (1) minimize the number of virtual method
        // calls (since subtypes are typically "final") and (2) ensure that the
        // correct versions of non-virtual methods such as cell() are called.
        template <
        class Iter>
        static bool LocateImpl(const S2Point& target, Iter* it)

        // The default implementation of Locate(S2CellId) (see comments above).
        template <
        class Iter>
        static CellRelation LocateImpl(S2CellId target, Iter* it)

        private :
        friend
        class Iterator

        // This method is "const" because it is used internally by "const" methods
        // in order to implement decoding on demand.
        void set_cell(const S2ShapeIndexCell* cell) const 

        S2CellId id_
        mutable std::atomic<const S2ShapeIndexCell*> cell_
    }

    
 */

    // Returns a new iterator positioned as specified.
    abstract fun iterator(pos: InitialPosition = InitialPosition.UNPOSITIONED): S2ShapeIndexIteratorBase

}
/*
//////////////////   Implementation details follow   ////////////////////


inline S2Shape* S2ShapeIndex::ShapeIterator::operator *() const {
    return index_->shape(shape_id_)
}

inline S2ShapeIndex::ShapeIterator& S2ShapeIndex::ShapeIterator::operator ++() {
    ++shape_id_
    return * this
}

inline S2ShapeIndex::ShapeIterator S2ShapeIndex::ShapeIterator::operator ++(
int) {
    return ShapeIterator(index_, shape_id_++)
}

inline bool S2ShapeIndex::ShapeIterator::operator ==(ShapeIterator it) const {
    S2_DCHECK_EQ(index_, it.index_)
    return shape_id_ == it.shape_id_
}

inline bool S2ShapeIndex::ShapeIterator::operator !=(ShapeIterator it) const {
    S2_DCHECK_EQ(index_, it.index_)
    return shape_id_ != it.shape_id_
}

inline S2ShapeIndex::ShapeIterator S2ShapeIndex::begin() const {
    return ShapeIterator(this, 0)
}

inline S2ShapeIndex::ShapeIterator S2ShapeIndex::end() const {
    return ShapeIterator(this, num_shape_ids())
}

inline S2ShapeIndex::IteratorBase::IteratorBase(const IteratorBase& other)
: id_(other.id_), cell_(other.raw_cell()) {
}

inline S2ShapeIndex::IteratorBase&
S2ShapeIndex::IteratorBase::operator =(const IteratorBase& other) {
    id_ = other.id_
    set_cell(other.raw_cell())
    return * this
}

inline S2CellId S2ShapeIndex::IteratorBase::id() const {
    return id_
}

inline const S2ShapeIndexCell& S2ShapeIndex::IteratorBase::cell() const {
    // Like other const methods, this method is thread-safe provided that it
    // does not overlap with calls to non-const methods.
    S2_DCHECK(!done())
    auto cell = raw_cell ()
    if (cell == nullptr) {
        cell = GetCell()
        set_cell(cell)
    }
    return * cell
}

inline bool S2ShapeIndex::IteratorBase::done() const {
    return id_ == S2CellId::Sentinel()
}

inline S2Point S2ShapeIndex::IteratorBase::center() const {
    S2_DCHECK(!done())
    return id().ToPoint()
}

inline void S2ShapeIndex::IteratorBase::set_state(
S2CellId id, const S2ShapeIndexCell* cell) {
    id_ = id
    set_cell(cell)
}

inline void S2ShapeIndex::IteratorBase::set_finished() {
    id_ = S2CellId::Sentinel()
    set_cell(nullptr)
}

inline const S2ShapeIndexCell* S2ShapeIndex::IteratorBase::raw_cell()
const {
    return cell_.load(std::memory_order_relaxed)
}

inline void S2ShapeIndex::IteratorBase::set_cell(
const S2ShapeIndexCell* cell) const {
    cell_.store(cell, std::memory_order_relaxed)
}

template <class Iter>
inline bool S2ShapeIndex::IteratorBase::LocateImpl(
const S2Point& target_point, Iter* it) {
    // Let I = cell_map_->lower_bound(T), where T is the leaf cell containing
    // "target_point".  Then if T is contained by an index cell, then the
    // containing cell is either I or I'.  We test for containment by comparing
    // the ranges of leaf cells spanned by T, I, and I'.

    S2CellId target (target_point)
    it->Seek(target)
    if (!it->done() && it->id().range_min() <= target) return true
    if (it->Prev() && it->id().range_max() >= target) return true
    return false
}

template <class Iter>
inline S2ShapeIndex::CellRelation
S2ShapeIndex::IteratorBase::LocateImpl(S2CellId target, Iter* it) {
    // Let T be the target, let I = cell_map_->lower_bound(T.range_min()), and
    // let I' be the predecessor of I.  If T contains any index cells, then T
    // contains I.  Similarly, if T is contained by an index cell, then the
    // containing cell is either I or I'.  We test for containment by comparing
    // the ranges of leaf cells spanned by T, I, and I'.

    it->Seek(target.range_min())
    if (!it->done()) {
        if (it->id() >= target && it->id().range_min() <= target) return INDEXED
        if (it->id() <= target.range_max()) return SUBDIVIDED
    }
    if (it->Prev() && it->id().range_max() >= target) return INDEXED
    return DISJOINT
}

inline S2ShapeIndex::Iterator::Iterator(const Iterator& other)
: iter_(other.iter_->Clone()) {
}

inline S2ShapeIndex::Iterator& S2ShapeIndex::Iterator::operator =(
const Iterator& other) {
    iter_->Copy(*other.iter_)
    return * this
}

inline S2ShapeIndex::Iterator::Iterator(Iterator&& other)
: iter_(std::move(other.iter_)) {
}

inline S2ShapeIndex::Iterator& S2ShapeIndex::Iterator::operator =(
Iterator&& other) {
    iter_ = std::move(other.iter_)
    return * this
}

#endif  // S2_S2SHAPE_INDEX_H_
*/