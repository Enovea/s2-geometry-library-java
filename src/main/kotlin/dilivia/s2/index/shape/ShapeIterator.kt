package dilivia.s2.index.shape

import dilivia.s2.shape.S2Shape

/**
 * Allows iterating over the indexed shapes.
 * CAVEAT: Returns null for shapes that have been removed from the index.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
class ShapeIterator(val index: S2ShapeIndex) {

    /** Current shape id. */
    private var shapeId: Int = 0

    internal constructor(index: S2ShapeIndex, shapeId: Int): this(index) {
        this.shapeId = shapeId
    }

    /**
     * Gets the shape at the current position, null if the shape has been removed.
     *
     * @return The current shape.
     */
    fun shape(): S2Shape? = index.shape(shapeId)

    /**
     * Increments the iterator.
     *
     * @return A new iterator instance that points to the next shape id.
     */
    operator fun inc(): ShapeIterator = ShapeIterator(index = index, shapeId = shapeId + 1)

    fun next(): ShapeIterator { ++shapeId; return this }

    /**
     * Decrements the iterator.
     *
     * @return A new iterator instance that points to the previous shape id.
     */
    operator fun dec(): ShapeIterator = ShapeIterator(index = index, shapeId = shapeId - 1)

    fun previous(): ShapeIterator { --shapeId; return this }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is ShapeIterator) return false

        if (index != other.index) return false
        if (shapeId != other.shapeId) return false

        return true
    }

    override fun hashCode(): Int {
        var result = index.hashCode()
        result = 31 * result + shapeId
        return result
    }

}
