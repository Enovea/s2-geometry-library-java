package dilivia.s2

class S2CellIdIterator(val parent: S2CellId) : Iterator<S2CellId> {

    private var currentId: S2CellId = parent.rangeMin()
    private val max = parent.rangeMax()

    override fun hasNext(): Boolean = currentId < max

    override fun next(): S2CellId {
        val next = currentId
        currentId = currentId.next()
        return next
    }
}