package dilivia.s2

@Strictfp
class MutableR2Rect(override var x: MutableR1Interval, override var y: MutableR1Interval) : R2Rect(x, y) {

    operator fun set(idx: Int, v: MutableR1Interval) {
        require(idx in 0..1)
        if (idx == 0) x = v else y = v
    }

    override fun addPoint(p: R2Point): MutableR2Rect {
        x.addPoint(p[0])
        y.addPoint(p[1])
        return this
    }

    // Expand the rectangle to include the given other rectangle.  This is the
    // same as replacing the rectangle by the union of the two rectangles, but
    // is somewhat more efficient.
    override fun addRect(other: R2Rect): MutableR2Rect {
        x.addInterval(other.x)
        y.addInterval(other.y)
        return this
    }

    fun expands(margin: R2Point): MutableR2Rect {
        x.expands(margin.x())
        y.expands(margin.y())
        if (x.isEmpty || y.isEmpty) {
            x.setEmpty()
            y.setEmpty()
        }
        return this
    }

    fun expands(margin: Double): MutableR2Rect = expands(R2Point(margin, margin))

    fun setUnion(other: R2Rect): MutableR2Rect {
        x.setUnion(other.x)
        y.setUnion(other.y)
        return this
    }

    fun setIntersection(other: R2Rect): MutableR2Rect {
        x.setIntersection(other.x)
        y.setIntersection(other.y)
        if (x.isEmpty || y.isEmpty) {
            x.setEmpty()
            y.setEmpty()
        }
        return this
    }
}