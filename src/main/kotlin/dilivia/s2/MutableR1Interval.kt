package dilivia.s2

import com.google.common.geometry.R2Vector
import kotlin.math.max
import kotlin.math.min

class MutableR1Interval(override var lo: Double, override var hi: Double): R1Interval(lo, hi) {

    constructor(lo: Int, hi: Int): this(lo.toDouble(), hi.toDouble())

    operator fun set(i: Int, v: Double) {
        when(i) {
            0 -> lo = v
            1 -> hi = v
            else -> throw ArrayIndexOutOfBoundsException("Index $i is out of bounds 0..1")
        }
    }

    override var bounds: R2Vector
        get() = super.bounds
        set(value) {
            lo = value[0]
            hi = value[1]
        }

    override fun addPoint(p: Double): R1Interval {
        when {
            isEmpty -> {
                lo = p
                hi = p
            }
            p < lo -> lo = p
            p > hi -> hi = p
        }
        return this
    }

    override fun addInterval(y: R1Interval): R1Interval {
        if (y.isNotEmpty) {
            if (isEmpty) {
                lo = y.lo
                hi = y.hi
            }
            else {
                lo = min(lo, y.lo)
                hi = max(hi, y.hi)
            }
        }
        return this
    }

    fun expands(radius: Double): R1Interval {
        // assert (radius >= 0);
        if (isNotEmpty) {
            lo -= radius
            hi += radius
        }
        return this
    }

    fun setUnion(y: R1Interval): R1Interval {
        if (isEmpty) {
            lo = y.lo
            hi = y.hi
        }
        else if (y.isNotEmpty) {
            lo = min(lo, y.lo)
            hi = max(hi, y.hi)
        }
        return this
    }

    fun setIntersection(y: R1Interval): R1Interval {
        lo = max(lo, y.lo)
        hi = min(hi, y.hi)
        return this
    }

}