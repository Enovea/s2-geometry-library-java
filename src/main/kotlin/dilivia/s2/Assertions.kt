package dilivia.s2

import com.google.common.geometry.S2

object Assertions {

    val enabled: Boolean = javaClass.desiredAssertionStatus()

    inline fun assert(assertion: () -> Boolean) {
        if (enabled) {
            assert(assertion()) { "" }
        }
    }

    inline fun assert(assertion: () -> Boolean, lazyMessage: () -> Any) {
        if (enabled) {
            assert(assertion(), lazyMessage)
        }
    }


    fun assertLatLngIsValid(point: S2LatLng, lazyMessage: (() -> Any)? = null) {
        this.assert({ point.isValid }, {
            if (lazyMessage != null) lazyMessage() else "The LatLng point $point is not valid."
        })
    }

    fun assertCapIsValid(cap: S2Cap, lazyMessage: (() -> Any)? = null) {
        this.assert({ cap.isValid }, {
            if (lazyMessage != null) lazyMessage() else "The cap $cap is not valid."
        })
    }

    fun assertLatLngRectIsValid(rect: S2LatLngRect, lazyMessage: (() -> Any)? = null) {
        this.assert({ rect.isValid }, {
            if (lazyMessage != null) lazyMessage() else "The rectangle $rect is not valid."
        })
    }

    fun assertPointIsUnitLength(p: S2Point, lazyMessage: (() -> Any)? = null) {
        this.assert({ p.isUnitLength() }, {
            if (lazyMessage != null) lazyMessage() else "The point $p is not on the unit sphere. (|p|^2=${p.norm2()}, 5*eps = ${S2.DBL_EPSILON * 5}"
        })
    }

    fun assertGreaterOrEquals(d1: Double, d2: Double, lazyMessage: (() -> Any)? = null) {
        this.assert({ d1 >= d2 }, {
            if (lazyMessage != null) lazyMessage() else "$d1 is not greater or equals to $d2."
        })
    }
}