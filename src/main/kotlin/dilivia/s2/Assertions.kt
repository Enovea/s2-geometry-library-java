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
package dilivia.s2

import com.google.common.geometry.S2
import dilivia.s2.region.S2Cap

/**
 * Assertion Helpers
 *
 * @author Fabien Meurisse
 * @since 1.0
 */
object Assertions {

    /** Indicates if the assertions are enabled. */
    val enabled: Boolean = javaClass.desiredAssertionStatus()

    /**
     *
     * @param assertion
     */
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

    fun assertNE(v1: Any, v2: Any, lazyMessage: (() -> Any)? = null) {
        this.assert( { v1 != v2 }, {
            if (lazyMessage != null) lazyMessage() else "The value $v1 is equals to $v2."
        })
    }

    fun assertEQ(v1: Any, v2: Any, lazyMessage: (() -> Any)? = null) {
        this.assert( { v1 == v2 }, {
            if (lazyMessage != null) lazyMessage() else "The value $v1 is not equals to $v2."
        })
    }

    fun <T: Comparable<T>> assertGE(v1: T, v2: T, lazyMessage: (() -> Any)? = null) {
        this.assert( { v1 >= v2 }, {
            if (lazyMessage != null) lazyMessage() else "The value $v1 is not greater than or equals value $v2."
        })
    }

    fun <T: Comparable<T>> assertGT(v1: T, v2: T, lazyMessage: (() -> Any)? = null) {
        this.assert( { v1 > v2 }, {
            if (lazyMessage != null) lazyMessage() else "The value $v1 is not greater than value $v2."
        })
    }
    fun <T: Comparable<T>> assertLE(v1: T, v2: T, lazyMessage: (() -> Any)? = null) {
        this.assert( { v1 <= v2 }, {
            if (lazyMessage != null) lazyMessage() else "The value $v1 is not lesser than or equals value $v2."
        })
    }

    fun <T: Comparable<T>> assertLT(v1: T, v2: T, lazyMessage: (() -> Any)? = null) {
        this.assert( { v1 < v2 }, {
            if (lazyMessage != null) lazyMessage() else "The value $v1 is not lesser than value $v2."
        })
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