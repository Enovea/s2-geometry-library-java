package dilivia.s2

import com.google.common.geometry.S2
import kotlin.math.IEEErem

@Strictfp
class MutableS1Angle(override var radians: Double) : S1Angle(radians) {

    /**
     * Normalize this angle to the range (-180, 180] degrees.
     */
    fun normalize() {
        radians = radians.IEEErem(2.0 * S2.M_PI)
        if (radians <= -S2.M_PI) radians = S2.M_PI
    }

}