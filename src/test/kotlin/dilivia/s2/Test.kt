package dilivia.s2

import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2LoopTest
import junit.framework.TestCase

class Test : TestCase() {

    fun test() {

        val loop_g: S2Loop = S2TextParser.makeLoop("0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")
    }
}