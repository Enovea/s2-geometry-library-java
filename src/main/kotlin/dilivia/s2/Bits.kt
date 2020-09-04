package dilivia.s2

import kotlin.random.Random
import kotlin.random.nextUInt

@ExperimentalUnsignedTypes
object Bits {

    fun findLSBSetNonZero(n: UInt): Int {
        var v = n
        var rc = 31
        var shift = 1 shl 4
        var i = 4
        while(i >= 0) {
            val x = v shl shift
            if (x != 0.toUInt()) {
                v = x;
                rc -= shift;
            }
            shift = shift shr 1
            --i
        }
        return rc
    }

}


@ExperimentalUnsignedTypes
fun main() {

    val n = 144U
    println(n.toString(2))
    println(Bits.findLSBSetNonZero(n))
    println("Leading 0 bits: ${n.countLeadingZeroBits()}")
    println("Trailing 0 bits: ${n.countTrailingZeroBits()}")

    val random = Random.Default
    repeat(30) {
        val n = random.nextUInt(4565445676L.toUInt())
        println("$n: ${n.toString(2)} ctzb = ${n.countTrailingZeroBits()} => ${Bits.findLSBSetNonZero(n)}")
    }
}