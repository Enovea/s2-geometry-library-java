package dilivia.s2.math

import java.math.BigDecimal
import java.math.MathContext

@Strictfp
open class R3VectorLongDouble constructor(coords: List<BigDecimal>) : R3Vector<R3VectorLongDouble, BigDecimal>(coords, LongDoubleType()) {

    init {
        require(coords.size == 3) { "Points must have exactly 3 coordinates" }
    }

    @JvmOverloads constructor(x: BigDecimal = 0.0.toBigDecimal(MathContext.DECIMAL128), y: BigDecimal =0.0.toBigDecimal(MathContext.DECIMAL128), z: BigDecimal = 0.0.toBigDecimal(MathContext.DECIMAL128)): this(listOf(x, y ,z))

    @JvmOverloads constructor(x: Double, y: Double, z: Double): this(listOf(x.toBigDecimal(MathContext.DECIMAL128), y.toBigDecimal(MathContext.DECIMAL128) ,z.toBigDecimal(MathContext.DECIMAL128)))

    /**
     * Create a R3Vector instance from Int values.
     *
     * @param x x coordinate as a integer.
     * @param y y coordinate as a integer.
     * @param z z coordinate as a integer.
     */
    constructor(x: Int, y: Int, z: Int) : this(x.toBigDecimal(MathContext.DECIMAL128), y.toBigDecimal(MathContext.DECIMAL128), z.toBigDecimal(MathContext.DECIMAL128))

    override fun newInstance(coords: List<BigDecimal>): R3VectorLongDouble = R3VectorLongDouble(coords)

    companion object {
        @JvmStatic
        fun plus(p1: R3VectorLongDouble, p2: R3VectorLongDouble): R3VectorLongDouble {
            return p1 + p2
        }

        @JvmStatic
        fun times(p: R3VectorLongDouble, m: BigDecimal): R3VectorLongDouble {
            return p * m
        }

        @JvmStatic
        fun dotProd(p1: R3VectorLongDouble, p2: R3VectorLongDouble): BigDecimal {
            return p1.dotProd(p2)
        }
    }
}