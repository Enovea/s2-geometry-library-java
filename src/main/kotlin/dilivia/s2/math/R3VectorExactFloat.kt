package dilivia.s2.math

import java.math.BigDecimal
import java.math.MathContext

@Strictfp
open class R3VectorExactFloat constructor(coords: List<BigDecimal>) : R3Vector<R3VectorExactFloat, BigDecimal>(coords, ExactFloatType()) {

    init {
        require(coords.size == 3) { "Points must have exactly 3 coordinates" }
    }

    @JvmOverloads constructor(x: BigDecimal = 0.0.toBigDecimal(MathContext.UNLIMITED), y: BigDecimal =0.0.toBigDecimal(MathContext.UNLIMITED), z: BigDecimal = 0.0.toBigDecimal(MathContext.UNLIMITED)): this(listOf(x, y ,z))

    @JvmOverloads constructor(x: Double, y: Double, z: Double): this(listOf(x.toBigDecimal(MathContext.UNLIMITED), y.toBigDecimal(MathContext.UNLIMITED) ,z.toBigDecimal(MathContext.UNLIMITED)))

    /**
     * Create a R3Vector instance from Int values.
     *
     * @param x x coordinate as a integer.
     * @param y y coordinate as a integer.
     * @param z z coordinate as a integer.
     */
    constructor(x: Int, y: Int, z: Int) : this(x.toBigDecimal(MathContext.UNLIMITED), y.toBigDecimal(MathContext.UNLIMITED), z.toBigDecimal(MathContext.UNLIMITED))

    override fun newInstance(coords: List<BigDecimal>): R3VectorExactFloat = R3VectorExactFloat(coords)

    companion object {
        @JvmStatic
        fun plus(p1: R3VectorExactFloat, p2: R3VectorExactFloat): R3VectorExactFloat {
            return p1 + p2
        }

        @JvmStatic
        fun times(p: R3VectorExactFloat, m: BigDecimal): R3VectorExactFloat {
            return p * m
        }

        @JvmStatic
        fun dotProd(p1: R3VectorExactFloat, p2: R3VectorExactFloat): BigDecimal {
            return p1.dotProd(p2)
        }
    }
}