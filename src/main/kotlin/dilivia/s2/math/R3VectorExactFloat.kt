package dilivia.s2.math

import ch.obermuhlner.math.big.BigDecimalMath
import java.math.BigDecimal
import java.math.MathContext
import java.math.RoundingMode


object ExactFloatType : FloatingPointType<BigDecimal> {

    val mathContext = MathContext.UNLIMITED

    override val zero: BigDecimal = fromDouble(0.0)

    override val one: BigDecimal = fromDouble(1.0)

    override fun fromDouble(v: Double): BigDecimal {
        val bd = BigDecimal(v, mathContext)
        bd.setScale(2, RoundingMode.HALF_UP)
        return bd
    }

    override fun sqrt(v: BigDecimal): BigDecimal = v.sqrt(mathContext)

    override fun plus(a: BigDecimal, b: BigDecimal): BigDecimal = a.add(b, mathContext)

    override fun minus(a: BigDecimal, b: BigDecimal): BigDecimal = a.subtract(b, mathContext)

    override fun inv(a: BigDecimal): BigDecimal = one.multiply(a, mathContext)

    override fun times(a: BigDecimal, b: BigDecimal): BigDecimal = a.multiply(b, mathContext)

    override fun div(a: BigDecimal, b: BigDecimal): BigDecimal = a.divide(b, mathContext)

    override fun abs(a: BigDecimal): BigDecimal = a.abs(mathContext)

    override fun sum(values: List<BigDecimal>): BigDecimal = values.reduce { acc, bigDecimal -> plus(acc, bigDecimal) }

    override fun equals(v1: BigDecimal, v2: BigDecimal): Boolean = v1 == v2

    override fun atan2(y: BigDecimal, x: BigDecimal): BigDecimal = BigDecimalMath.atan2(y, x, mathContext)
}

@Strictfp
open class R3VectorExactFloat constructor(coords: List<BigDecimal>) : R3Vector<R3VectorExactFloat, BigDecimal>(coords, ExactFloatType) {

    init {
        require(coords.size == 3) { "Points must have exactly 3 coordinates" }
    }

    @JvmOverloads constructor(x: BigDecimal = ExactFloatType.zero, y: BigDecimal = ExactFloatType.zero, z: BigDecimal = ExactFloatType.zero): this(listOf(x, y ,z))

    @JvmOverloads constructor(x: Double, y: Double, z: Double): this(listOf(ExactFloatType.fromDouble(x), ExactFloatType.fromDouble(y) , ExactFloatType.fromDouble(z)))

    /**
     * Create a R3Vector instance from Int values.
     *
     * @param x x coordinate as a integer.
     * @param y y coordinate as a integer.
     * @param z z coordinate as a integer.
     */
    constructor(x: Int, y: Int, z: Int) : this(ExactFloatType.fromDouble(x.toDouble()), ExactFloatType.fromDouble(y.toDouble()), ExactFloatType.fromDouble(z.toDouble()))

    override fun newInstance(coords: List<BigDecimal>): R3VectorExactFloat = R3VectorExactFloat(coords)

    override fun toString(): String {
        return coords.map { "$it (p = ${it.precision()})" }.joinToString(prefix = "(", postfix = ")")
    }
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

fun Int.toXF(): BigDecimal = ExactFloatType.fromDouble(this.toDouble())
fun Double.toXF(): BigDecimal = ExactFloatType.fromDouble(this)