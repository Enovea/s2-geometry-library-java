package dilivia.s2.math

import ch.obermuhlner.math.big.BigDecimalMath
import java.math.BigDecimal
import java.math.MathContext

    class LongDoubleType : FloatingPointType<BigDecimal> {

    val mathContext = MathContext.DECIMAL128

    override val zero: BigDecimal = BigDecimal(0.0, mathContext)

    override val one: BigDecimal = BigDecimal(1.0, mathContext)

    override fun fromDouble(v: Double): BigDecimal = BigDecimal(v, mathContext)

    override fun sqrt(v: BigDecimal): BigDecimal = v.sqrt(mathContext)

    override fun plus(a: BigDecimal, b: BigDecimal): BigDecimal = a + b

    override fun minus(a: BigDecimal, b: BigDecimal): BigDecimal = a - b

    override fun inv(a: BigDecimal): BigDecimal = one / a

    override fun times(a: BigDecimal, b: BigDecimal): BigDecimal = a * b

    override fun div(a: BigDecimal, b: BigDecimal): BigDecimal = a / b

    override fun abs(a: BigDecimal): BigDecimal = a.abs(mathContext)

    override fun sum(values: List<BigDecimal>): BigDecimal = values.sumOf { it }

    override fun equals(v1: BigDecimal, v2: BigDecimal): Boolean = v1 == v2

    override fun atan2(y: BigDecimal, x: BigDecimal): BigDecimal = BigDecimalMath.atan2(y, x, mathContext)
}