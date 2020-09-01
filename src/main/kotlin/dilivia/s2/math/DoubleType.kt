package dilivia.s2.math

@Strictfp
class DoubleType() : FloatingPointType<Double> {

    override val zero: Double = 0.0

    override val one: Double = 1.0

    override fun fromDouble(v: Double): Double  = v

    override fun sqrt(v: Double): Double = kotlin.math.sqrt(v)

    override fun plus(a: Double, b: Double): Double = a + b

    override fun minus(a: Double, b: Double): Double = a - b

    override fun inv(a: Double): Double = 1.0 / a

    override fun times(a: Double, b: Double): Double = a * b

    override fun div(a: Double, b: Double): Double = a / b

    override fun abs(a: Double): Double = kotlin.math.abs(a)

    override fun sum(values: List<Double>): Double = values.sum()

    override fun equals(v1: Double, v2: Double): Boolean {
        return v1 == v2
    }

    override fun atan2(y: Double, x: Double): Double = kotlin.math.atan2(y, x)
}