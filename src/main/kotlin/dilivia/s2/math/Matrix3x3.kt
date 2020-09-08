import dilivia.s2.Assertions.assertGE
import dilivia.s2.Assertions.assertLT
import dilivia.s2.Assertions.assertNE
import dilivia.s2.math.R2Vector
import dilivia.s2.math.R3Vector
import dilivia.s2.math.R3VectorDouble
import kotlin.math.cos
import kotlin.math.sin
import kotlin.math.sqrt

//
// A simple class to handle 3x3 matrices
// The aim of this class is to be able to manipulate 3x3 matrices
// and 3D vectors as naturally as possible and make calculations
// readable.
// For that reason, the operators +, -, * are overloaded.
// (Reading a = a + b*2 - c is much easier to read than
// a = Sub(Add(a, Mul(b,2)),c)   )
//
// Please be careful about overflows when using those matrices wth integer types
// The calculations are carried with VType. eg : if you are using uint8 as the
// base type, all values will be modulo 256.


class Matrix3x3(val m: Array<DoubleArray>) {

    // Initialize the matrix to 0
    constructor() : this(m = Array(3) { doubleArrayOf(0.0, 0.0, 0.0) })

    // Constructor explicitly setting the values of all the coefficient of
    // the matrix
    constructor(m00: Double, m01: Double, m02: Double,
                m10: Double, m11: Double, m12: Double,
                m20: Double, m21: Double, m22: Double) : this(
            arrayOf(
                    doubleArrayOf(m00, m01, m02),
                    doubleArrayOf(m10, m11, m12),
                    doubleArrayOf(m20, m21, m22)
            )
    )

    constructor(other: Matrix3x3) : this((0..2).map { other.m[it].copyOf() }.toTypedArray())

    // Change the value of all the coefficients of the matrix
    fun set(m00: Double, m01: Double, m02: Double,
            m10: Double, m11: Double, m12: Double,
            m20: Double, m21: Double, m22: Double): Matrix3x3 {
        m[0][0] = m00;
        m[0][1] = m01;
        m[0][2] = m02;

        m[1][0] = m10;
        m[1][1] = m11;
        m[1][2] = m12;

        m[2][0] = m20;
        m[2][1] = m21;
        m[2][2] = m22;
        return this
    }

    // Matrix addition
    operator fun plusAssign(mb: Matrix3x3) {
        m[0][0] += mb.m[0][0];
        m[0][1] += mb.m[0][1];
        m[0][2] += mb.m[0][2];

        m[1][0] += mb.m[1][0];
        m[1][1] += mb.m[1][1];
        m[1][2] += mb.m[1][2];

        m[2][0] += mb.m[2][0];
        m[2][1] += mb.m[2][1];
        m[2][2] += mb.m[2][2];
    }

    // Matrix subtration
    operator fun minusAssign(mb: Matrix3x3) {
        m[0][0] -= mb.m[0][0];
        m[0][1] -= mb.m[0][1];
        m[0][2] -= mb.m[0][2];

        m[1][0] -= mb.m[1][0];
        m[1][1] -= mb.m[1][1];
        m[1][2] -= mb.m[1][2];

        m[2][0] -= mb.m[2][0];
        m[2][1] -= mb.m[2][1];
        m[2][2] -= mb.m[2][2];
    }

    // Matrix multiplication by a scalar
    operator fun timesAssign(k: Double) {
        m[0][0] *= k;
        m[0][1] *= k;
        m[0][2] *= k;

        m[1][0] *= k;
        m[1][1] *= k;
        m[1][2] *= k;

        m[2][0] *= k;
        m[2][1] *= k;
        m[2][2] *= k;
    }

    // Matrix addition
    operator fun plus(mb: Matrix3x3): Matrix3x3 {
        val result = Matrix3x3(this)
        result += mb
        return result
    }

    // Matrix subtraction
    operator fun minus(mb: Matrix3x3): Matrix3x3 {
        val result = Matrix3x3(this)
        result -= mb
        return result
    }

    // Change the sign of all the coefficients in the matrix
    operator fun unaryMinus(): Matrix3x3 = Matrix3x3(
            -m[0][0], -m[0][1], -m[0][2],
            -m[1][0], -m[1][1], -m[1][2],
            -m[2][0], -m[2][1], -m[2][2]
    )

    // Matrix multiplication by a scalar
    operator fun times(k: Double): Matrix3x3 {
        val result = Matrix3x3(this)
        result *= k
        return result
    }

    // Matrix multiplication
    operator fun times(mb: Matrix3x3): Matrix3x3 = Matrix3x3(
            m[0][0] * mb.m[0][0] + m[0][1] * mb.m[1][0] + m[0][2] * mb.m[2][0],
            m[0][0] * mb.m[0][1] + m[0][1] * mb.m[1][1] + m[0][2] * mb.m[2][1],
            m[0][0] * mb.m[0][2] + m[0][1] * mb.m[1][2] + m[0][2] * mb.m[2][2],

            m[1][0] * mb.m[0][0] + m[1][1] * mb.m[1][0] + m[1][2] * mb.m[2][0],
            m[1][0] * mb.m[0][1] + m[1][1] * mb.m[1][1] + m[1][2] * mb.m[2][1],
            m[1][0] * mb.m[0][2] + m[1][1] * mb.m[1][2] + m[1][2] * mb.m[2][2],

            m[2][0] * mb.m[0][0] + m[2][1] * mb.m[1][0] + m[2][2] * mb.m[2][0],
            m[2][0] * mb.m[0][1] + m[2][1] * mb.m[1][1] + m[2][2] * mb.m[2][1],
            m[2][0] * mb.m[0][2] + m[2][1] * mb.m[1][2] + m[2][2] * mb.m[2][2]
    )

    // Multiplication of a matrix by a vector
    operator fun <V> times(v: V): V where V : R3Vector<V, Double> {
        return v.newInstance(listOf(
                m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
                m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
                m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
        ))
    }

    // Return the determinant of the matrix
    fun det(): Double =
            m[0][0] * m[1][1] * m[2][2] +
                    m[0][1] * m[1][2] * m[2][0] +
                    m[0][2] * m[1][0] * m[2][1] -
                    m[2][0] * m[1][1] * m[0][2] -
                    m[2][1] * m[1][2] * m[0][0] -
                    m[2][2] * m[1][0] * m[0][1]

    // Return the trace of the matrix
    fun trace(): Double {
        return m[0][0] + m[1][1] + m[2][2]
    }

    // Return matrix element (i,j) with 0<=i<=2 0<=j<=2
    fun get(i: Int, j: Int): Double {
        assertGE(i, 0)
        assertLT(i, 3)
        assertGE(j, 0)
        assertLT(j, 3)
        return m[i][j]
    }

    // Return matrix element (i/3,i%3) with 0<=i<=8 (access concatenated rows).
    operator fun get(i: Int): Double {
        assertGE(i, 0)
        assertLT(i, 9)
        return get(i / 3, i % 3)
    }

    // Return the transposed matrix
    fun transpose(): Matrix3x3 = Matrix3x3(
            m[0][0], m[1][0], m[2][0],
            m[0][1], m[1][1], m[2][1],
            m[0][2], m[1][2], m[2][2]
    )

    // Return the transposed of the matrix of the cofactors
    // (Useful for inversion for example)
    fun comatrixTransposed(): Matrix3x3 = Matrix3x3(
            m[1][1] * m[2][2] - m[2][1] * m[1][2],
            m[2][1] * m[0][2] - m[0][1] * m[2][2],
            m[0][1] * m[1][2] - m[1][1] * m[0][2],

            m[1][2] * m[2][0] - m[2][2] * m[1][0],
            m[2][2] * m[0][0] - m[0][2] * m[2][0],
            m[0][2] * m[1][0] - m[1][2] * m[0][0],

            m[1][0] * m[2][1] - m[2][0] * m[1][1],
            m[2][0] * m[0][1] - m[0][0] * m[2][1],
            m[0][0] * m[1][1] - m[1][0] * m[0][1])

    // Matrix inversion
    fun inverse(): Matrix3x3 {
        val det = det()
        assertNE(det, 0.0) { "Can't inverse. Determinant = 0." }
        return comatrixTransposed() * (1.0 / det)
    }

    // Return the vector 3D at row i
    fun row(i: Int): R3VectorDouble {
        assertGE(i, 0);
        assertLT(i, 3);
        return R3VectorDouble(m[i][0], m[i][1], m[i][2])
    }

    // Return the vector 3D at col i
    fun col(i: Int): R3VectorDouble {
        assertGE(i, 0);
        assertLT(i, 3);
        return R3VectorDouble(m[0][i], m[1][i], m[2][i])
    }

    // Set the vector in row i to be v1
    fun <V> setRow(i: Int, v : V)  where V : R3Vector<V, Double> {
        assertGE(i, 0);
        assertLT(i, 3);
        m[i][0] = v[0];
        m[i][1] = v[1];
        m[i][2] = v[2];
    }

    // Set the vector in column i to be v1
    fun <V> setCol(i: Int, v: V) where V : R3Vector<V, Double> {
        assertGE(i, 0);
        assertLT(i, 3);
        m[0][i] = v[0];
        m[1][i] = v[1];
        m[2][i] = v[2];
    }

    // Return a matrix M close to the original but verifying MtM = I
    // (useful to compensate for errors in a rotation matrix)
    fun orthogonalize(): Matrix3x3 {
        val r1 = row(0).normalize()
        val r2 = (row(2).crossProd(r1)).normalize()
        val r3 = (r1.crossProd(r2)).normalize()
        return fromRows(r1, r2, r3)
    }

    // Returns v.Transpose() * (*this) * u
    fun <V> mulBothSides(v: V, u: V): Double  where V : R3Vector<V, Double>{
        return (this * u).dotProd(v)
    }

    // Use the 3x3 matrix as a projective transform for 2d points
    fun project(v: R2Vector): R2Vector {
        val temp = this * R3VectorDouble(v[0], v[1], 1.0)
        return R2Vector(temp[0] / temp[2], temp[1] / temp[2])
    }

    // Return the Frobenius norm of the matrix: sqrt(sum(aij^2))
    fun frobeniusNorm(): Double {
        var sum = 0.0
        for (i in 0..2) {
            for (j in 0..2) {
            sum += m[i][j] * m[i][j];
        }
        }
        return sqrt(sum)
    }

    // Return true is one of the elements of the matrix is NaN
    fun isNaN(): Boolean = (0..8).any { get(it).isNaN() }


    override fun equals(other: Any?): Boolean = when {
        this === other -> true
        other !is Matrix3x3 -> false
        else -> (0..8).all { i -> this[i] == other[i] }
    }

    override fun hashCode(): Int {
        return m.contentDeepHashCode()
    }

    override fun toString(): String {
        var str = ""
        for (i in 0..2) {
            str += if (i == 0) "[" else " "

            for (j in  0..2) {
                str += "${get(i, j)} "
            }

            if (i == 2) str += "]"
        }

        return str
    }

    companion object {


        // Create a matrix from 3 row vectors
        fun <V> fromRows(v1: V, v2: V, v3: V) : Matrix3x3 where V : R3Vector<V, Double> = Matrix3x3(
                v1[0], v1[1], v1[2],
                v2[0], v2[1], v2[2],
                v3[0], v3[1], v3[2]
        )

        // Create a matrix from 3 column vectors
        fun <V> fromCols(v1: V, v2: V, v3: V): Matrix3x3  where V : R3Vector<V, Double> = Matrix3x3(
                v1[0], v2[0], v3[0],
                v1[1], v2[1], v3[1],
                v1[2], v2[2], v3[2]
        )

        // Return the identity matrix
        fun identity(): Matrix3x3 = Matrix3x3(
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0
        )

        // Return a matrix full of zeros
        fun zero(): Matrix3x3 = Matrix3x3()

        // Return a diagonal matrix with the coefficients in v
        fun <V> diagonal(v: V): Matrix3x3 where V : R3Vector<V, Double> {
            return Matrix3x3(
                    v[0], 0.0, 0.0,
                    0.0, v[1], 0.0,
                    0.0, 0.0, v[2]);
        }

        // Return the matrix vvT
       fun <V> sym3(v: V): Matrix3x3 where V : R3Vector<V, Double> = Matrix3x3(
               v[0] * v[0], v[0] * v[1], v[0] * v[2],
               v[1] * v[0], v[1] * v[1], v[1] * v[2],
               v[2] * v[0], v[2] * v[1], v[2] * v[2]
        )

        // Return a matrix M such that:
        // for each u,  M * u = v.CrossProd(u)
        fun <V> antiSym3(v: V): Matrix3x3 where V : R3Vector<V, Double> = Matrix3x3(
                0.0, -v[2], v[1],
                v[2], 0.0, -v[0],
                -v[1], v[0], 0.0
        )

        // Returns matrix that rotates |rot| radians around axis rot.
        fun <V> rodrigues(rot: V): Matrix3x3  where V : R3Vector<V, Double> {
            val theta = rot.norm()
            val w = rot.normalize()
            val wv = antiSym3(w)
            val ident = identity()
            val a = sym3(w)
            return (1 - cos(theta)) * a + sin(theta) * wv + cos(theta) * ident;

        }

    }
}

operator fun Double.times(m: Matrix3x3): Matrix3x3 = m * this

