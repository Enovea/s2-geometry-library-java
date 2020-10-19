package dilivia.s2.index.shape

/**
 * Enum that represents the different status of a shape index.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
enum class S2ShapeIndexStatus {
    /** There are pending updates. */
    STALE,
    /** Updates are currently being applied. */
    UPDATING,
    /** There are no pending updates. */
    FRESH,
}
