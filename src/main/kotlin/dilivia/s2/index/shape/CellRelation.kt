package dilivia.s2.index.shape

/**
 * The possible relationships between a "target" cell and the cells of the S2ShapeIndex.
 *
 * - If the target is an index cell or is contained by an index cell, it is "INDEXED".
 * - If the target is subdivided into one or more index cells, it is "SUBDIVIDED".
 * - Otherwise it is "DISJOINT".
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
enum class CellRelation {
    INDEXED,       // Target is contained by an index cell
    SUBDIVIDED,    // Target is subdivided into one or more index cells
    DISJOINT       // Target does not intersect any index cells
}
