package dilivia.s2.index.shape

/**
 * When passed to an Iterator constructor, specifies whether the iterator should be positioned :
 * - at the beginning of the index (BEGIN),
 * - the end of the index (END),
 * - or arbitrarily (UNPOSITIONED).
 *
 * By default iterators are unpositioned, since this avoids an extra seek in this situation where one of the seek
 * methods (such as Locate) is immediately called.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
enum class InitialPosition { BEGIN, END, UNPOSITIONED }
