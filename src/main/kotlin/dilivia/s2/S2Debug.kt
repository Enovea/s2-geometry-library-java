package dilivia.s2

// Class that allows the --s2debug validity checks to be enabled or disabled
// for specific objects (e.g., see S2Polygon).
enum class S2Debug {
    ALLOW,    // Validity checks are controlled by --s2debug
    DISABLE   // No validity checks even when --s2debug is true
}
