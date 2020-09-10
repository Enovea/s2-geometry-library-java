package dilivia.s2

import com.google.common.geometry.S2.DBL_EPSILON

class S2LatLngRectBounder {

    companion object {

        fun maxErrorForTests(): S2LatLng {
            // The maximum error in the latitude calculation is
            //    3.84 * DBL_EPSILON   for the RobustCrossProd calculation
            //    0.96 * DBL_EPSILON   for the Latitude() calculation
            //    5    * DBL_EPSILON   added by AddPoint/GetBound to compensate for error
            //    ------------------
            //    9.80 * DBL_EPSILON   maximum error in result
            //
            // The maximum error in the longitude calculation is DBL_EPSILON.  GetBound
            // does not do any expansion because this isn't necessary in order to
            // bound the *rounded* longitudes of contained points.
            return S2LatLng.fromRadians(10 * DBL_EPSILON, 1 * DBL_EPSILON);
        }
    }
}
