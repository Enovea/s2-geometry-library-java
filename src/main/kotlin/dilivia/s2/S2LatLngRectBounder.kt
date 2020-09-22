/**
 * This project is a kotlin port of the Google s2 geometry library (Copyright 2005 Google Inc. All Rights Reserved.):
 *                                 https://github.com/google/s2geometry.git
 *
 * Copyright © 2020 Dilivia (contact@dilivia.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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
