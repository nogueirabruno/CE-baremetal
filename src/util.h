//
// Created by optlab on 17/01/25.
//

#ifndef UTIL_H
#define UTIL_H

#include <cmath>

// Compare if two doubles are approximately equal
inline bool areAlmostEqual(double a, double b, double epsilon = 1e-9) {
    return std::fabs(a - b) <= epsilon;
}

// Compare if 'a' is less than or equal to 'b' with a tolerance
inline bool isLessThanOrAlmostEqual(double a, double b, double epsilon = 1e-9) {
    return (a < b) || areAlmostEqual(a, b, epsilon);
}

// Compare if 'a' is less than 'b' with a tolerance
inline bool isLessThan(double a, double b, double epsilon = 1e-9) {
    return (b - a) > epsilon;
}

#endif //UTIL_H
