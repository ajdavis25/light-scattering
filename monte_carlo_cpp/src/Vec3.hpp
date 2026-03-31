#ifndef VEC3_HPP
#define VEC3_HPP

#include <cmath>

struct Vec3 {
    double x, y, z;
};

inline Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs)
{
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

inline Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs)
{
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

inline Vec3 operator*(const Vec3 &vec, double scalar)
{
    return {vec.x * scalar, vec.y * scalar, vec.z * scalar};
}

inline Vec3 operator*(double scalar, const Vec3 &vec)
{
    return vec * scalar;
}

inline Vec3 operator/(const Vec3 &vec, double scalar)
{
    return {vec.x / scalar, vec.y / scalar, vec.z / scalar};
}

inline double dot(const Vec3 &lhs, const Vec3 &rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline Vec3 cross(const Vec3 &lhs, const Vec3 &rhs)
{
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x,
    };
}

inline double norm(const Vec3 &vec)
{
    return std::sqrt(dot(vec, vec));
}

inline Vec3 normalize(const Vec3 &vec)
{
    const double magnitude = norm(vec);
    if (magnitude <= 0.0) {
        return {0.0, 0.0, 0.0};
    }
    return vec / magnitude;
}

#endif
