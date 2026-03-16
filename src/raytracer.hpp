#pragma once

using Float = double;
//using Float = float;

template <typename T>
class Vector2;
template <typename T>
class Vector3;

template <typename T>
class Point2;
template <typename T>
class Point3;

template <typename T>
class Normal3;

template <typename T>
class Bounds3;

using Vector2f = Vector2<Float>;
using Vector3f = Vector3<Float>;
using Vector3i = Vector3<int>;

using Point2f = Point2<Float>;
using Point3f = Point3<Float>;
using Point3i = Point3<int>;

using Normal3f = Normal3<Float>;
using Normal3i = Normal3<int>;

using Bounds3f = Bounds3<Float>;

class SampledSpectrum;
class SampledWavelengths;

void init();
void cleanup();
