#ifndef MATERIAL_H
#define MATERIAL_H

#include "hittable.h"
#include "util/profiler.hpp"
#include "util/vecmath.hpp"
#include "util/color.hpp"
#include "raytracer.hpp"

struct material {
    virtual ~material() = default;

    virtual bool scatter(const Ray &r_in, 
                         const hit_record &rec, 
                         color &attenuation, 
                         Ray &scattered) const {
        return false;
    }
};

class lambertian : public material {
public:
    lambertian(const color &albedo) : albedo(albedo) {}

    bool scatter(const Ray &r_in, const hit_record &rec, color &attenuation, Ray &scattered)
    const override {
        PROFILE_SCOPE("lambertian::scatter");
        Vector3f scatter_direction = Vector3f(rec.normal) + random_unit_vector<Float>();

        if (scatter_direction.near_zero())
            scatter_direction = Vector3f(rec.normal);

        scattered = Ray(rec.p, scatter_direction);
        attenuation = albedo;
        return true;
    }

private:
    color albedo;
};

class metal : public material {
public:
    metal(const color &albedo, Float fuzz) : albedo(albedo), fuzz(fuzz < 1 ? fuzz : 1) {}

    bool scatter(const Ray &r_in, const hit_record &rec, color &attenuation, Ray &scattered)
    const override {
        PROFILE_SCOPE("metal::scatter");
        Vector3f reflected = reflect(normalize(r_in.d), rec.normal);
        reflected = reflected + (fuzz * random_unit_vector<Float>());
        scattered = Ray(rec.p, reflected);
        attenuation = albedo;
        return (dot(scattered.d, rec.normal) > 0);
    }

private:
    color albedo;
    Float fuzz;

};

class dielectric : public material {
public:
    dielectric(Float refraction_index) : refraction_index(refraction_index) {}

    bool scatter(const Ray &r_in, const hit_record &rec, color &attenuation, Ray &scattered)
    const override {
        PROFILE_SCOPE("dielectric::scatter");
        attenuation = color(1.0, 1.0, 1.0);
        Float ri = rec.front_face ? (1.0/refraction_index) : refraction_index;

        Vector3f unit_direction = normalize(r_in.d);
        Float cos_theta = std::fmin(dot(-unit_direction, rec.normal), 1.0);
        Float sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);

        bool cannnot_refract = ri * sin_theta > 1.0;
        Vector3f direction;
        if (cannnot_refract || reflectance(cos_theta, ri) > Rand::random<Float>())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, ri);

        scattered = Ray(rec.p, direction);
        return true;
    }

private:
    Float refraction_index;

    static Float reflectance(double cosine, double refraction_index) {
        auto r0 = (1 - refraction_index) / (1 + refraction_index);
        r0 = r0*r0;
        return r0 + (1 - r0) * std::pow((1 - cosine), 5);
    }
};

#endif
