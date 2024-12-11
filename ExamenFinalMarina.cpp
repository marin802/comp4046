
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>

using namespace std;

// escoge un numero entero entre 1 a 360 para ver distintos angulos
const float frame = 35;

// constantes para el tamaño de imagen y fondo
const int Cw = 600, Ch = 660;
const double Vw = 1.0, Vh = 1.0;
const double d = 0.5;
const double BACKGROUND_COLOR[3] = {0, 0, 0}; // color de fondo azul claro

// clase para representar vectores 3D
class Vector3D {
public:
    double x, y, z;
    
    Vector3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    // calcula la longitud del vector
    double length() const {
        return sqrt(x * x + y * y + z * z);
    }

    // normaliza vector
    Vector3D normalize() const {
        double len = length();
        return Vector3D(x / len, y / len, z / len);
    }

    // sobrecarga del operador de resta para vectores
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }

    // sobrecarga del operador de suma para vectores
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }

    // sobrecarga del operador de multiplicación por un escalar
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    Vector3D operator-() const {
        return Vector3D(-x, -y, -z);
    }

    // producto punto entre dos vectores
    double dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // producto cruz entre dos vectores
    Vector3D cross(const Vector3D& other) const {
        return Vector3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
};

class Object {
public:
    virtual ~Object() = default;
};

// clase pra esferas
class Sphere : public Object {
public:
    Vector3D center;
    double radius;
    double color[3];
    double specular;
    double reflective;
    double refractive;

    // constructor de esferas
    Sphere(Vector3D c, double r, double col[3], double s, double refl, double refr)
        : center(c), radius(r), specular(s), reflective(refl), refractive(refr) {
        color[0] = col[0];
        color[1] = col[1];
        color[2] = col[2];
    }
};

// clase para luces en la escena
class Light{
public:
    string type;
    double intensity;
    Vector3D position, direction;

    // constructor de luces
    Light(string t, double i, Vector3D p = Vector3D(), Vector3D d = Vector3D())
        : type(t), intensity(i), position(p), direction(d) {}
};

class Cone : public Object {
public:
    Vector3D apex, base_center;
    double radius, height;
    double color[3];
    double specular, reflective;

    Cone(Vector3D apex, Vector3D base_center, double radius, double height, double col[3], double specular, double reflective)
        : apex(apex), base_center(base_center), radius(radius), height(height), specular(specular), reflective(reflective) {
        color[0] = col[0];
        color[1] = col[1];
        color[2] = col[2];
    }

    // Compute the axis direction from apex to base center
    Vector3D axis() const {
        return (base_center - apex).normalize();
    }
};


// configuracion de la escena con esferas, luces y triangulos
vector<Sphere> spheres;
vector<Light> lights;
vector<Cone> cones;

// Rotation matrix around the Y-axis (in the xy-plane)
const double angle = (frame/360) * M_PI *2;
const double cos_angle = cos(angle);
const double sin_angle = sin(angle);

auto rotate_coordinates = [](Vector3D& vec) {
    double x_new = cos_angle * vec.x - sin_angle * vec.z;
    double z_new = sin_angle * vec.x + cos_angle * vec.z;
    vec.x = x_new;
    vec.z = z_new;
};

void setup_scene() {
    double white[3] = {255, 255, 255};
    double pink[3] = {255, 192, 203};
    double gold[3] = {200, 150, 100};
    double black[3] = {0, 0, 0};

    // Cuerno
    cones.push_back(Cone(Vector3D(0, 3.6, -5.7), Vector3D(0, 2.7, -5.2), 0.25, 1.7, gold, 150, 0.5));

    // Cabeza
    spheres.push_back(Sphere(Vector3D(0, 2.2, -5), 0.75, white, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 1.85, -5.65), 0.27, white, 100, 0.15, 0));


    // Malena
    spheres.push_back(Sphere(Vector3D(0, 1.74, -4.12), 0.2, pink, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 2.11, -4.22), 0.2, pink, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 2.37, -4.32), 0.2, pink, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 2.66, -4.5), 0.2, pink, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 2.87, -4.75), 0.2, pink, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 2.95, -5.1), 0.2, pink, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 2.8, -5.35), 0.2, pink, 100, 0.15, 0));
    spheres.push_back(Sphere(Vector3D(0, 2.6, -5.65), 0.2, pink, 100, 0.15, 0));

    // Cuello
    cones.push_back(Cone(Vector3D(0, 2.9, -5), Vector3D(0, 0.6, -4), 0.75, 2.5, white, 100, 0.15));

    // Cuerpo
    spheres.push_back(Sphere(Vector3D(0, 0.3, -3.85), 1.2, white, 100, 0.15, 0));

    // Patas
    cones.push_back(Cone(Vector3D(0.6, -2, -4.5), Vector3D(0.6, 0, -4.5), 0.15, 2, white, 100, 0.15));
    cones.push_back(Cone(Vector3D(-0.6, -2, -4.5), Vector3D(-0.6, 0, -4.5), 0.15, 2, white, 100, 0.15));
    cones.push_back(Cone(Vector3D(0.6, -2, -3.25), Vector3D(0.6, 0, -3.45), 0.15, 2.3, white, 100, 0.15));
    cones.push_back(Cone(Vector3D(-0.6, -2, -3.25), Vector3D(-0.6, 0, -3.45), 0.15, 2.3, white, 100, 0.15));

    // Cola
    cones.push_back(Cone(Vector3D(0, 0.3, -1.88), Vector3D(0, 0.4, -5), 0.8, 2, pink, 50, 0.1));

    // Orejas
    cones.push_back(Cone(Vector3D(0.65, 2.95, -4.8), Vector3D(0.4, 2.2, -4.8), 0.2, 0.45, white, 100, 0.15));
    cones.push_back(Cone(Vector3D(-0.65, 2.95, -4.8), Vector3D(-0.4, 2.2, -4.8), 0.2, 0.45, white, 100, 0.15));

    // Ojos
    spheres.push_back(Sphere(Vector3D(0.56, 2.25, -5.3), 0.15, black, 200, 0.45, 0));
    spheres.push_back(Sphere(Vector3D(-0.56, 2.25, -5.3), 0.15, black, 200, 0.45, 0));

    // Suelo
    spheres.push_back(Sphere(Vector3D(0, -102, 0), 100.4, gold, 150, .5, 0.9));

    lights.push_back(Light("ambient", 0.7));  // luz ambiental
    lights.push_back(Light("point", 0.4, Vector3D(-1, 3, -5.5)));  // luz point
    lights.push_back(Light("directional", 0.2, Vector3D(-1,4, 1), Vector3D(-0.56, 2.25, -5.3)));  // luz direccional

    for (auto& cone : cones) {
        rotate_coordinates(cone.apex);
        rotate_coordinates(cone.base_center);
    }

    for (auto& sphere : spheres) {
        rotate_coordinates(sphere.center);
    }
}

Vector3D reflect_ray(const Vector3D& L, const Vector3D& N) {
    return N * 2 * N.dot(L) - L;
}

Vector3D refract_ray(const Vector3D& D, const Vector3D& N, double n1, double n2) {
    double eta = n1 / n2;
    double cos_theta_i = -D.dot(N);
    double sin2_theta_t = eta * eta * (1 - cos_theta_i * cos_theta_i);

    // total internal reflection check
    if (sin2_theta_t > 1.0) {
        return Vector3D(0, 0, 0);  // return zero vector for no refraction
    }

    double cos_theta_t = sqrt(1.0 - sin2_theta_t);
    return D * eta + N * (eta * cos_theta_i - cos_theta_t);
}

pair<double, double> intersect_ray_sphere(const Vector3D& O, const Vector3D& D, const Sphere& sphere) {
    Vector3D CO = O - sphere.center;
    double a = D.dot(D);
    double b = 2 * CO.dot(D);
    double c = CO.dot(CO) - sphere.radius * sphere.radius;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()}; // No intersection
    }

    double sqrt_discriminant = sqrt(discriminant);
    double t1 = (-b - sqrt_discriminant) / (2 * a);
    double t2 = (-b + sqrt_discriminant) / (2 * a);

    return {t1, t2};
}

pair<double, double> intersect_ray_cone(const Vector3D& O, const Vector3D& D, const Cone& cone) {
    Vector3D CO = O - cone.apex; // Vector from origin to cone apex
    Vector3D V = (cone.base_center - cone.apex).normalize(); // Cone's axis direction

    // Compute cone's angle squared
    double cos_theta = cone.height / sqrt(cone.height * cone.height + cone.radius * cone.radius);
    double cos_theta_sq = cos_theta * cos_theta;

    // Quadratic equation coefficients
    double a = D.dot(V) * D.dot(V) - cos_theta_sq * D.dot(D);
    double b = 2 * (D.dot(V) * CO.dot(V) - cos_theta_sq * D.dot(CO));
    double c = CO.dot(V) * CO.dot(V) - cos_theta_sq * CO.dot(CO);

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
    }

    double sqrt_discriminant = sqrt(discriminant);
    double t1 = (-b - sqrt_discriminant) / (2 * a);
    double t2 = (-b + sqrt_discriminant) / (2 * a);

    auto within_height = [&](double t) -> bool {
        Vector3D P = O + D * t; // Intersection point
        double height = (P - cone.apex).dot(V); // Projection onto cone axis
        return height >= 0 && height <= cone.height;
    };

    if (!within_height(t1)) t1 = numeric_limits<double>::infinity();
    if (!within_height(t2)) t2 = numeric_limits<double>::infinity();

    return {t1, t2};
}

Object* closest_intersection(const Vector3D& O, const Vector3D& D, double t_min, double t_max, double& closest_t) {
    closest_t = numeric_limits<double>::infinity();
    Object* closest_object = nullptr;

    for (auto& sphere : spheres) {
        auto [t1, t2] = intersect_ray_sphere(O, D, sphere);
        if (t_min < t1 && t1 < t_max && t1 < closest_t) {
            closest_t = t1;
            closest_object = &sphere;
        }
        if (t_min < t2 && t2 < t_max && t2 < closest_t) {
            closest_t = t2;
            closest_object = &sphere;
        }
    }

    for (auto& cone : cones) {
        auto [t1, t2] = intersect_ray_cone(O, D, cone);
        if (t_min < t1 && t1 < t_max && t1 < closest_t) {
            closest_t = t1;
            closest_object = &cone;
        }
        if (t_min < t2 && t2 < t_max && t2 < closest_t) {
            closest_t = t2;
            closest_object = &cone;
        }
    }

    return closest_object;
}


double compute_lighting(const Vector3D& P, const Vector3D& N, const Vector3D& V, double specular) {
    double intensity = 0.0;
    for (const auto& light : lights) {
        if (light.type == "ambient") {
            intensity += light.intensity;  // luz ambiental afecta a toda la escena por igual
        } else {
            Vector3D L;
            if (light.type == "point") {
                L = light.position - P;  // luz puntual desde una posicion específica
            } else {
                L = light.direction;  // luz direccional
            }

            L = L.normalize();
            double n_dot_l = N.dot(L);

            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (N.length() * L.length());
            }

            // cslculo de iluminacioon especular
            if (specular != -1) {
                Vector3D R = reflect_ray(L, N);
                double r_dot_v = R.dot(V);
                if (r_dot_v > 0) {
                    intensity += light.intensity * pow(r_dot_v / (R.length() * V.length()), specular);
                }
            }
        }
    }
    return intensity;
}

Vector3D canvas_to_viewport(int x, int y) {
    return Vector3D(x * Vw / Cw, y * Vh / Ch, d);
}

Vector3D trace_ray(const Vector3D& O, const Vector3D& D, int depth, Vector3D color) {
    double closest_t;
    Object* closest_object = closest_intersection(O, D, 1.0, numeric_limits<double>::infinity(), closest_t);

    if (!closest_object) {
        // Return background color
        color.x = BACKGROUND_COLOR[0];
        color.y = BACKGROUND_COLOR[1];
        color.z = BACKGROUND_COLOR[2];
        return color;
    }

    Vector3D P = O + D * closest_t;  // Intersection point
    Vector3D N;  // Normal at the intersection
    Vector3D V = D * -1;
    double lighting = 0.0;

    if (Sphere* sphere = dynamic_cast<Sphere*>(closest_object)) {
        // Sphere intersection
        N = (P - sphere->center).normalize();
        lighting = compute_lighting(P, N, V, sphere->specular);
        color = Vector3D(sphere->color[0], sphere->color[1], sphere->color[2]);

        if (depth > 0 && sphere->reflective > 0) {
            Vector3D R = reflect_ray(V, N);
            Vector3D P_offset = P + N * 0.001;  // Offset intersection point
            Vector3D reflected_color = trace_ray(P_offset, R, depth - 1, reflected_color);
            double reflective = sphere->reflective;
            if (sphere->refractive > 0) {
                Vector3D T = refract_ray(V, N, 1.0, 1.0);
                Vector3D P_offset = P - N * 0.001;  // Offset to avoid self-intersections
                Vector3D refracted_color = trace_ray(P_offset, T, depth - 1, refracted_color);
                double refractivity = sphere->refractive;
                color = color * lighting * (1 - reflective - refractivity)
                        + reflected_color * reflective
                        + refracted_color * refractivity;
            }else{
                color = color * lighting * (1 - reflective) + reflected_color * reflective; 
            }
        } else {
            color = color * lighting;
        }
    } else if (Cone* cone = dynamic_cast<Cone*>(closest_object)) {
        // Cone intersection
        Vector3D axis = cone->axis();
        N = (P - cone->apex).normalize() - axis * (P - cone->apex).dot(axis);
        N = N.normalize();
        lighting = compute_lighting(P, N, V, cone->specular);
        color = Vector3D(cone->color[0], cone->color[1], cone->color[2]);

        if (depth > 0 && cone->reflective > 0) {
            Vector3D R = reflect_ray(V, N);
            Vector3D P_offset = P + N * 0.001;  // Offset intersection point
            Vector3D reflected_color = trace_ray(P_offset, R, depth - 1, reflected_color);
            double reflective = cone->reflective;
            color = color * lighting * (1 - reflective) + reflected_color * reflective;
        } else {
            color = color * lighting;
        }
    }

    // Clamp color values
    color.x = min(255.0, max(0.0, color.x));
    color.y = min(255.0, max(0.0, color.y));
    color.z = min(255.0, max(0.0, color.z));

    return color;
}

void render_scene(ofstream& outFile) {
    outFile << "P3\n" << Cw << " " << Ch << "\n255\n";
    for (int y = Ch / 2; y > -Ch / 2; --y) {
        for (int x = -Cw / 2; x < Cw / 2; ++x) {
            Vector3D D = canvas_to_viewport(x, y);
            Vector3D color = trace_ray(Vector3D(0, 1.2, -9), D, 5, Vector3D(0,0,0));

            outFile << color.x << " " << color.y << " " << color.z << "\n";
        }
    }
}

int main() {
    ofstream outFile("Unicornio.ppm");
    setup_scene();
    render_scene(outFile);
    outFile.close();
    cout << "Imagen renderizada y guardada como 'Unicornio.ppm'." << endl;
    return 0;
}