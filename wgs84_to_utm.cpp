#include <iostream>
#include <cmath>

// Constants for WGS84 ellipsoid parameters
const double a = 6378137.0;  // semi-major axis in meters
const double f_inv = 298.257223563;  // inverse flattening

// Function to convert WGS84 coordinates to UTM
void convertWGS84toUTM(double lat, double lon, int& zone, double& easting, double& northing) {
    // Convert latitude and longitude to radians
    lat = lat * M_PI / 180.0;
    lon = lon * M_PI / 180.0;

    // Constants for UTM
    double k0 = 0.9996;  // scale factor
    double lon0 = -183.0;  // central meridian for zone 1

    // Calculate UTM zone
    zone = static_cast<int>((lon - lon0 * M_PI / 180.0) / 6.0) + 1;

    // Calculate central meridian for the UTM zone
    double lonc = lon0 + (zone - 1) * 6.0;

    // Calculate parameters for the Transverse Mercator projection
    double ecc = sqrt(1.0 - (1.0 - 1.0 / f_inv) * (1.0 - 1.0 / f_inv));
    double N = a / sqrt(1.0 - ecc * ecc * sin(lat) * sin(lat));
    double T = tan(lat) * tan(lat);
    double C = ecc * ecc * cos(lat) * cos(lat);
    double A = cos(lat) * (lon - lonc);
    double M = a * ((1.0 - ecc * ecc / 4.0 - 3.0 * ecc * ecc * ecc * ecc / 64.0 - 5.0 * ecc * ecc * ecc * ecc * ecc * ecc / 256.0) * lat
                  - (3.0 * ecc * ecc / 8.0 + 3.0 * ecc * ecc * ecc * ecc / 32.0 + 45.0 * ecc * ecc * ecc * ecc * ecc * ecc / 1024.0) * sin(2.0 * lat)
                  + (15.0 * ecc * ecc * ecc * ecc / 256.0 + 45.0 * ecc * ecc * ecc * ecc * ecc * ecc / 1024.0) * sin(4.0 * lat)
                  - (35.0 * ecc * ecc * ecc * ecc * ecc * ecc / 3072.0) * sin(6.0 * lat));

    // Calculate UTM coordinates
    easting = k0 * N * (A + (1.0 - T + C) * A * A * A / 6.0
                       + (5.0 - 18.0 * T + T * T + 72.0 * C - 58.0 * ecc * ecc) * A * A * A * A * A / 120.0)
              + 500000.0;

    northing = k0 * (M + N * tan(lat) * (A * A / 2.0
                                        + (5.0 - T + 9.0 * C + 4.0 * C * C) * A * A * A * A / 24.0
                                        + (61.0 - 58.0 * T + T * T + 600.0 * C - 330.0 * ecc * ecc) * A * A * A * A * A * A / 720.0));

    // Adjust northing for southern hemisphere
    if (lat < 0.0) {
        northing += 10000000.0;
    }
}

int main() {
    // Example WGS84 coordinates (latitude, longitude) in decimal degrees
    double latitude = 42.3309794;
    double longitude = 119.3169481;

    // Variables to store UTM coordinates
    int zone;
    double easting, northing;

    // Convert WGS84 to UTM
    convertWGS84toUTM(latitude, longitude, zone, easting, northing);

    // Display results
    std::cout << "UTM Zone: " << zone << std::endl;
    std::cout << "Easting: " << easting << " meters" << std::endl;
    std::cout << "Northing: " << northing << " meters" << std::endl;

    return 0;
}
