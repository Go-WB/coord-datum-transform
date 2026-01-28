/*
 * =====================================================================================
 *
 * Copyright (c) 2026 Zepp Health. All Rights Reserved. This computer program includes
 * Confidential, Proprietary Information and is a Trade Secret of Zepp Health Ltd.
 * All use, disclosure, and/or reproduction is prohibited unless authorized in writing.
 * Licensed under the MIT License. You can contact below email if need.
 *
 * version: 0.0.1
 * Author: wangwenbing@zepp.com
 *
 * =====================================================================================
 */

#ifndef COORD_TRANSFORM_H
#define COORD_TRANSFORM_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

// Forward declaration of GeographicLib struct
struct geod_geodesic;

// Coordinate format enum
typedef enum
{
    COORD_FORMAT_DD = 0,        // Decimal degrees DD.ddddd°
    COORD_FORMAT_DMM,           // Degrees and minutes DD°MM.mmm'
    COORD_FORMAT_DMS,           // Degrees, minutes, seconds DD°MM'SS"
    COORD_FORMAT_UTM,           // UTM coordinates
    COORD_FORMAT_MGRS,          // MGRS coordinates (default)
    COORD_FORMAT_BRITISH_GRID,  // British National Grid
    COORD_FORMAT_JAPAN_GRID,    // Japan grid
    COORD_FORMAT_MAX
} CoordFormat;

// Map datum enum
typedef enum
{
    DATUM_WGS84 = 0,            // World Geodetic System 1984
    DATUM_MGRS_GRID,           // MGRS Grid
    DATUM_UTM_GRID,            // UTM Grid
    DATUM_NAD83,               // North American Datum 1983
    DATUM_NAD27,               // North American Datum 1927
    DATUM_ED50,                // European Datum 1950
    DATUM_TOKYO,               // Tokyo Datum
    DATUM_OSGB36,              // Ordnance Survey of Great Britain 1936
    DATUM_MAX
} MapDatum;

// Ellipsoid parameters
typedef struct
{
    double a;           // Semi-major axis (meters)
    double f;           // Flattening
    double b;           // Semi-minor axis (meters)
    double e2;          // First eccentricity squared
    double ep2;         // Second eccentricity squared
    const char *name;   // Ellipsoid name
} Ellipsoid;

// 7-parameter datum transform
typedef struct
{
    double dx, dy, dz;          // Translation parameters (meters)
    double rx, ry, rz;          // Rotation parameters (arc-seconds)
    double scale;               // Scale factor (ppm)
} DatumTransform;

// Geographic coordinate
typedef struct
{
    double latitude;            // Latitude (degrees)
    double longitude;           // Longitude (degrees)
    double altitude;            // Altitude (meters)
    MapDatum datum;             // Coordinate datum
} GeoCoord;

// UTM coordinate
typedef struct
{
    int zone;                   // UTM zone (1-60)
    char band;                  // Latitude band (C-X)
    double easting;             // Easting (meters)
    double northing;            // Northing (meters)
    double convergence;         // Meridian convergence (degrees)
    double scale_factor;        // Scale factor
    MapDatum datum;             // Datum
} UTMPoint;

// MGRS coordinate
typedef struct
{
    int zone;                   // UTM zone (1-60)
    char band;                  // Latitude band (C-X)
    char square[3];             // 100km grid square (2 chars)
    double easting;             // Easting (meters, within grid square)
    double northing;            // Northing (meters, within grid square)
    MapDatum datum;             // Datum
} MGRSPoint;

// British National Grid coordinate
typedef struct
{
    char letters[3];           // Two-letter code
    double easting;             // Easting (meters)
    double northing;            // Northing (meters)
    MapDatum datum;             // Datum
} BritishGridPoint;

// Japan grid coordinate
typedef struct
{
    int zone;                   // Zone number
    double x;                   // X coordinate
    double y;                   // Y coordinate
    MapDatum datum;             // Datum
} JapanGridPoint;

// Parse result
typedef struct
{
    int success;                // Success flag
    GeoCoord coord;             // Parsed coordinate
    CoordFormat format;         // Detected format
    MapDatum datum;             // Detected datum
    char error_msg[256];        // Error message
} ParseResult;

// Geodesic result
typedef struct
{
    double distance;            // Distance (meters)
    double azimuth1;            // Forward azimuth (degrees)
    double azimuth2;            // Reverse azimuth (degrees)
} GeodesicResult;

// Coordinate transform context
typedef struct
{
    struct geod_geodesic *geod;  // Pointer to GeographicLib geodesic object
    Ellipsoid ellipsoid;        // Current ellipsoid
    DatumTransform transforms[DATUM_MAX][DATUM_MAX]; // Transform parameter table
} CoordContext;

// ============================ Public API ============================

// Error codes
#define COORD_SUCCESS 0
#define COORD_ERROR_INVALID_INPUT 1
#define COORD_ERROR_OUT_OF_RANGE 2
#define COORD_ERROR_PARSE_FAILED 3
#define COORD_ERROR_FORMAT 4
#define COORD_ERROR_MEMORY 5
#define COORD_ERROR_INVALID_COORD 6
#define COORD_ERROR_INVALID_UTM_ZONE 7
#define COORD_ERROR_DATUM_TRANSFORM 8
#define COORD_ERROR_CALCULATION 9
#define COORD_ERROR_UNSUPPORTED_FORMAT 10

// ==================== Initialization and cleanup ====================
CoordContext *coord_create_context(MapDatum datum);
void coord_destroy_context(CoordContext *ctx);
int coord_set_datum(CoordContext *ctx, MapDatum datum);

// ==================== Parsing functions ====================
ParseResult coord_parse_string(const char *str, CoordFormat format,
                               MapDatum datum);
ParseResult coord_auto_parse(const char *str);

// ==================== Formatting functions ====================
int coord_format_to_string(const GeoCoord *coord, CoordFormat format,
                           char *buffer, size_t buffer_size);
int coord_format_dd(const GeoCoord *coord, char *buffer, size_t buffer_size);
int coord_format_dmm(const GeoCoord *coord, char *buffer, size_t buffer_size);
int coord_format_dms(const GeoCoord *coord, char *buffer, size_t buffer_size);
int coord_format_utm(const UTMPoint *utm, char *buffer, size_t buffer_size);
int coord_format_mgrs(const MGRSPoint *mgrs, char *buffer, size_t buffer_size);
int coord_format_british_grid(const BritishGridPoint *bg, char *buffer,
                              size_t buffer_size);
int coord_format_japan_grid(const JapanGridPoint *jg, char *buffer,
                            size_t buffer_size);

// ==================== Coordinate conversion functions ====================
// Geographic coordinate to other formats
int coord_to_utm(CoordContext *ctx, const GeoCoord *geo, UTMPoint *utm);
int coord_from_utm(CoordContext *ctx, const UTMPoint *utm, GeoCoord *geo);
int coord_to_mgrs(CoordContext *ctx, const GeoCoord *geo, MGRSPoint *mgrs);
int coord_from_mgrs(CoordContext *ctx, const MGRSPoint *mgrs, GeoCoord *geo);
int coord_to_british_grid(CoordContext *ctx, const GeoCoord *geo,
                          BritishGridPoint *bg);
int coord_from_british_grid(CoordContext *ctx, const BritishGridPoint *bg,
                            GeoCoord *geo);
int coord_to_japan_grid(CoordContext *ctx, const GeoCoord *geo,
                        JapanGridPoint *jg);
int coord_from_japan_grid(CoordContext *ctx, const JapanGridPoint *jg,
                          GeoCoord *geo);

// Datum conversion
int coord_convert_datum(CoordContext *ctx, const GeoCoord *src,
                        MapDatum target_datum, GeoCoord *dst);

// ==================== Geodesic calculations ====================
int coord_distance(CoordContext *ctx, const GeoCoord *p1, const GeoCoord *p2,
                   double *distance, double *azi1, double *azi2);
int coord_direct(CoordContext *ctx, const GeoCoord *start,
                 double distance, double azimuth, GeoCoord *end);
int coord_inverse(CoordContext *ctx, const GeoCoord *p1, const GeoCoord *p2,
                  GeodesicResult *result);

// ==================== Utility functions ====================
int coord_get_utm_zone(double longitude, double latitude);
char coord_get_utm_band(double latitude);
int coord_validate_point(const GeoCoord *coord);
int coord_validate_utm(const UTMPoint *utm);
int coord_is_valid_latitude(double lat);
int coord_is_valid_longitude(double lon);
double coord_normalize_latitude(double lat);
double coord_normalize_longitude(double lon);
double coord_deg_to_rad(double deg);
double coord_rad_to_deg(double rad);
double coord_meters_to_feet(double meters);
double coord_feet_to_meters(double feet);

// ==================== Datum transform utilities ====================
int coord_set_transform_params(CoordContext *ctx, MapDatum from, MapDatum to,
                               const DatumTransform *params);
int coord_get_transform_params(CoordContext *ctx, MapDatum from, MapDatum to,
                               DatumTransform *params);

// ==================== Ellipsoid utilities ====================
const Ellipsoid *coord_get_ellipsoid(MapDatum datum);
int coord_set_custom_ellipsoid(CoordContext *ctx, double a, double f);

// ==================== Error handling ====================
const char *coord_get_error_string(int error_code);
void coord_set_error_callback(void (*callback)(int, const char *));

// ==================== Main format conversion function ====================
int coord_convert(CoordContext *ctx, const GeoCoord *src,
                  CoordFormat target_format, MapDatum target_datum,
                  char *result_buffer, size_t buffer_size);

#endif // COORD_TRANSFORM_H
