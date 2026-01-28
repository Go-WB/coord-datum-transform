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

#include "coord_datum_transform.h"
#include "geodesic.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

// Constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 / M_PI)
#define ARC_SEC_TO_RAD (M_PI / (180.0 * 3600.0))
#define PPM_TO_SCALE 1e-6
#define METERS_TO_FEET 3.280839895
#define FEET_TO_METERS 0.3048

// Ellipsoid definitions
static const Ellipsoid ELLIPSOIDS[] =
{
    // WGS84
    {
        6378137.0, 1.0 / 298.257223563, 6356752.314245,
        0.0066943799901413165, 0.0067394967422764341, "WGS84"
    },
    // MGRS Grid (uses WGS84 ellipsoid)
    {
        6378137.0, 1.0 / 298.257223563, 6356752.314245,
        0.0066943799901413165, 0.0067394967422764341, "WGS84"
    },
    // UTM Grid (uses WGS84 ellipsoid)
    {
        6378137.0, 1.0 / 298.257223563, 6356752.314245,
        0.0066943799901413165, 0.0067394967422764341, "WGS84"
    },
    // GRS80 (NAD83)
    {
        6378137.0, 1.0 / 298.257222101, 6356752.314140,
        0.006694380022903416, 0.006739496775478858, "GRS80"
    },
    // Clarke 1866 (NAD27)
    {
        6378206.4, 1.0 / 294.9786982, 6356583.8,
        0.006768658, 0.006814784, "Clarke1866"
    },
    // International 1924 (ED50)
    {
        6378388.0, 1.0 / 297.0, 6356911.946,
        0.006722670, 0.006768170, "Intl1924"
    },
    // Bessel 1841 (Tokyo)
    {
        6377397.155, 1.0 / 299.1528128, 6356078.963,
        0.006674372, 0.006719219, "Bessel1841"
    },
    // Airy 1830 (OSGB36 - British National Grid)
    {
        6377563.396, 1.0 / 299.3249646, 6356256.909,
        0.0066705397616, 0.006715826523, "Airy1830"
    }
};

// British National Grid parameters
static const double OSGB36_A = 6377563.396;    // Airy 1830 semi-major axis
static const double OSGB36_F = 1.0 / 299.3249646; // Airy 1830 flattening
static const double OSGB36_N0 = -100000.0;     // Northing offset
static const double OSGB36_E0 = 400000.0;     // Easting offset
static const double OSGB36_F0 = 0.9996012717; // Central meridian scale factor
static const double OSGB36_LAT0 = 49.0 * DEG_TO_RAD; // True origin latitude
static const double OSGB36_LON0 = -2.0 * DEG_TO_RAD; // True origin longitude

// Japan grid parameters (Tokyo Datum, Bessel 1841 ellipsoid)
static const double JAPAN_GRID_A = 6377397.155;
static const double JAPAN_GRID_F = 1.0 / 299.1528128;

// Error messages
static const char *ERROR_MESSAGES[] =
{
    "Success",
    "Invalid parameter",
    "Out of range",
    "Parse error",
    "Format error",
    "Memory allocation failed",
    "Invalid coordinate",
    "Invalid UTM zone",
    "Datum transformation failed",
    "Calculation error",
    "Unsupported format"
};

// Global error callback
static void (*error_callback)(int, const char *) = NULL;

// Set error
static void set_error(int code, const char *message)
{
    if (error_callback)
    {
        error_callback(code, message);
    }
}

// ==================== Basic utility functions ====================
int coord_is_valid_latitude(double lat)
{
    return lat >= -90.0 && lat <= 90.0;
}

int coord_is_valid_longitude(double lon)
{
    return lon >= -180.0 && lon <= 180.0;
}

double coord_normalize_latitude(double lat)
{
    if (lat > 90.0)
    {
        lat = 90.0;
    }
    if (lat < -90.0)
    {
        lat = -90.0;
    }
    return lat;
}

double coord_normalize_longitude(double lon)
{
    while (lon > 180.0)
    {
        lon -= 360.0;
    }
    while (lon < -180.0)
    {
        lon += 360.0;
    }
    return lon;
}

double coord_deg_to_rad(double deg)
{
    return deg * DEG_TO_RAD;
}

double coord_rad_to_deg(double rad)
{
    return rad * RAD_TO_DEG;
}

double coord_meters_to_feet(double meters)
{
    return meters * METERS_TO_FEET;
}

double coord_feet_to_meters(double feet)
{
    return feet * FEET_TO_METERS;
}

// ==================== Context management ====================
CoordContext *coord_create_context(MapDatum datum)
{
    if (datum >= DATUM_MAX)
    {
        return NULL;
    }
    CoordContext *ctx = (CoordContext *)malloc(sizeof(CoordContext));
    if (!ctx)
    {
        set_error(COORD_ERROR_MEMORY, "Memory allocation failed");
        return NULL;
    }
    memset(ctx, 0, sizeof(CoordContext));
    // Set ellipsoid
    ctx->ellipsoid = ELLIPSOIDS[datum];
    // Initialize GeographicLib geodesic object
    ctx->geod = (struct geod_geodesic *)malloc(sizeof(struct geod_geodesic));
    if (!ctx->geod)
    {
        free(ctx);
        set_error(COORD_ERROR_MEMORY, "Failed to create geodesic object");
        return NULL;
    }
    geod_init(ctx->geod, ctx->ellipsoid.a, ctx->ellipsoid.f);
    // Initialize transform parameter table
    memset(ctx->transforms, 0, sizeof(ctx->transforms));
    // Set default transform parameters
    // WGS84 <-> NAD83 (nearly identical)
    ctx->transforms[DATUM_WGS84][DATUM_NAD83].dx = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD83].dy = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD83].dz = 0.0;
    ctx->transforms[DATUM_NAD83][DATUM_WGS84] =
        ctx->transforms[DATUM_WGS84][DATUM_NAD83];
    // WGS84 <-> MGRS Grid (same)
    ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID].dx = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID].dy = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID].dz = 0.0;
    ctx->transforms[DATUM_MGRS_GRID][DATUM_WGS84] =
        ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID];
    // WGS84 <-> UTM Grid (same)
    ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID].dx = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID].dy = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID].dz = 0.0;
    ctx->transforms[DATUM_UTM_GRID][DATUM_WGS84] =
        ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID];
    // WGS84 -> NAD27 (NADCON parameters, CONUS)
    // Source: National Geodetic Survey
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].dx = -8.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].dy = 160.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].dz = 176.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].rx = -0.25;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].ry = 0.75;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].rz = -0.06;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].scale = -0.34;
    // WGS84 -> ED50 (EPSG parameters)
    // Source: EPSG Dataset
    ctx->transforms[DATUM_WGS84][DATUM_ED50].dx = -87.0;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].dy = -98.0;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].dz = -121.0;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].rx = -0.59;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].ry = -0.32;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].rz = -1.12;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].scale = -3.72;
    // WGS84 -> Tokyo (approximate parameters)
    ctx->transforms[DATUM_WGS84][DATUM_TOKYO].dx = -148.0;
    ctx->transforms[DATUM_WGS84][DATUM_TOKYO].dy = 507.0;
    ctx->transforms[DATUM_WGS84][DATUM_TOKYO].dz = 685.0;
    // WGS84 -> OSGB36 (OSTN15 parameters)
    // Source: Ordnance Survey National Grid (OSTN15)
    ctx->transforms[DATUM_WGS84][DATUM_OSGB36].dx = -446.448;
    ctx->transforms[DATUM_WGS84][DATUM_OSGB36].dy = 125.157;
    ctx->transforms[DATUM_WGS84][DATUM_OSGB36].dz = -542.060;
    ctx->transforms[DATUM_WGS84][DATUM_OSGB36].rx = -0.1502;
    ctx->transforms[DATUM_WGS84][DATUM_OSGB36].ry = -0.2470;
    ctx->transforms[DATUM_WGS84][DATUM_OSGB36].rz = -0.8421;
    ctx->transforms[DATUM_WGS84][DATUM_OSGB36].scale = 20.4894;
    return ctx;
}

void coord_destroy_context(CoordContext *ctx)
{
    if (ctx)
    {
        if (ctx->geod)
        {
            free(ctx->geod);
        }
        free(ctx);
    }
}

int coord_set_datum(CoordContext *ctx, MapDatum datum)
{
    if (!ctx || datum >= DATUM_MAX)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    ctx->ellipsoid = ELLIPSOIDS[datum];
    geod_init(ctx->geod, ctx->ellipsoid.a, ctx->ellipsoid.f);
    return COORD_SUCCESS;
}

// ==================== UTM zone calculation ====================
int coord_get_utm_zone(double longitude, double latitude)
{
    if (!coord_is_valid_longitude(longitude) || !coord_is_valid_latitude(latitude))
    {
        return 0;
    }
    // Normalize longitude
    double lon_norm = longitude;
    while (lon_norm < -180.0)
    {
        lon_norm += 360.0;
    }
    while (lon_norm >= 180.0)
    {
        lon_norm -= 360.0;
    }
    // Special-case regions
    if (latitude >= 56.0 && latitude < 64.0)
    {
        if (lon_norm >= 3.0 && lon_norm < 12.0)
        {
            return 32;
        }
    }
    if (latitude >= 72.0 && latitude < 84.0)
    {
        if (lon_norm >= 0.0 && lon_norm < 9.0)
        {
            return 31;
        }
        else if (lon_norm >= 9.0 && lon_norm < 21.0)
        {
            return 33;
        }
        else if (lon_norm >= 21.0 && lon_norm < 33.0)
        {
            return 35;
        }
        else if (lon_norm >= 33.0 && lon_norm < 42.0)
        {
            return 37;
        }
    }
    // Standard UTM zone
    int zone = (int)((lon_norm + 180.0) / 6.0) + 1;
    if (zone < 1)
    {
        zone = 1;
    }
    if (zone > 60)
    {
        zone = 60;
    }
    return zone;
}

char coord_get_utm_band(double latitude)
{
    if (latitude < -80.0)
    {
        return 'C';
    }
    if (latitude > 84.0)
    {
        return 'X';
    }
    // UTM latitude band table (8° per band, skipping I and O)
    struct
    {
        double min, max;
        char band;
    } bands[] =
    {
        {80.0, 84.0, 'X'},
        {72.0, 80.0, 'X'},
        {64.0, 72.0, 'W'},
        {56.0, 64.0, 'V'},
        {48.0, 56.0, 'U'},
        {40.0, 48.0, 'T'},
        {32.0, 40.0, 'S'},
        {24.0, 32.0, 'R'},
        {16.0, 24.0, 'Q'},
        {8.0, 16.0, 'P'},
        {0.0, 8.0, 'N'},
        {-8.0, 0.0, 'M'},
        {-16.0, -8.0, 'L'},
        {-24.0, -16.0, 'K'},
        {-32.0, -24.0, 'J'},
        {-40.0, -32.0, 'H'},
        {-48.0, -40.0, 'G'},
        {-56.0, -48.0, 'F'},
        {-64.0, -56.0, 'E'},
        {-72.0, -64.0, 'D'},
        {-80.0, -72.0, 'C'}
    };
    for (int i = 0; i < 21; i++)
    {
        if (latitude >= bands[i].min && latitude < bands[i].max)
        {
            return bands[i].band;
        }
    }
    return 'Z'; // Should not reach here
}

// ==================== Coordinate validation ====================
int coord_validate_point(const GeoCoord *coord)
{
    if (!coord)
    {
        return 0;
    }
    return coord_is_valid_latitude(coord->latitude) &&
           coord_is_valid_longitude(coord->longitude);
}

int coord_validate_utm(const UTMPoint *utm)
{
    if (!utm)
    {
        return 0;
    }
    if (utm->zone < 1 || utm->zone > 60)
    {
        return 0;
    }
    if (utm->band < 'C' || utm->band > 'X' ||
            utm->band == 'I' || utm->band == 'O')
    {
        return 0;
    }
    // Relax easting check range
    if (utm->easting < 100000.0 || utm->easting > 900000.0)
    {
        return 0;
    }
    // Relax northing check range, consider hemispheres
    if (utm->band >= 'N' && utm->band <= 'X')
    {
        // Northern hemisphere: 0-10,000,000 m
        if (utm->northing < 0.0 || utm->northing > 10000000.0)
        {
            return 0;
        }
    }
    else
    {
        // Southern hemisphere: 10,000,000-20,000,000 m (false northing)
        if (utm->northing < 10000000.0 || utm->northing > 20000000.0)
        {
            return 0;
        }
    }
    return 1;
}

static int coord_validate_mgrs(const MGRSPoint *mgrs)
{
    if (!mgrs)
    {
        return 0;
    }
    // Validate UTM zone
    if (mgrs->zone < 1 || mgrs->zone > 60)
    {
        return 0;
    }
    // Validate latitude band
    if (mgrs->band < 'C' || mgrs->band > 'X' ||
            mgrs->band == 'I' || mgrs->band == 'O')
    {
        return 0;
    }
    // Validate grid square identifier
    if (strlen(mgrs->square) != 2)
    {
        return 0;
    }
    // Validate grid square letters (skip I and O)
    if (mgrs->square[0] < 'A' || mgrs->square[0] > 'Z' ||
            mgrs->square[0] == 'I' || mgrs->square[0] == 'O' ||
            mgrs->square[1] < 'A' || mgrs->square[1] > 'Z' ||
            mgrs->square[1] == 'I' || mgrs->square[1] == 'O')
    {
        return 0;
    }
    // Validate easting and northing
    if (mgrs->easting < 0.0 || mgrs->easting > 99999.0)
    {
        return 0;
    }
    if (mgrs->northing < 0.0 || mgrs->northing > 99999.0)
    {
        return 0;
    }
    return 1;
}

// ==================== Coordinate parsing ====================
ParseResult coord_parse_string(const char *str, CoordFormat format,
                               MapDatum datum)
{
    ParseResult result = {0};
    result.success = 0;
    result.format = format;
    result.datum = datum;
    result.coord.altitude = 0.0;
    result.coord.datum = datum;
    if (!str)
    {
        strcpy(result.error_msg, "Input string is NULL");
        return result;
    }
    // Skip leading whitespace
    while (isspace((unsigned char)*str))
    {
        str++;
    }
    switch (format)
    {
        case COORD_FORMAT_DD:
        {
            // Format: "31.230416°N, 121.473701°E" or "31.230416, 121.473701"
            double lat, lon;
            char lat_dir = 'N', lon_dir = 'E';
            int count = sscanf(str, "%lf%*[ °]%c%*[ ,]%lf%*[ °]%c",
                               &lat, &lat_dir, &lon, &lon_dir);
            if (count == 4)
            {
                if (lat_dir == 'S' || lat_dir == 's')
                {
                    lat = -lat;
                }
                if (lon_dir == 'W' || lon_dir == 'w')
                {
                    lon = -lon;
                }
            }
            else
            {
                // Try format without direction letters
                count = sscanf(str, "%lf%*[ ,]%lf", &lat, &lon);
                if (count != 2)
                {
                    strcpy(result.error_msg, "Failed to parse DD format");
                    return result;
                }
            }
            if (!coord_is_valid_latitude(lat) || !coord_is_valid_longitude(lon))
            {
                strcpy(result.error_msg, "Coordinate out of range");
                return result;
            }
            result.coord.latitude = coord_normalize_latitude(lat);
            result.coord.longitude = coord_normalize_longitude(lon);
            result.success = 1;
            break;
        }
        case COORD_FORMAT_DMS:
        {
            // Format: "31°13'49.5\"N, 121°28'25.32\"E"
            int lat_deg, lat_min, lon_deg, lon_min;
            double lat_sec, lon_sec;
            char lat_dir, lon_dir;
            int count = sscanf(str, "%d°%d'%lf\"%c%*[ ,]%d°%d'%lf\"%c",
                               &lat_deg, &lat_min, &lat_sec, &lat_dir,
                               &lon_deg, &lon_min, &lon_sec, &lon_dir);
            if (count != 8)
            {
                strcpy(result.error_msg, "Failed to parse DMS format");
                return result;
            }
            double lat = lat_deg + lat_min / 60.0 + lat_sec / 3600.0;
            double lon = lon_deg + lon_min / 60.0 + lon_sec / 3600.0;
            if (lat_dir == 'S' || lat_dir == 's')
            {
                lat = -lat;
            }
            if (lon_dir == 'W' || lon_dir == 'w')
            {
                lon = -lon;
            }
            if (!coord_is_valid_latitude(lat) || !coord_is_valid_longitude(lon))
            {
                strcpy(result.error_msg, "Coordinate out of range");
                return result;
            }
            result.coord.latitude = coord_normalize_latitude(lat);
            result.coord.longitude = coord_normalize_longitude(lon);
            result.success = 1;
            break;
        }
        case COORD_FORMAT_DMM:
        {
            // Format: "31°13.825'N, 121°28.422'E"
            int lat_deg, lon_deg;
            double lat_min, lon_min;
            char lat_dir, lon_dir;
            int count = sscanf(str, "%d°%lf'%c%*[ ,]%d°%lf'%c",
                               &lat_deg, &lat_min, &lat_dir,
                               &lon_deg, &lon_min, &lon_dir);
            if (count != 6)
            {
                strcpy(result.error_msg, "Failed to parse DMM format");
                return result;
            }
            double lat = lat_deg + lat_min / 60.0;
            double lon = lon_deg + lon_min / 60.0;
            if (lat_dir == 'S' || lat_dir == 's')
            {
                lat = -lat;
            }
            if (lon_dir == 'W' || lon_dir == 'w')
            {
                lon = -lon;
            }
            if (!coord_is_valid_latitude(lat) || !coord_is_valid_longitude(lon))
            {
                strcpy(result.error_msg, "Coordinate out of range");
                return result;
            }
            result.coord.latitude = coord_normalize_latitude(lat);
            result.coord.longitude = coord_normalize_longitude(lon);
            result.success = 1;
            break;
        }
        case COORD_FORMAT_UTM:
        {
            // Format: "50N 447600E 4419300N" or "50N 447600 4419300"
            int zone;
            char band;
            double easting, northing;
            char east_dir = 'E', north_dir = 'N';
            int count = sscanf(str, "%d%c %lf%c %lf%c",
                               &zone, &band, &easting, &east_dir, &northing, &north_dir);
            if (count != 6)
            {
                // Try format without direction letters
                count = sscanf(str, "%d%c %lf %lf", &zone, &band, &easting, &northing);
                if (count != 4)
                {
                    strcpy(result.error_msg, "Failed to parse UTM format");
                    return result;
                }
            }
            // Create UTM point
            UTMPoint utm = {zone, band, easting, northing, 0.0, 0.9996, datum};
            if (!coord_validate_utm(&utm))
            {
                strcpy(result.error_msg, "Invalid UTM coordinate");
                return result;
            }
            // Convert to geographic coordinate
            CoordContext *ctx = coord_create_context(datum);
            if (!ctx)
            {
                strcpy(result.error_msg, "Failed to create context for UTM parsing");
                return result;
            }
            int ret = coord_from_utm(ctx, &utm, &result.coord);
            coord_destroy_context(ctx);
            if (ret != COORD_SUCCESS)
            {
                snprintf(result.error_msg, sizeof(result.error_msg),
                         "Failed to convert UTM to geographic: %s",
                         coord_get_error_string(ret));
                return result;
            }
            result.success = 1;
            break;
        }
        case COORD_FORMAT_MGRS:
        {
            // Format: "51Q SB 54634 56142" or "51QSB 54634 56142"
            int zone;
            char band;
            char square[3] = {0};  // Grid square (2 letters + null terminator)
            double easting, northing;
            // Try multiple formats
            int count = sscanf(str, "%d%c%2s %lf %lf",
                               &zone, &band, square, &easting, &northing);
            if (count != 5)
            {
                // Try spaced format: "51Q SB 54634 56142"
                count = sscanf(str, "%d%c %2s %lf %lf",
                               &zone, &band, square, &easting, &northing);
                if (count != 5)
                {
                    strcpy(result.error_msg, "Failed to parse MGRS format");
                    return result;
                }
            }
            // Validate MGRS parameters
            if (zone < 1 || zone > 60)
            {
                strcpy(result.error_msg, "Invalid MGRS zone (1-60)");
                return result;
            }
            if (band < 'C' || band > 'X' || band == 'I' || band == 'O')
            {
                strcpy(result.error_msg, "Invalid MGRS band");
                return result;
            }
            if (strlen(square) != 2)
            {
                strcpy(result.error_msg, "Invalid MGRS square (must be 2 letters)");
                return result;
            }
            // Validate grid square letters
            if (square[0] < 'A' || square[0] > 'Z' || square[0] == 'I' || square[0] == 'O'
                    ||
                    square[1] < 'A' || square[1] > 'Z' || square[1] == 'I' || square[1] == 'O')
            {
                strcpy(result.error_msg, "Invalid MGRS square letters");
                return result;
            }
            // Validate easting and northing
            if (easting < 0.0 || easting > 100000.0)
            {
                strcpy(result.error_msg, "MGRS easting must be 0-100000 meters");
                return result;
            }
            if (northing < 0.0 || northing > 100000.0)
            {
                strcpy(result.error_msg, "MGRS northing must be 0-100000 meters");
                return result;
            }
            // Create MGRS point
            MGRSPoint mgrs;
            mgrs.zone = zone;
            mgrs.band = band;
            mgrs.square[0] = square[0];
            mgrs.square[1] = square[1];
            mgrs.square[2] = '\0';
            mgrs.easting = easting;
            mgrs.northing = northing;
            mgrs.datum = datum;
            // Validate with coord_validate_mgrs
            if (!coord_validate_mgrs(&mgrs))
            {
                strcpy(result.error_msg, "Invalid MGRS coordinate");
                return result;
            }
            // Convert to geographic coordinate
            CoordContext *ctx = coord_create_context(datum);
            if (!ctx)
            {
                strcpy(result.error_msg, "Failed to create context for MGRS parsing");
                return result;
            }
            int ret = coord_from_mgrs(ctx, &mgrs, &result.coord);
            coord_destroy_context(ctx);
            if (ret != COORD_SUCCESS)
            {
                snprintf(result.error_msg, sizeof(result.error_msg),
                         "Failed to convert MGRS to geographic: %s",
                         coord_get_error_string(ret));
                return result;
            }
            result.success = 1;
            break;
        }
        case COORD_FORMAT_BRITISH_GRID:
        {
            // Format: "TQ 12345 67890" or "TQ1234567890"
            char letters[3] = {0};
            double easting, northing;
            int count = sscanf(str, "%2s %lf %lf", letters, &easting, &northing);
            if (count != 3)
            {
                // Try format without spaces: "TQ1234567890"
                char buffer[32];
                strncpy(buffer, str, sizeof(buffer) - 1);
                buffer[sizeof(buffer) - 1] = '\0';
                if (strlen(buffer) >= 2)
                {
                    letters[0] = buffer[0];
                    letters[1] = buffer[1];
                    letters[2] = '\0';
                    // Parse numeric part
                    if (sscanf(buffer + 2, "%lf%lf", &easting, &northing) == 2)
                    {
                        count = 3;
                    }
                }
            }
            if (count != 3)
            {
                strcpy(result.error_msg, "Failed to parse British Grid format");
                return result;
            }
            // Create British Grid point
            BritishGridPoint bg;
            strncpy(bg.letters, letters, 3);
            bg.easting = easting;
            bg.northing = northing;
            bg.datum = datum;
            // Convert to geographic coordinate
            CoordContext *ctx = coord_create_context(datum);
            if (!ctx)
            {
                strcpy(result.error_msg, "Failed to create context for British Grid parsing");
                return result;
            }
            int ret = coord_from_british_grid(ctx, &bg, &result.coord);
            coord_destroy_context(ctx);
            if (ret != COORD_SUCCESS)
            {
                snprintf(result.error_msg, sizeof(result.error_msg),
                         "Failed to convert British Grid to geographic: %s",
                         coord_get_error_string(ret));
                return result;
            }
            result.success = 1;
            break;
        }
        case COORD_FORMAT_JAPAN_GRID:
        {
            // Format: "Zone 3: 12345.6, 67890.1" or "3 12345.6 67890.1"
            int zone;
            double x, y;
            int count = sscanf(str, "Zone %d: %lf, %lf", &zone, &x, &y);
            if (count != 3)
            {
                count = sscanf(str, "%d %lf %lf", &zone, &x, &y);
            }
            if (count != 3)
            {
                strcpy(result.error_msg, "Failed to parse Japan Grid format");
                return result;
            }
            // Create Japan Grid point
            JapanGridPoint jg;
            jg.zone = zone;
            jg.x = x;
            jg.y = y;
            jg.datum = datum;
            // Convert to geographic coordinate
            CoordContext *ctx = coord_create_context(datum);
            if (!ctx)
            {
                strcpy(result.error_msg, "Failed to create context for Japan Grid parsing");
                return result;
            }
            int ret = coord_from_japan_grid(ctx, &jg, &result.coord);
            coord_destroy_context(ctx);
            if (ret != COORD_SUCCESS)
            {
                snprintf(result.error_msg, sizeof(result.error_msg),
                         "Failed to convert Japan Grid to geographic: %s",
                         coord_get_error_string(ret));
                return result;
            }
            result.success = 1;
            break;
        }
        default:
            snprintf(result.error_msg, sizeof(result.error_msg),
                     "Unsupported format: %d", format);
            break;
    }
    return result;
}

// In coord_auto_parse, add auto-detection for MGRS format
ParseResult coord_auto_parse(const char *str)
{
    ParseResult result = {0};
    if (!str)
    {
        strcpy(result.error_msg, "Input string is NULL");
        return result;
    }
    // Skip leading whitespace
    const char *s = str;
    while (isspace((unsigned char)*s))
    {
        s++;
    }
    // Check for MGRS format
    int zone;
    char band;
    char square[3];
    double easting, northing;
    // Try parsing MGRS format
    int count = sscanf(s, "%d%c%2s %lf %lf", &zone, &band, square, &easting,
                       &northing);
    if (count == 5)
    {
        // Validate MGRS parameters
        if (zone >= 1 && zone <= 60 &&
                band >= 'C' && band <= 'X' && band != 'I' && band != 'O' &&
                strlen(square) == 2 &&
                square[0] >= 'A' && square[0] <= 'Z' && square[0] != 'I' && square[0] != 'O' &&
                square[1] >= 'A' && square[1] <= 'Z' && square[1] != 'I' && square[1] != 'O' &&
                easting >= 0.0 && easting <= 100000.0 &&
                northing >= 0.0 && northing <= 100000.0)
        {
            // Looks like MGRS format
            result = coord_parse_string(str, COORD_FORMAT_MGRS, DATUM_WGS84);
            if (result.success)
            {
                return result;
            }
        }
    }
    // Check for UTM format
    char east_dir, north_dir;
    count = sscanf(s, "%d%c %lf%c %lf%c", &zone, &band, &easting, &east_dir,
                   &northing, &north_dir);
    if (count == 6
            || (count = sscanf(s, "%d%c %lf %lf", &zone, &band, &easting, &northing)) == 4)
    {
        if (zone >= 1 && zone <= 60 &&
                band >= 'C' && band <= 'X' && band != 'I' && band != 'O')
        {
            result = coord_parse_string(str, COORD_FORMAT_UTM, DATUM_WGS84);
            if (result.success)
            {
                return result;
            }
        }
    }
    // Check for British Grid format
    char letters[3];
    double east, north;
    if (sscanf(s, "%2s %lf %lf", letters, &east, &north) == 3)
    {
        if (isalpha((unsigned char)letters[0]) && isalpha((unsigned char)letters[1]))
        {
            result = coord_parse_string(str, COORD_FORMAT_BRITISH_GRID, DATUM_ED50);
            if (result.success)
            {
                return result;
            }
        }
    }
    // Check for Japan Grid format
    int j_zone;
    double x, y;
    if (sscanf(s, "Zone %d: %lf, %lf", &j_zone, &x, &y) == 3 ||
            sscanf(s, "%d %lf %lf", &j_zone, &x, &y) == 3)
    {
        result = coord_parse_string(str, COORD_FORMAT_JAPAN_GRID, DATUM_TOKYO);
        if (result.success)
        {
            return result;
        }
    }
    // Try other formats
    CoordFormat formats[] = {COORD_FORMAT_DD, COORD_FORMAT_DMS, COORD_FORMAT_DMM};
    MapDatum datum = DATUM_WGS84;
    for (int i = 0; i < 3; i++)
    {
        result = coord_parse_string(str, formats[i], datum);
        if (result.success)
        {
            return result;
        }
    }
    strcpy(result.error_msg, "Failed to auto-parse coordinate string");
    return result;
}


// ==================== Coordinate formatting functions ====================
int coord_format_to_string(const GeoCoord *coord, CoordFormat format,
                           char *buffer, size_t buffer_size)
{
    if (!coord || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(coord))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    switch (format)
    {
        case COORD_FORMAT_DD:
            return coord_format_dd(coord, buffer, buffer_size);
        case COORD_FORMAT_DMM:
            return coord_format_dmm(coord, buffer, buffer_size);
        case COORD_FORMAT_DMS:
            return coord_format_dms(coord, buffer, buffer_size);
        default:
            return COORD_ERROR_UNSUPPORTED_FORMAT;
    }
}

int coord_format_dd(const GeoCoord *coord, char *buffer, size_t buffer_size)
{
    if (!coord || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    char lat_dir = (coord->latitude >= 0.0) ? 'N' : 'S';
    char lon_dir = (coord->longitude >= 0.0) ? 'E' : 'W';
    double lat_abs = fabs(coord->latitude);
    double lon_abs = fabs(coord->longitude);
    int written = snprintf(buffer, buffer_size, "%.6f°%c, %.6f°%c",
                           lat_abs, lat_dir, lon_abs, lon_dir);
    return (written < 0
            || (size_t)written >= buffer_size) ? COORD_ERROR_FORMAT : COORD_SUCCESS;
}

int coord_format_dmm(const GeoCoord *coord, char *buffer, size_t buffer_size)
{
    if (!coord || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    char lat_dir = (coord->latitude >= 0.0) ? 'N' : 'S';
    char lon_dir = (coord->longitude >= 0.0) ? 'E' : 'W';
    double lat_abs = fabs(coord->latitude);
    double lon_abs = fabs(coord->longitude);
    int lat_deg = (int)lat_abs;
    double lat_min = (lat_abs - lat_deg) * 60.0;
    int lon_deg = (int)lon_abs;
    double lon_min = (lon_abs - lon_deg) * 60.0;
    int written = snprintf(buffer, buffer_size, "%d°%.3f'%c, %d°%.3f'%c",
                           lat_deg, lat_min, lat_dir,
                           lon_deg, lon_min, lon_dir);
    return (written < 0
            || (size_t)written >= buffer_size) ? COORD_ERROR_FORMAT : COORD_SUCCESS;
}

int coord_format_dms(const GeoCoord *coord, char *buffer, size_t buffer_size)
{
    if (!coord || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    char lat_dir = (coord->latitude >= 0.0) ? 'N' : 'S';
    char lon_dir = (coord->longitude >= 0.0) ? 'E' : 'W';
    double lat_abs = fabs(coord->latitude);
    double lon_abs = fabs(coord->longitude);
    int lat_deg = (int)lat_abs;
    double lat_remainder = (lat_abs - lat_deg) * 60.0;
    int lat_min = (int)lat_remainder;
    double lat_sec = (lat_remainder - lat_min) * 60.0;
    int lon_deg = (int)lon_abs;
    double lon_remainder = (lon_abs - lon_deg) * 60.0;
    int lon_min = (int)lon_remainder;
    double lon_sec = (lon_remainder - lon_min) * 60.0;
    int written = snprintf(buffer, buffer_size, "%d°%d'%.2f\"%c, %d°%d'%.2f\"%c",
                           lat_deg, lat_min, lat_sec, lat_dir,
                           lon_deg, lon_min, lon_sec, lon_dir);
    return (written < 0
            || (size_t)written >= buffer_size) ? COORD_ERROR_FORMAT : COORD_SUCCESS;
}

int coord_format_utm(const UTMPoint *utm, char *buffer, size_t buffer_size)
{
    if (!utm || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    int written = snprintf(buffer, buffer_size, "%d%c %.0fE %.0fN",
                           utm->zone, utm->band, utm->easting, utm->northing);
    return (written < 0
            || (size_t)written >= buffer_size) ? COORD_ERROR_FORMAT : COORD_SUCCESS;
}

int coord_format_mgrs(const MGRSPoint *mgrs, char *buffer, size_t buffer_size)
{
    if (!mgrs || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    int written = snprintf(buffer, buffer_size, "%d%c %s %05.0f %05.0f",
                           mgrs->zone, mgrs->band, mgrs->square,
                           mgrs->easting, mgrs->northing);
    return (written < 0
            || (size_t)written >= buffer_size) ? COORD_ERROR_FORMAT : COORD_SUCCESS;
}

int coord_format_british_grid(const BritishGridPoint *bg, char *buffer,
                              size_t buffer_size)
{
    if (!bg || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    int written = snprintf(buffer, buffer_size, "%s %.0f %.0f",
                           bg->letters, bg->easting, bg->northing);
    return (written < 0
            || (size_t)written >= buffer_size) ? COORD_ERROR_FORMAT : COORD_SUCCESS;
}

int coord_format_japan_grid(const JapanGridPoint *jg, char *buffer,
                            size_t buffer_size)
{
    if (!jg || !buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    int written = snprintf(buffer, buffer_size, "Zone %d: %.3f, %.3f",
                           jg->zone, jg->x, jg->y);
    return (written < 0
            || (size_t)written >= buffer_size) ? COORD_ERROR_FORMAT : COORD_SUCCESS;
}

// ==================== Coordinate conversion functions ====================
// Geographic coordinate to UTM
int coord_to_utm(CoordContext *ctx, const GeoCoord *geo, UTMPoint *utm)
{
    if (!ctx || !geo || !utm)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(geo))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    // Calculate UTM zone
    int zone = coord_get_utm_zone(geo->longitude, geo->latitude);
    if (zone < 1 || zone > 60)
    {
        return COORD_ERROR_INVALID_UTM_ZONE;
    }
    // Calculate central meridian
    double lon_center = (zone - 1) * 6.0 - 180.0 + 3.0;
    // Convert to radians
    double lat_rad = coord_deg_to_rad(geo->latitude);
    double lon_rad = coord_deg_to_rad(geo->longitude);
    double lon_center_rad = coord_deg_to_rad(lon_center);
    // UTM conversion parameters
    double k0 = 0.9996;  // UTM scale factor
    double a = ctx->ellipsoid.a;
    double f = ctx->ellipsoid.f;
    double e2 = 2 * f - f * f;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double tan_lat = sin_lat / cos_lat;
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double T = tan_lat * tan_lat;
    double C = e2 * cos_lat * cos_lat / (1.0 - e2);
    double A = (lon_rad - lon_center_rad) * cos_lat;
    // Compute M (meridional arc length)
    double M = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 /
                     256.0) * lat_rad
                    - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(
                        2.0 * lat_rad)
                    + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(4.0 * lat_rad)
                    - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * lat_rad));
    // Compute UTM coordinates
    double A2 = A * A;
    double A3 = A2 * A;
    double A4 = A3 * A;
    double A5 = A4 * A;
    double A6 = A5 * A;
    // Easting
    utm->easting = k0 * N * (A + (1.0 - T + C) * A3 / 6.0
                             + (5.0 - 18.0 * T + T * T + 72.0 * C - 58.0 * e2) * A5 / 120.0)
                   + 500000.0;  // False easting
    // Northing
    utm->northing = k0 * (M + N * tan_lat *
                          (A2 / 2.0 + (5.0 - T + 9.0 * C + 4.0 * C * C) * A4 / 24.0
                           + (61.0 - 58.0 * T + T * T + 600.0 * C - 330.0 * e2) * A6 / 720.0));
    // If southern hemisphere, add false northing
    if (geo->latitude < 0.0)
    {
        utm->northing += 10000000.0;
    }
    utm->zone = zone;
    utm->band = coord_get_utm_band(geo->latitude);
    utm->convergence = atan(tan_lat * sin(lon_rad - lon_center_rad));
    utm->scale_factor = k0;
    utm->datum = geo->datum;
    return COORD_SUCCESS;
}

int coord_from_utm(CoordContext *ctx, const UTMPoint *utm, GeoCoord *geo)
{
    if (!ctx || !utm || !geo)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_utm(utm))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    // Calculate central meridian
    double lon_center = (utm->zone - 1) * 6.0 - 180.0 + 3.0;
    double k0 = 0.9996;
    double a = ctx->ellipsoid.a;
    double f = ctx->ellipsoid.f;
    double e2 = 2 * f - f * f;
    // Remove false easting
    double x = utm->easting - 500000.0;
    double y = utm->northing;
    // If southern hemisphere, remove false northing
    if (utm->band < 'N')
    {
        y -= 10000000.0;
    }
    // Compute footpoint latitude
    double M = y / k0;
    double mu = M / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2
                          / 256.0));
    double e1 = (1.0 - sqrt(1.0 - e2)) / (1.0 + sqrt(1.0 - e2));
    double J1 = 3.0 * e1 / 2.0 - 27.0 * e1 * e1 * e1 / 32.0;
    double J2 = 21.0 * e1 * e1 / 16.0 - 55.0 * e1 * e1 * e1 * e1 / 32.0;
    double J3 = 151.0 * e1 * e1 * e1 / 96.0;
    double J4 = 1097.0 * e1 * e1 * e1 * e1 / 512.0;
    double fp = mu + J1 * sin(2.0 * mu) + J2 * sin(4.0 * mu) + J3 * sin(
                    6.0 * mu) + J4 * sin(8.0 * mu);
    double sin_fp = sin(fp);
    double cos_fp = cos(fp);
    double tan_fp = sin_fp / cos_fp;
    double C1 = e2 * cos_fp * cos_fp;
    double T1 = tan_fp * tan_fp;
    double R1 = a * (1.0 - e2) / pow(1.0 - e2 * sin_fp * sin_fp, 1.5);
    double N1 = a / sqrt(1.0 - e2 * sin_fp * sin_fp);
    double D = x / (N1 * k0);
    double Q1 = N1 * tan_fp / R1;
    double Q2 = D * D / 2.0;
    double Q3 = (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C1 * C1 - 9.0 * e2) * pow(D,
                4) / 24.0;
    double Q4 = (61.0 + 90.0 * T1 + 298.0 * C1 + 45.0 * T1 * T1 - 252.0 * e2 - 3.0 *
                 C1 * C1) * pow(D, 6) / 720.0;
    double lat_rad = fp - Q1 * (Q2 - Q3 + Q4);
    double Q5 = D;
    double Q6 = (1.0 + 2.0 * T1 + C1) * pow(D, 3) / 6.0;
    double Q7 = (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C1 * C1 + 8.0 * e2 + 24.0 * T1 *
                 T1) * pow(D, 5) / 120.0;
    double lon_rad = coord_deg_to_rad(lon_center) + (Q5 - Q6 + Q7) / cos_fp;
    geo->latitude = coord_normalize_latitude(coord_rad_to_deg(lat_rad));
    geo->longitude = coord_normalize_longitude(coord_rad_to_deg(lon_rad));
    geo->altitude = 0.0;
    geo->datum = utm->datum;
    return COORD_SUCCESS;
}

// Helper: get index of letter in MGRS alphabet (skip I and O)
static int get_mgrs_letter_index(char letter)
{
    if (letter < 'A' || letter > 'Z' || letter == 'I' || letter == 'O')
    {
        return -1;
    }
    int index = letter - 'A';
    if (letter > 'I')
    {
        index--;
    }
    if (letter > 'O')
    {
        index--;
    }
    return index;
}

// Helper: get MGRS letter from index (skip I and O)
static char get_mgrs_letter_from_index(int index)
{
    if (index < 0 || index > 23)
    {
        return 'Z';
    }
    char letter = 'A' + index;
    if (letter >= 'I')
    {
        letter++;
    }
    if (letter >= 'O')
    {
        letter++;
    }
    if (letter > 'Z')
    {
        letter = 'Z';
    }
    return letter;
}

// ==================== Key functions to fix MGRS conversion ====================
int coord_from_mgrs(CoordContext *ctx, const MGRSPoint *mgrs, GeoCoord *geo)
{
    if (!ctx || !mgrs || !geo)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_mgrs(mgrs))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    int zone = mgrs->zone;
    char band = mgrs->band;
    char col_letter = mgrs->square[0];
    char row_letter = mgrs->square[1];
    double easting = mgrs->easting;
    double northing = mgrs->northing;

    // Reverse compute MGRS 100k grid column letter (WGS84/modern datums)
    // Use the same 6-set system as coord_to_mgrs for reverse computation
    //
    // Logic to derive col_100k from column letter:
    // 1. Determine set (zone % 6)
    // 2. Get the start letter for the set
    // 3. Count from the start letter until the target column letter

    // Determine set (1-6)
    int set = zone % 6;
    if (set == 0)
        set = 6;

    // Get the start column letter for the set
    char col_origin;
    switch (set)
    {
        case 1:
        case 4:
            col_origin = 'A';
            break;
        case 2:
        case 5:
            col_origin = 'J';
            break;
        case 3:
        case 6:
            col_origin = 'S';
            break;
        default:
            col_origin = 'A';
            break;
    }

    // Count from start letter to target column letter, skipping I and O
    int col_100k = 0;
    int temp_code = col_origin;
    while (temp_code != col_letter)
    {
        temp_code++;
        if (temp_code == 'I' || temp_code == 'O')
            temp_code++;
        if (temp_code > 'Z')
        {
            temp_code = temp_code - 'Z' + 'A' - 1;
            if (temp_code == 'I' || temp_code == 'O')
                temp_code++;
        }
        col_100k++;
    }

    // Reverse compute row letter
    int row_idx = get_mgrs_letter_index(row_letter);
    if (row_idx < 0)
    {
        return COORD_ERROR_INVALID_COORD;
    }

    int row_offset = 0;
    int zone_parity = zone % 2;
    if (band >= 'N' && band <= 'X')
    {
        if (zone_parity == 0)
        {
            row_offset = 5;
        }
        else
        {
            row_offset = 0;
        }
    }
    else
    {
        if (zone_parity == 0)
        {
            row_offset = 0;
        }
        else
        {
            row_offset = 5;
        }
    }
    int row_100k = row_idx - row_offset;
    if (row_100k < 0)
    {
        row_100k += 20;
    }
    UTMPoint utm;
    utm.zone = zone;
    utm.band = band;
    utm.datum = mgrs->datum;
    utm.convergence = 0.0;
    utm.scale_factor = 0.9996;
    utm.easting = col_100k * 100000.0 + easting;
    if (band >= 'N' && band <= 'X')
    {
        utm.northing = row_100k * 100000.0 + northing;
    }
    else
    {
        double base_northing = row_100k * 100000.0 + northing;
        if (base_northing < 0)
        {
            base_northing += 2000000.0;
        }
        utm.northing = base_northing + 10000000.0;
    }
    if (!coord_validate_utm(&utm))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    int ret = coord_from_utm(ctx, &utm, geo);
    if (ret != COORD_SUCCESS)
    {
        return ret;
    }
    return COORD_SUCCESS;
}

int coord_to_mgrs(CoordContext *ctx, const GeoCoord *geo, MGRSPoint *mgrs)
{
    if (!ctx || !geo || !mgrs)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(geo))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    UTMPoint utm;
    int ret = coord_to_utm(ctx, geo, &utm);
    if (ret != COORD_SUCCESS)
    {
        return ret;
    }
    int zone = utm.zone;
    char band = utm.band;
    int col_100k = (int)(utm.easting / 100000.0);

    // Handle southern hemisphere northing
    // Southern hemisphere UTM northing = 10000000 + true northing (positive south of equator)
    // But MGRS uses true northing to compute 100km grid
    double utm_northing_for_mgrs = utm.northing;
    if (band < 'N')
    {
        // Southern hemisphere: remove false northing to get true northing (relative to southern origin)
        // Note: in the southern hemisphere, true northing may be negative (relative to equator)
        // But MGRS in the southern hemisphere is computed relative to the poleward direction
        utm_northing_for_mgrs -= 10000000.0;
    }

    // Compute 100km grid indices
    int row_100k = (int)(utm_northing_for_mgrs / 100000.0);

    // Southern hemisphere special case: row_100k may be negative
    // In the southern hemisphere, row_100k starts at 0 and increases southward
    // If negative, the coordinate is in high southern latitudes
    if (row_100k < 0)
    {
        // Polar region in southern hemisphere needs special handling
        // Convert negative index to positive index
        row_100k = 20 + (row_100k % 20);
    }

    // MGRS 100k grid column letter calculation (WGS84/modern datums)
    // Reference: GeoTrellis MGRS implementation (based on proj4js/mgrs)
    // Use 6-set cycle system instead of simple (zone-1)%3
    //
    // SET_ORIGIN_COLUMN_LETTERS = "AJSAJS" (6-set cycle)
    // - Zone % 6 = 1: A (ASCII 65)
    // - Zone % 6 = 2: J (ASCII 74)
    // - Zone % 6 = 3: S (ASCII 83)
    // - Zone % 6 = 4: A (ASCII 65)
    // - Zone % 6 = 5: J (ASCII 74)
    // - Zone % 6 = 0: S (ASCII 65)
    //
    // Formula: col_letter = origin + col_100k - 1
    // Skip I (73) and O (79)
    //
    // Example: Zone 50, col_100k=5
    // - set = 50 % 6 = 2
    // - origin = 'J' (ASCII 74)
    // - col = 'J' + 5 - 1 = 78 = 'N' ✓

    // Determine set (1-6)
    int set = zone % 6;
    if (set == 0)
        set = 6;

    // Get the start column letter for the set
    char col_origin;
    switch (set)
    {
        case 1:
        case 4:
            col_origin = 'A';
            break;
        case 2:
        case 5:
            col_origin = 'J';
            break;
        case 3:
        case 6:
            col_origin = 'S';
            break;
        default:
            col_origin = 'A';
            break;
    }

    // Compute column letter, skipping I and O
    int col_letter_code = col_origin + col_100k - 1;
    char col_letter;

    // Handle I and O skipping
    int temp_code = col_origin;
    for (int i = 1; i < col_100k; i++)
    {
        temp_code++;
        if (temp_code == 'I' || temp_code == 'O')
            temp_code++;
    }
    col_letter_code = temp_code;

    // Handle wrap-around (after Z back to A)
    if (col_letter_code > 'Z')
    {
        col_letter_code = col_letter_code - 'Z' + 'A' - 1;
        // Re-check I and O
        if (col_letter_code == 'I' || col_letter_code == 'O')
            col_letter_code++;
    }

    col_letter = (char)col_letter_code;

    // Compute row letter (using existing odd/even zone logic)
    int zone_parity = zone % 2;
    int row_offset = 0;
    if (band >= 'N' && band <= 'X')
    {
        if (zone_parity == 0)
        {
            row_offset = 5;
        }
        else
        {
            row_offset = 0;
        }
    }
    else
    {
        if (zone_parity == 0)
        {
            row_offset = 0;
        }
        else
        {
            row_offset = 5;
        }
    }
    int row_idx = (row_100k + row_offset) % 20;
    char row_letter = get_mgrs_letter_from_index(row_idx);
    mgrs->zone = zone;
    mgrs->band = band;
    mgrs->square[0] = col_letter;
    mgrs->square[1] = row_letter;
    mgrs->square[2] = '\0';
    mgrs->easting = fmod(utm.easting, 100000.0);
    mgrs->northing = fmod(utm_northing_for_mgrs, 100000.0);
    // Ensure northing is positive
    if (mgrs->northing < 0)
    {
        mgrs->northing += 100000.0;
    }
    mgrs->datum = utm.datum;
    return COORD_SUCCESS;
}

// Geographic coordinate to British Grid
int coord_to_british_grid(CoordContext *ctx, const GeoCoord *geo,
                          BritishGridPoint *bg)
{
    if (!ctx || !geo || !bg)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(geo))
    {
        return COORD_ERROR_INVALID_COORD;
    }

    // British National Grid must use OSGB36 datum and Airy 1830 ellipsoid
    // If input is not OSGB36, convert datum first
    GeoCoord osgb_geo;
    if (geo->datum != DATUM_OSGB36)
    {
        int ret = coord_convert_datum(ctx, geo, DATUM_OSGB36, &osgb_geo);
        if (ret != COORD_SUCCESS)
        {
            return ret;
        }
    }
    else
    {
        osgb_geo = *geo;
    }

    // Use OSGB36/Airy 1830 ellipsoid parameters
    double a = OSGB36_A;
    double f = OSGB36_F;
    double e2 = 2 * f - f * f;

    double lat_rad = coord_deg_to_rad(osgb_geo.latitude);
    double lon_rad = coord_deg_to_rad(osgb_geo.longitude);
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double tan_lat = sin_lat / cos_lat;
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double T = tan_lat * tan_lat;
    double C = e2 * cos_lat * cos_lat / (1.0 - e2);
    double A = (lon_rad - OSGB36_LON0) * cos_lat;
    // Compute M
    double M = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 /
                     256.0) * lat_rad
                    - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(
                        2.0 * lat_rad)
                    + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(4.0 * lat_rad)
                    - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * lat_rad));
    // Compute M0
    double M0 = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 /
                      256.0) * OSGB36_LAT0
                     - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(
                         2.0 * OSGB36_LAT0)
                     + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(
                         4.0 * OSGB36_LAT0)
                     - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * OSGB36_LAT0));
    double A2 = A * A;
    double A3 = A2 * A;
    double A4 = A3 * A;
    double A5 = A4 * A;
    double A6 = A5 * A;
    bg->easting = OSGB36_E0 + OSGB36_F0 * N * (A + (1.0 - T + C) * A3 / 6.0
                  + (5.0 - 18.0 * T + T * T + 72.0 * C - 58.0 * e2) * A5 / 120.0);
    bg->northing = OSGB36_N0 + OSGB36_F0 * (M - M0 + N * tan_lat *
                                            (A2 / 2.0 + (5.0 - T + 9.0 * C + 4.0 * C * C) * A4 / 24.0
                                                    + (61.0 - 58.0 * T + T * T + 600.0 * C - 330.0 * e2) * A6 / 720.0));

    // Compute British National Grid letters
    // British Grid uses a special 500km square letter system
    // Note: for coordinates outside the UK, letters are not standard
    // Here we use an extended cyclic computation

    // Compute 500km square index
    int e500k = (int)(bg->easting / 500000.0);
    int n500k = (int)(bg->northing / 500000.0);

    // Handle negative indices
    while (e500k < 0) e500k += 25;  // Ensure positive
    while (n500k < 0) n500k += 25;

    // Letters within 100km square
    int e100k = (int)(fmod(fabs(bg->easting), 500000.0) / 100000.0);
    int n100k = (int)(fmod(fabs(bg->northing), 500000.0) / 100000.0);

    // Use British Grid alphabet (skip I)
    const char *bg_letters = "ABCDEFGHJKLMNPQRSTUVWXYZ";

    // Compute easting letter (cycle alphabet)
    int e_idx = (e500k * 5 + e100k) % 25;
    bg->letters[0] = bg_letters[e_idx];

    // Compute northing letter (cycle alphabet)
    int n_idx = (n500k * 5 + n100k) % 25;
    bg->letters[1] = bg_letters[n_idx];

    bg->letters[2] = '\0';
    bg->datum = DATUM_OSGB36;  // British Grid always uses OSGB36 datum
    return COORD_SUCCESS;
}

int coord_from_british_grid(CoordContext *ctx, const BritishGridPoint *bg,
                            GeoCoord *geo)
{
    if (!ctx || !bg || !geo)
    {
        return COORD_ERROR_INVALID_INPUT;
    }

    // British National Grid (OSGB36) uses the Airy 1830 ellipsoid
    // First convert British Grid coordinates to OSGB36 lat/lon
    // Then convert to WGS84

    // Remove easting and northing offsets
    double E = bg->easting;
    double N = bg->northing;

    // Temporarily use WGS84 ellipsoid for simplified computation
    // Note: full OSGB36 conversion requires Airy 1830 ellipsoid parameters
    double a = OSGB36_A;
    double f = OSGB36_F;
    double e2 = 2 * f - f * f;

    // Compute central meridian and latitude of origin
    double lon0 = OSGB36_LON0;
    double E0 = OSGB36_E0;
    double N0 = OSGB36_N0;
    double F0 = OSGB36_F0;

    // Initial latitude estimate
    double M_prime = (N - N0) / F0;
    double mu_prime = M_prime / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0
                                   - 5.0 * e2 * e2 * e2 / 256.0));

    // Compute initial latitude
    double e1 = (1.0 - sqrt(1.0 - e2)) / (1.0 + sqrt(1.0 - e2));
    double J1 = 3.0 * e1 / 2.0 - 27.0 * e1 * e1 * e1 / 32.0;
    double J2 = 21.0 * e1 * e1 / 16.0 - 55.0 * e1 * e1 * e1 * e1 / 32.0;
    double J3 = 151.0 * e1 * e1 * e1 / 96.0;
    double J4 = 1097.0 * e1 * e1 * e1 * e1 / 512.0;

    double lat_prime = mu_prime + J1 * sin(2.0 * mu_prime) + J2 * sin(4.0 * mu_prime)
                      + J3 * sin(6.0 * mu_prime) + J4 * sin(8.0 * mu_prime);

    // Iteratively compute precise latitude and longitude
    double lat_rad = lat_prime;
    double lon_rad = lon0;

    for (int iter = 0; iter < 10; iter++)
    {
        double sin_lat = sin(lat_rad);
        double cos_lat = cos(lat_rad);
        double tan_lat = sin_lat / cos_lat;
        double sec_lat = 1.0 / cos_lat;

        double nu = a * F0 / sqrt(1.0 - e2 * sin_lat * sin_lat);
        double rho = a * F0 * (1.0 - e2) / pow(1.0 - e2 * sin_lat * sin_lat, 1.5);
        double eta2 = nu / rho - 1.0;

        double M = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0
                        - 5.0 * e2 * e2 * e2 / 256.0) * lat_rad
                       - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0
                          + 45.0 * e2 * e2 * e2 / 1024.0) * sin(2.0 * lat_rad)
                       + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(4.0 * lat_rad)
                       - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * lat_rad));

        double dM = M - M_prime;

        double delta_E = E - E0;
        double delta_E2 = delta_E * delta_E;
        double delta_E3 = delta_E2 * delta_E;
        double delta_E4 = delta_E3 * delta_E;

        double VII = tan_lat / (2.0 * rho * nu);
        double VIII = tan_lat / (24.0 * rho * nu * nu * nu)
                      * (5.0 + 3.0 * tan_lat * tan_lat + eta2 - 9.0 * tan_lat * tan_lat * eta2);
        double IX = tan_lat / (720.0 * rho * pow(nu, 5))
                    * (61.0 + 90.0 * tan_lat * tan_lat + 45.0 * pow(tan_lat, 4));
        double X = sec_lat / nu;
        double XI = sec_lat / (6.0 * nu * nu * nu) * (nu / rho + 2.0 * tan_lat * tan_lat);
        double XII = sec_lat / (120.0 * pow(nu, 5))
                     * (5.0 + 28.0 * tan_lat * tan_lat + 24.0 * pow(tan_lat, 4));
        double XIIA = sec_lat / (5040.0 * pow(nu, 7))
                      * (61.0 + 662.0 * tan_lat * tan_lat + 1320.0 * pow(tan_lat, 4)
                         + 720.0 * pow(tan_lat, 6));

        double lat_delta = VII * dM + VIII * dM * dM * dM + IX * dM * dM * dM * dM * dM;
        lat_rad -= lat_delta;

        double lon_delta = X * delta_E + XI * delta_E3 + XII * delta_E4 * delta_E
                          + XIIA * delta_E4 * delta_E3;
        lon_rad = lon0 + lon_delta;

        // Check convergence
        if (fabs(lat_delta) < 1e-12)
        {
            break;
        }
    }

    // Convert to degrees
    geo->latitude = coord_rad_to_deg(lat_rad);
    geo->longitude = coord_rad_to_deg(lon_rad);
    geo->altitude = 0.0;

    // Convert from OSGB36 to WGS84 (simplified approximation)
    // Approximate Helmert transform parameters for the UK
    double tx = -446.448;  // m
    double ty = 125.157;   // m
    double tz = -542.060;  // m
    double rx = -0.1502;   // arc-seconds
    double ry = -0.2470;   // arc-seconds
    double rz = -0.8421;   // arc-seconds
    double s = 20.4894;    // ppm

    // Use 7-parameter transform
    double lat_rad_osgb = lat_rad;
    double lon_rad_osgb = lon_rad;

    double sin_lat = sin(lat_rad_osgb);
    double cos_lat = cos(lat_rad_osgb);
    double sin_lon = sin(lon_rad_osgb);
    double cos_lon = cos(lon_rad_osgb);

    double nu = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double X = (nu + 0.0) * cos_lat * cos_lon;
    double Y = (nu + 0.0) * cos_lat * sin_lon;
    double Z = (nu * (1.0 - e2) + 0.0) * sin_lat;

    // Apply 7-parameter transform
    double rx_rad = rx * ARC_SEC_TO_RAD;
    double ry_rad = ry * ARC_SEC_TO_RAD;
    double rz_rad = rz * ARC_SEC_TO_RAD;
    double scale_factor = 1.0 + s * PPM_TO_SCALE;

    double X2 = tx + X * scale_factor + Y * rz_rad - Z * ry_rad;
    double Y2 = ty - X * rz_rad + Y * scale_factor + Z * rx_rad;
    double Z2 = tz + X * ry_rad - Y * rx_rad + Z * scale_factor;

    // Convert back to WGS84 geodetic coordinates
    const Ellipsoid *wgs84 = &ELLIPSOIDS[DATUM_WGS84];
    double p = sqrt(X2 * X2 + Y2 * Y2);
    double theta = atan2(Z2 * wgs84->a, p * wgs84->b);
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double lat_rad_wgs84 = atan2(Z2 + wgs84->ep2 * wgs84->b * sin_theta * sin_theta * sin_theta,
                                 p - wgs84->e2 * wgs84->a * cos_theta * cos_theta * cos_theta);
    double lon_rad_wgs84 = atan2(Y2, X2);

    geo->latitude = coord_rad_to_deg(lat_rad_wgs84);
    geo->longitude = coord_rad_to_deg(lon_rad_wgs84);
    geo->altitude = 0.0;
    geo->datum = DATUM_WGS84;

    return COORD_SUCCESS;
}

// Japan plane rectangular coordinate system zone parameters
static const struct
{
    int zone;
    double lat0;
    double lon0;
    double false_e;   // Easting offset (false easting)
    double false_n;   // Northing offset (false northing)
    double scale;
} japan_zones[] =
{
    {1,  33.0,  129.5,   0,       0,         0.9999},
    {2,  33.0,  131.0,   0,       0,         0.9999},
    {3,  36.0,  132.1667, 0,       0,         0.9999},
    {4,  33.0,  133.5,   0,       0,         0.9999},
    {5,  36.0,  134.3333, 0,       0,         0.9999},
    {6,  36.0,  136.0,   0,       0,         0.9999},
    {7,  36.0,  137.1667, 0,       0,         0.9999},
    {8,  36.0,  138.5,   0,       0,         0.9999},
    {9,  36.0,  139.8333, 0,       0,         0.9999},
    {10, 40.0,  140.8333, 0,       0,         0.9999},
    {11, 44.0,  140.25,  0,       0,         0.9999},
    {12, 44.0,  142.25,  0,       0,         0.9999},
    {13, 44.0,  144.25,  0,       0,         0.9999},
    {14, 26.0,  142.0,   0,       0,         0.9999},
    {15, 26.0,  127.5,   0,       0,         0.9999},
    {16, 26.0,  124.0,   0,       0,         0.9999},
    {17, 26.0,  131.0,   0,       0,         0.9999},
    {18, 20.0,  136.0,   0,       0,         0.9999},
    {19, 26.0,  154.0,   0,       0,         0.9999}
};

// Geographic coordinate to Japan Grid
int coord_to_japan_grid(CoordContext *ctx, const GeoCoord *geo,
                        JapanGridPoint *jg)
{
    if (!ctx || !geo || !jg)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(geo))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    // Convert to Tokyo datum
    GeoCoord tokyo_geo;
    int ret = coord_convert_datum(ctx, geo, DATUM_TOKYO, &tokyo_geo);
    if (ret != COORD_SUCCESS)
    {
        return ret;
    }

    double lat = tokyo_geo.latitude;
    double lon = tokyo_geo.longitude;

    // Select the appropriate zone based on lat/lon (nearest central meridian)
    // No geographic bounds; support any coordinates
    int zone_idx = -1;
    double min_dist = 1e308;
    for (size_t i = 0; i < sizeof(japan_zones) / sizeof(japan_zones[0]); i++)
    {
        double dx = (lon - japan_zones[i].lon0);
        double dy = (lat - japan_zones[i].lat0);
        double dist = dx * dx + dy * dy;
        if (dist < min_dist)
        {
            min_dist = dist;
            zone_idx = i;
        }
    }

    if (zone_idx < 0)
    {
        return COORD_ERROR_OUT_OF_RANGE;
    }

    // Get zone parameters
    double lat0 = japan_zones[zone_idx].lat0 * DEG_TO_RAD;
    double lon0 = japan_zones[zone_idx].lon0 * DEG_TO_RAD;
    double false_e = japan_zones[zone_idx].false_e;
    double false_n = japan_zones[zone_idx].false_n;
    double k0 = japan_zones[zone_idx].scale;

    double lat_rad = tokyo_geo.latitude * DEG_TO_RAD;
    double lon_rad = tokyo_geo.longitude * DEG_TO_RAD;

    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double tan_lat = sin_lat / cos_lat;

    // Bessel 1841 ellipsoid parameters
    double a = JAPAN_GRID_A;
    double f = JAPAN_GRID_F;
    double e2 = 2 * f - f * f;

    // Compute meridional arc length M
    double M = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 / 256.0) * lat_rad
                   - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(2.0 * lat_rad)
                   + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(4.0 * lat_rad)
                   - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * lat_rad));

    // Compute auxiliary parameters
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double T = tan_lat * tan_lat;
    double C = e2 * cos_lat * cos_lat / (1.0 - e2);
    double A = (lon_rad - lon0) * cos_lat;
    double A2 = A * A;
    double A3 = A2 * A;
    double A4 = A3 * A;
    double A5 = A4 * A;
    double A6 = A5 * A;

    // Compute X (northing)
    jg->x = k0 * (M + N * tan_lat * (A2 / 2.0 + (5.0 - T + 9.0 * C + 4.0 * C * C) * A4 / 24.0
               + (61.0 - 58.0 * T + T * T + 600.0 * C - 330.0 * e2) * A6 / 720.0)) + false_n;

    // Compute Y (easting)
    jg->y = k0 * N * (A + (1.0 - T + C) * A3 / 6.0
               + (5.0 - 18.0 * T + T * T + 72.0 * C - 58.0 * e2) * A5 / 120.0) + false_e;

    jg->zone = japan_zones[zone_idx].zone;
    jg->datum = DATUM_TOKYO;
    return COORD_SUCCESS;
}

int coord_from_japan_grid(CoordContext *ctx, const JapanGridPoint *jg,
                          GeoCoord *geo)
{
    if (!ctx || !jg || !geo)
    {
        return COORD_ERROR_INVALID_INPUT;
    }

    // Get Bessel 1841 ellipsoid parameters
    double a = JAPAN_GRID_A;
    double f = JAPAN_GRID_F;
    double e2 = 2 * f - f * f;

    // Lookup zone parameters (using externally defined japan_zones array)
    int zone_idx = -1;
    for (size_t i = 0; i < sizeof(japan_zones) / sizeof(japan_zones[0]); i++)
    {
        if (japan_zones[i].zone == jg->zone)
        {
            zone_idx = i;
            break;
        }
    }

    if (zone_idx < 0)
    {
        return COORD_ERROR_INVALID_INPUT;
    }

    double lat0 = japan_zones[zone_idx].lat0 * DEG_TO_RAD;
    double lon0 = japan_zones[zone_idx].lon0 * DEG_TO_RAD;
    double false_e = japan_zones[zone_idx].false_e;
    double false_n = japan_zones[zone_idx].false_n;
    double k0 = japan_zones[zone_idx].scale;

    // Gauss-Kruger inverse projection
    // Remove false offsets
    // Note: jg->x is northing, jg->y is easting
    double northing = jg->x - false_n;   // Northing
    double easting = jg->y - false_e;     // Easting

    // Compute auxiliary parameters
    // M is the meridional arc length, computed from northing
    double M = northing / k0;
    double mu = M / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0
                        - 5.0 * e2 * e2 * e2 / 256.0));

    // Compute footpoint latitude
    double e1 = (1.0 - sqrt(1.0 - e2)) / (1.0 + sqrt(1.0 - e2));
    double J1 = 3.0 * e1 / 2.0 - 27.0 * e1 * e1 * e1 / 32.0;
    double J2 = 21.0 * e1 * e1 / 16.0 - 55.0 * e1 * e1 * e1 * e1 / 32.0;
    double J3 = 151.0 * e1 * e1 * e1 / 96.0;
    double J4 = 1097.0 * e1 * e1 * e1 * e1 / 512.0;

    double fp = mu + J1 * sin(2.0 * mu) + J2 * sin(4.0 * mu)
                + J3 * sin(6.0 * mu) + J4 * sin(8.0 * mu);

    double sin_fp = sin(fp);
    double cos_fp = cos(fp);
    double tan_fp = sin_fp / cos_fp;

    double C1 = e2 * cos_fp * cos_fp / (1.0 - e2);
    double T1 = tan_fp * tan_fp;
    double R1 = a * (1.0 - e2) / pow(1.0 - e2 * sin_fp * sin_fp, 1.5);
    double N1 = a / sqrt(1.0 - e2 * sin_fp * sin_fp);

    // D computed from easting
    double D = easting / (N1 * k0);
    double D2 = D * D;
    double D3 = D2 * D;
    double D4 = D3 * D;
    double D5 = D4 * D;
    double D6 = D5 * D;

    // Compute latitude correction
    double Q1 = N1 * tan_fp / R1;
    double Q2 = 0.5 * D2;
    double Q3 = (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C1 * C1
                - 9.0 * e2) * D4 / 24.0;
    double Q4 = (61.0 + 90.0 * T1 + 298.0 * C1 + 45.0 * T1 * T1
                - 252.0 * e2 - 3.0 * C1 * C1) * D6 / 720.0;

    double lat_rad = fp - Q1 * (Q2 - Q3 + Q4);

    // Compute longitude correction
    double Q5 = D;
    double Q6 = (1.0 + 2.0 * T1 + C1) * D3 / 6.0;
    double Q7 = (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C1 * C1
                + 8.0 * e2 + 24.0 * T1 * T1) * D5 / 120.0;

    double lon_rad = lon0 + (Q5 - Q6 + Q7) / cos_fp;

    // Convert to degrees
    geo->latitude = coord_rad_to_deg(lat_rad);
    geo->longitude = coord_rad_to_deg(lon_rad);
    geo->altitude = 0.0;
    geo->datum = DATUM_TOKYO;

    // If WGS84 coordinates are needed, perform datum conversion
    if (geo->datum == DATUM_WGS84)
    {
        GeoCoord tokyo_coord = *geo;
        tokyo_coord.datum = DATUM_TOKYO;
        int ret = coord_convert_datum(ctx, &tokyo_coord, DATUM_WGS84, geo);
        if (ret != COORD_SUCCESS)
        {
            return ret;
        }
    }

    return COORD_SUCCESS;
}

// ==================== Datum conversion functions ====================
int coord_convert_datum(CoordContext *ctx, const GeoCoord *src,
                        MapDatum target_datum, GeoCoord *dst)
{
    if (!ctx || !src || !dst)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (src->datum == target_datum)
    {
        *dst = *src;
        return COORD_SUCCESS;
    }
    if (!coord_validate_point(src))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    // Get transform parameters
    DatumTransform *params = &ctx->transforms[src->datum][target_datum];
    if (params->dx == 0.0 && params->dy == 0.0 && params->dz == 0.0 &&
            params->rx == 0.0 && params->ry == 0.0 && params->rz == 0.0 &&
            params->scale == 0.0)
    {
        // No transform parameters; return directly
        *dst = *src;
        dst->datum = target_datum;
        return COORD_SUCCESS;
    }
    // Get source ellipsoid parameters
    const Ellipsoid *src_ell = &ELLIPSOIDS[src->datum];
    const Ellipsoid *dst_ell = &ELLIPSOIDS[target_datum];
    // Convert lat/lon to geocentric Cartesian coordinates
    double lat_rad = coord_deg_to_rad(src->latitude);
    double lon_rad = coord_deg_to_rad(src->longitude);
    double alt = src->altitude;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);
    double N = src_ell->a / sqrt(1.0 - src_ell->e2 * sin_lat * sin_lat);
    double X = (N + alt) * cos_lat * cos_lon;
    double Y = (N + alt) * cos_lat * sin_lon;
    double Z = (N * (1.0 - src_ell->e2) + alt) * sin_lat;
    // Apply 7-parameter transform
    double rx_rad = params->rx * ARC_SEC_TO_RAD;
    double ry_rad = params->ry * ARC_SEC_TO_RAD;
    double rz_rad = params->rz * ARC_SEC_TO_RAD;
    double scale_factor = 1.0 + params->scale * PPM_TO_SCALE;
    double X2 = params->dx + X * scale_factor + Y * rz_rad - Z * ry_rad;
    double Y2 = params->dy - X * rz_rad + Y * scale_factor + Z * rx_rad;
    double Z2 = params->dz + X * ry_rad - Y * rx_rad + Z * scale_factor;
    // Convert back to geodetic coordinates
    double p = sqrt(X2 * X2 + Y2 * Y2);
    double theta = atan2(Z2 * dst_ell->a, p * dst_ell->b);
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double lat_rad_out = atan2(Z2 + dst_ell->ep2 * dst_ell->b * sin_theta *
                               sin_theta * sin_theta,
                               p - dst_ell->e2 * dst_ell->a * cos_theta * cos_theta * cos_theta);
    double lon_rad_out = atan2(Y2, X2);
    double N2 = dst_ell->a / sqrt(1.0 - dst_ell->e2 * sin(lat_rad_out) * sin(
                                      lat_rad_out));
    double alt_out = p / cos(lat_rad_out) - N2;
    dst->latitude = coord_normalize_latitude(coord_rad_to_deg(lat_rad_out));
    dst->longitude = coord_normalize_longitude(coord_rad_to_deg(lon_rad_out));
    dst->altitude = alt_out;
    dst->datum = target_datum;
    return COORD_SUCCESS;
}

// ==================== Geodesic calculation functions ====================
int coord_distance(CoordContext *ctx, const GeoCoord *p1, const GeoCoord *p2,
                   double *distance, double *azi1, double *azi2)
{
    if (!ctx || !p1 || !p2)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(p1) || !coord_validate_point(p2))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    double s12, a1, a2;
    // If datums differ, convert
    if (p1->datum != p2->datum)
    {
        GeoCoord p2_same_datum;
        int ret = coord_convert_datum(ctx, p2, p1->datum, &p2_same_datum);
        if (ret != COORD_SUCCESS)
        {
            return ret;
        }
        geod_inverse(ctx->geod, p1->latitude, p1->longitude,
                     p2_same_datum.latitude, p2_same_datum.longitude,
                     &s12, &a1, &a2);
    }
    else
    {
        geod_inverse(ctx->geod, p1->latitude, p1->longitude,
                     p2->latitude, p2->longitude,
                     &s12, &a1, &a2);
    }
    if (distance)
    {
        *distance = s12;
    }
    if (azi1)
    {
        *azi1 = a1;
    }
    if (azi2)
    {
        *azi2 = a2;
    }
    return COORD_SUCCESS;
}

int coord_direct(CoordContext *ctx, const GeoCoord *start,
                 double distance, double azimuth, GeoCoord *end)
{
    if (!ctx || !start || !end)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(start))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    if (distance < 0.0)
    {
        return COORD_ERROR_OUT_OF_RANGE;
    }
    double lat2, lon2, azi2;
    geod_direct(ctx->geod, start->latitude, start->longitude,
                azimuth, distance, &lat2, &lon2, &azi2);
    end->latitude = coord_normalize_latitude(lat2);
    end->longitude = coord_normalize_longitude(lon2);
    end->altitude = 0.0;
    end->datum = start->datum;
    return COORD_SUCCESS;
}

int coord_inverse(CoordContext *ctx, const GeoCoord *p1, const GeoCoord *p2,
                  GeodesicResult *result)
{
    if (!ctx || !p1 || !p2 || !result)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(p1) || !coord_validate_point(p2))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    // If datums differ, convert
    if (p1->datum != p2->datum)
    {
        GeoCoord p2_same_datum;
        int ret = coord_convert_datum(ctx, p2, p1->datum, &p2_same_datum);
        if (ret != COORD_SUCCESS)
        {
            return ret;
        }
        geod_inverse(ctx->geod, p1->latitude, p1->longitude,
                     p2_same_datum.latitude, p2_same_datum.longitude,
                     &result->distance, &result->azimuth1, &result->azimuth2);
    }
    else
    {
        geod_inverse(ctx->geod, p1->latitude, p1->longitude,
                     p2->latitude, p2->longitude,
                     &result->distance, &result->azimuth1, &result->azimuth2);
    }
    return COORD_SUCCESS;
}

// ==================== Datum transform utility functions ====================
int coord_set_transform_params(CoordContext *ctx, MapDatum from, MapDatum to,
                               const DatumTransform *params)
{
    if (!ctx || from >= DATUM_MAX || to >= DATUM_MAX || !params)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    ctx->transforms[from][to] = *params;
    // Set inverse transform parameters (correct 7-parameter inverse)
    if (from != to)
    {
        // Correct formula for 7-parameter inverse transform:
        // Inverse transform requires inversion via rotation matrix
        // If forward: X2 = T + s*R*X1
        // Then inverse: X1 = (-1/s)*R^T*(X2-T) = (-1/s)*R^T*X2 + (1/s)*R^T*T
        //
        // For small-angle approximation (arc-seconds), inverse params are approximately:
        // Inverse translation = -R^T * T / (1+s)
        // Inverse rotation ≈ -rotation (but more accurate calc needed)

        double s = params->scale * PPM_TO_SCALE;

        // Build rotation matrix R
        // R = Rz(θz) * Ry(θy) * Rx(θx)
        // For small angles, can be approximated as:
        // R ≈ [ 1    -rz    ry
        //       rz     1   -rx
        //      -ry    rx     1 ]

        // Compute inverse transform parameters
        // Inverse scale factor
        ctx->transforms[to][from].scale = -params->scale;

        // Inverse rotation parameters (approximate, for small angles)
        ctx->transforms[to][from].rx = -params->rx;
        ctx->transforms[to][from].ry = -params->ry;
        ctx->transforms[to][from].rz = -params->rz;

        // Inverse translation parameters: T_back = -(dx,dy,dz) / (1+s)
        // More accurate calculation requires considering rotation matrix transpose
        double factor = 1.0 / (1.0 + s);
        ctx->transforms[to][from].dx = -params->dx * factor;
        ctx->transforms[to][from].dy = -params->dy * factor;
        ctx->transforms[to][from].dz = -params->dz * factor;

        // For small angles, add rotation correction term
        // Correction: T_back ≈ -(T + R×T) / (1+s)
        // Here we use first-order approximation
        double dx_corr = params->ry * params->dz - params->rz * params->dy;
        double dy_corr = params->rz * params->dx - params->rx * params->dz;
        double dz_corr = params->rx * params->dy - params->ry * params->dx;
        // Adjustment for converting arc-seconds to radians
        dx_corr *= ARC_SEC_TO_RAD;
        dy_corr *= ARC_SEC_TO_RAD;
        dz_corr *= ARC_SEC_TO_RAD;

        ctx->transforms[to][from].dx -= dx_corr * factor;
        ctx->transforms[to][from].dy -= dy_corr * factor;
        ctx->transforms[to][from].dz -= dz_corr * factor;
    }
    return COORD_SUCCESS;
}

int coord_get_transform_params(CoordContext *ctx, MapDatum from, MapDatum to,
                               DatumTransform *params)
{
    if (!ctx || from >= DATUM_MAX || to >= DATUM_MAX || !params)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    *params = ctx->transforms[from][to];
    return COORD_SUCCESS;
}

// ==================== Ellipsoid utility functions ====================
const Ellipsoid *coord_get_ellipsoid(MapDatum datum)
{
    if (datum >= DATUM_MAX)
    {
        return NULL;
    }
    return &ELLIPSOIDS[datum];
}

int coord_set_custom_ellipsoid(CoordContext *ctx, double a, double f)
{
    if (!ctx || a <= 0.0 || f <= 0.0)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    ctx->ellipsoid.a = a;
    ctx->ellipsoid.f = f;
    ctx->ellipsoid.b = a * (1.0 - f);
    ctx->ellipsoid.e2 = 2 * f - f * f;
    ctx->ellipsoid.ep2 = ctx->ellipsoid.e2 / (1.0 - ctx->ellipsoid.e2);
    ctx->ellipsoid.name = "Custom";
    // Reinitialize GeographicLib geodesic object
    geod_init(ctx->geod, a, f);
    return COORD_SUCCESS;
}

// ==================== Error handling functions ====================
const char *coord_get_error_string(int error_code)
{
    if (error_code < 0
            || error_code >= (int)(sizeof(ERROR_MESSAGES) / sizeof(ERROR_MESSAGES[0])))
    {
        return "Unknown error";
    }
    return ERROR_MESSAGES[error_code];
}

void coord_set_error_callback(void (*callback)(int, const char *))
{
    error_callback = callback;
}

// ==================== Main format conversion function ====================
int coord_convert(CoordContext *ctx, const GeoCoord *src,
                  CoordFormat target_format, MapDatum target_datum,
                  char *result_buffer, size_t buffer_size)
{
    if (!ctx || !src || !result_buffer)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    if (!coord_validate_point(src))
    {
        return COORD_ERROR_INVALID_COORD;
    }
    // Convert to target datum
    GeoCoord target_geo;
    if (src->datum != target_datum)
    {
        int ret = coord_convert_datum(ctx, src, target_datum, &target_geo);
        if (ret != COORD_SUCCESS)
        {
            return ret;
        }
    }
    else
    {
        target_geo = *src;
    }
    // Format according to target format
    switch (target_format)
    {
        case COORD_FORMAT_DD:
        case COORD_FORMAT_DMM:
        case COORD_FORMAT_DMS:
            return coord_format_to_string(&target_geo, target_format, result_buffer,
                                          buffer_size);
        case COORD_FORMAT_UTM:
        {
            UTMPoint utm;
            int ret = coord_to_utm(ctx, &target_geo, &utm);
            if (ret != COORD_SUCCESS)
            {
                return ret;
            }
            return coord_format_utm(&utm, result_buffer, buffer_size);
        }
        case COORD_FORMAT_MGRS:
        {
            MGRSPoint mgrs;
            int ret = coord_to_mgrs(ctx, &target_geo, &mgrs);
            if (ret != COORD_SUCCESS)
            {
                return ret;
            }
            return coord_format_mgrs(&mgrs, result_buffer, buffer_size);
        }
        case COORD_FORMAT_BRITISH_GRID:
        {
            BritishGridPoint bg;
            int ret = coord_to_british_grid(ctx, &target_geo, &bg);
            if (ret != COORD_SUCCESS)
            {
                return ret;
            }
            return coord_format_british_grid(&bg, result_buffer, buffer_size);
        }
        case COORD_FORMAT_JAPAN_GRID:
        {
            JapanGridPoint jg;
            int ret = coord_to_japan_grid(ctx, &target_geo, &jg);
            if (ret != COORD_SUCCESS)
            {
                return ret;
            }
            return coord_format_japan_grid(&jg, result_buffer, buffer_size);
        }
        default:
            return COORD_ERROR_UNSUPPORTED_FORMAT;
    }
}
