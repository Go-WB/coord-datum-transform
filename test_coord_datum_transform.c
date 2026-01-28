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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Error callback
static void error_handler(int code, const char *message)
{
    fprintf(stderr, "Coordinate transform error [%d]: %s\n", code, message);
}

// Utility function: compare doubles
static int compare_double(double a, double b, double epsilon)
{
    return fabs(a - b) < epsilon;
}

// Test context creation and destruction
void test_context_creation()
{
    printf("=== Test context creation and destruction ===\n");
    // Test creating WGS84 context
    CoordContext *ctx1 = coord_create_context(DATUM_WGS84);
    if (ctx1)
    {
        printf("Created WGS84 context successfully\n");
        // Test setting datum
        int ret = coord_set_datum(ctx1, DATUM_NAD83);
        if (ret == COORD_SUCCESS)
        {
            printf("Set datum to NAD83 successfully\n");
        }
        else
        {
            printf("Failed to set datum to NAD83: %s\n", coord_get_error_string(ret));
        }
        coord_destroy_context(ctx1);
        printf("Destroyed WGS84 context successfully\n");
    }
    else
    {
        printf("Failed to create WGS84 context\n");
    }
    // Test creating Tokyo datum context
    CoordContext *ctx2 = coord_create_context(DATUM_TOKYO);
    if (ctx2)
    {
        printf("Created Tokyo datum context successfully\n");
        coord_destroy_context(ctx2);
        printf("Destroyed Tokyo datum context successfully\n");
    }
    else
    {
        printf("Failed to create Tokyo datum context\n");
    }
    printf("\n");
}

// Test utility functions
void test_utility_functions()
{
    printf("=== Test utility functions ===\n");
    // Test latitude validation
    printf("Latitude validation:\n");
    printf("  90.0 -> %s\n", coord_is_valid_latitude(90.0) ? "valid" : "invalid");
    printf("  -90.0 -> %s\n", coord_is_valid_latitude(-90.0) ? "valid" : "invalid");
    printf("  91.0 -> %s\n", coord_is_valid_latitude(91.0) ? "valid" : "invalid");
    printf("  -91.0 -> %s\n", coord_is_valid_latitude(-91.0) ? "valid" : "invalid");
    // Test longitude validation
    printf("Longitude validation:\n");
    printf("  180.0 -> %s\n",
           coord_is_valid_longitude(180.0) ? "valid" : "invalid");
    printf("  -180.0 -> %s\n",
           coord_is_valid_longitude(-180.0) ? "valid" : "invalid");
    printf("  181.0 -> %s\n",
           coord_is_valid_longitude(181.0) ? "valid" : "invalid");
    printf("  -181.0 -> %s\n",
           coord_is_valid_longitude(-181.0) ? "valid" : "invalid");
    // Test UTM zone calculation
    printf("UTM zone calculation:\n");
    printf("  Shanghai (31.23, 121.47) -> zone %d\n", coord_get_utm_zone(121.47,
            31.23));
    printf("  New York (40.71, -74.01) -> zone %d\n", coord_get_utm_zone(-74.01,
            40.71));
    printf("  London (51.51, -0.13) -> zone %d\n", coord_get_utm_zone(-0.13,
            51.51));
    printf("  Sydney (-33.87, 151.21) -> zone %d\n", coord_get_utm_zone(151.21,
            -33.87));
    // Test UTM latitude band calculation
    printf("UTM latitude band:\n");
    printf("  31.23° -> band %c\n", coord_get_utm_band(31.23));
    printf("  40.71° -> band %c\n", coord_get_utm_band(40.71));
    printf("  -33.87° -> band %c\n", coord_get_utm_band(-33.87));
    printf("  51.51° -> band %c\n", coord_get_utm_band(51.51));
    // Test coordinate validation
    GeoCoord valid_coord = {31.23, 121.47, 0.0, DATUM_WGS84};
    GeoCoord invalid_coord = {100.0, 200.0, 0.0, DATUM_WGS84};
    printf("Coordinate validation:\n");
    printf("  Valid coordinate: %s\n",
           coord_validate_point(&valid_coord) ? "pass" : "fail");
    printf("  Invalid coordinate: %s\n",
           coord_validate_point(&invalid_coord) ? "pass" : "fail");
    printf("\n");
}

// Test coordinate parsing
void test_coord_parsing()
{
    printf("=== Test coordinate parsing ===\n");
    // Test DD parsing
    printf("DD parsing:\n");
    ParseResult result1 = coord_parse_string("31.230416°N, 121.473701°E",
                          COORD_FORMAT_DD, DATUM_WGS84);
    if (result1.success)
    {
        printf("  Parsed successfully: %.6f, %.6f (datum: %d)\n",
               result1.coord.latitude, result1.coord.longitude, result1.coord.datum);
    }
    else
    {
        printf("  Parse failed: %s\n", result1.error_msg);
    }
    // Test DMM parsing
    printf("DMM parsing:\n");
    ParseResult result2 = coord_parse_string("31°13.825'N, 121°28.422'E",
                          COORD_FORMAT_DMM, DATUM_WGS84);
    if (result2.success)
    {
        printf("  Parsed successfully: %.6f, %.6f\n", result2.coord.latitude,
               result2.coord.longitude);
    }
    else
    {
        printf("  Parse failed: %s\n", result2.error_msg);
    }
    // Test DMS parsing
    printf("DMS parsing:\n");
    ParseResult result3 = coord_parse_string("31°13'49.50\"N, 121°28'25.32\"E",
                          COORD_FORMAT_DMS, DATUM_WGS84);
    if (result3.success)
    {
        printf("  Parsed successfully: %.6f, %.6f\n", result3.coord.latitude,
               result3.coord.longitude);
    }
    else
    {
        printf("  Parse failed: %s\n", result3.error_msg);
    }
    // Test UTM parsing
    printf("UTM parsing (Zone 50N, 447600E 4419300N):\n");
    ParseResult result4 = coord_parse_string("50N 447600E 4419300N",
                          COORD_FORMAT_UTM, DATUM_WGS84);
    if (result4.success)
    {
        printf("  Parsed successfully: %.6f, %.6f\n", result4.coord.latitude,
               result4.coord.longitude);
    }
    else
    {
        printf("  Parse failed: %s\n", result4.error_msg);
    }
    // Adjusted MGRS parsing in the test program
    printf("MGRS parsing (51Q DQ 54634 56142):\n");
    ParseResult result5 = coord_parse_string("51Q DQ 54634 56142",
                          COORD_FORMAT_MGRS, DATUM_WGS84);
    if (result5.success)
    {
        printf("  Parsed successfully: %.6f, %.6f\n", result5.coord.latitude,
               result5.coord.longitude);
    }
    else
    {
        printf("  Parse failed: %s\n", result5.error_msg);
        // Debug: manually parse MGRS parameters
        int zone;
        char band;
        char square[3];
        double easting, northing;
        int count = sscanf("51Q SB 54634 56142", "%d%c%2s %lf %lf",
                           &zone, &band, square, &easting, &northing);
        if (count == 5)
        {
            printf("  Manual parse: zone=%d, band=%c, square=%s, easting=%.0f, northing=%.0f\n",
                   zone, band, square, easting, northing);
        }
    }
    // Test auto-parse
    printf("Auto-parse:\n");
    ParseResult result6 = coord_auto_parse("31.230416, 121.473701");
    if (result6.success)
    {
        printf("  Auto-parse success: format=%d, datum=%d, coord=(%.6f, %.6f)\n",
               result6.format, result6.datum, result6.coord.latitude, result6.coord.longitude);
    }
    else
    {
        printf("  Auto-parse failed: %s\n", result6.error_msg);
    }
    // Test auto-parse UTM
    printf("Auto-parse UTM (50N 447600 4419300):\n");
    ParseResult result7 = coord_auto_parse("50N 447600 4419300");
    if (result7.success)
    {
        printf("  Auto-parse success: format=%d, datum=%d, coord=(%.6f, %.6f)\n",
               result7.format, result7.datum, result7.coord.latitude, result7.coord.longitude);
    }
    else
    {
        printf("  Auto-parse failed: %s\n", result7.error_msg);
    }
    printf("\n");
}

// Test coordinate formatting
void test_coord_formatting()
{
    printf("=== Test coordinate formatting ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("Failed to create context; cannot run formatting tests\n");
        return;
    }
    GeoCoord test_coord = {31.230416, 121.473701, 0.0, DATUM_WGS84};
    char buffer[256];
    // Test DD format
    int ret = coord_format_dd(&test_coord, buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("DD format: %s\n", buffer);
    }
    else
    {
        printf("DD formatting failed: %s\n", coord_get_error_string(ret));
    }
    // Test DMM format
    ret = coord_format_dmm(&test_coord, buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("DMM format: %s\n", buffer);
    }
    else
    {
        printf("DMM formatting failed: %s\n", coord_get_error_string(ret));
    }
    // Test DMS format
    ret = coord_format_dms(&test_coord, buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("DMS format: %s\n", buffer);
    }
    else
    {
        printf("DMS formatting failed: %s\n", coord_get_error_string(ret));
    }
    // Test UTM format
    UTMPoint utm;
    ret = coord_to_utm(ctx, &test_coord, &utm);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_utm(&utm, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("UTM format: %s\n", buffer);
        }
        else
        {
            printf("UTM formatting failed: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("UTM conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test MGRS format
    MGRSPoint mgrs;
    ret = coord_to_mgrs(ctx, &test_coord, &mgrs);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_mgrs(&mgrs, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("MGRS format: %s\n", buffer);
        }
        else
        {
            printf("MGRS formatting failed: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("MGRS conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test British Grid format
    BritishGridPoint bg;
    ret = coord_to_british_grid(ctx, &test_coord, &bg);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_british_grid(&bg, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("British Grid format: %s\n", buffer);
        }
        else
        {
            printf("British Grid formatting failed: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("British Grid conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test Japan Grid format
    JapanGridPoint jg;
    ret = coord_to_japan_grid(ctx, &test_coord, &jg);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_japan_grid(&jg, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("Japan Grid format: %s\n", buffer);
        }
        else
        {
            printf("Japan Grid formatting failed: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("Japan Grid conversion failed: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// Test coordinate conversion
void test_coord_conversion()
{
    printf("=== Test coordinate conversion ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("Failed to create context; cannot run conversion tests\n");
        return;
    }
    GeoCoord test_coord = {31.230416, 121.473701, 0.0, DATUM_WGS84};
    char buffer[256];
    // Test main conversion function
    printf("Main conversion function:\n");
    // Test DD output
    int ret = coord_convert(ctx, &test_coord, COORD_FORMAT_DD, DATUM_WGS84, buffer,
                            sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  DD format: %s\n", buffer);
    }
    else
    {
        printf("  DD conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test UTM output
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_UTM, DATUM_UTM_GRID, buffer,
                        sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  UTM format: %s\n", buffer);
    }
    else
    {
        printf("  UTM conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test MGRS output
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_MGRS, DATUM_MGRS_GRID,
                        buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  MGRS format: %s\n", buffer);
    }
    else
    {
        printf("  MGRS conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test British Grid output
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_BRITISH_GRID, DATUM_ED50,
                        buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  British Grid format: %s\n", buffer);
    }
    else
    {
        printf("  British Grid conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test Japan Grid output
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_JAPAN_GRID, DATUM_TOKYO,
                        buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  Japan Grid format: %s\n", buffer);
    }
    else
    {
        printf("  Japan Grid conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test datum conversion
    printf("\nDatum conversion:\n");
    // WGS84 -> NAD83
    GeoCoord converted_nad83;
    ret = coord_convert_datum(ctx, &test_coord, DATUM_NAD83, &converted_nad83);
    if (ret == COORD_SUCCESS)
    {
        printf("  WGS84 -> NAD83: (%.6f, %.6f)\n",
               converted_nad83.latitude, converted_nad83.longitude);
        // Convert back to verify
        GeoCoord back_coord;
        ret = coord_convert_datum(ctx, &converted_nad83, DATUM_WGS84, &back_coord);
        if (ret == COORD_SUCCESS)
        {
            double lat_diff = fabs(back_coord.latitude - test_coord.latitude);
            double lon_diff = fabs(back_coord.longitude - test_coord.longitude);
            printf("  Round-trip error: Δlat=%.8f°, Δlon=%.8f°\n", lat_diff, lon_diff);
        }
    }
    else
    {
        printf("  WGS84 -> NAD83 conversion failed: %s\n", coord_get_error_string(ret));
    }
    // WGS84 -> NAD27
    GeoCoord converted_nad27;
    ret = coord_convert_datum(ctx, &test_coord, DATUM_NAD27, &converted_nad27);
    if (ret == COORD_SUCCESS)
    {
        printf("  WGS84 -> NAD27: (%.6f, %.6f)\n",
               converted_nad27.latitude, converted_nad27.longitude);
    }
    else
    {
        printf("  WGS84 -> NAD27 conversion failed: %s\n", coord_get_error_string(ret));
    }
    // WGS84 -> ED50
    GeoCoord converted_ed50;
    ret = coord_convert_datum(ctx, &test_coord, DATUM_ED50, &converted_ed50);
    if (ret == COORD_SUCCESS)
    {
        printf("  WGS84 -> ED50: (%.6f, %.6f)\n",
               converted_ed50.latitude, converted_ed50.longitude);
    }
    else
    {
        printf("  WGS84 -> ED50 conversion failed: %s\n", coord_get_error_string(ret));
    }
    // WGS84 -> Tokyo
    GeoCoord converted_tokyo;
    ret = coord_convert_datum(ctx, &test_coord, DATUM_TOKYO, &converted_tokyo);
    if (ret == COORD_SUCCESS)
    {
        printf("  WGS84 -> Tokyo: (%.6f, %.6f)\n",
               converted_tokyo.latitude, converted_tokyo.longitude);
    }
    else
    {
        printf("  WGS84 -> Tokyo conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test UTM conversion
    printf("\nUTM conversion:\n");
    UTMPoint utm;
    ret = coord_to_utm(ctx, &test_coord, &utm);
    if (ret == COORD_SUCCESS)
    {
        printf("  Geographic -> UTM: %d%c %.3fE %.3fN\n",
               utm.zone, utm.band, utm.easting, utm.northing);
        // Convert back
        GeoCoord geo_back;
        ret = coord_from_utm(ctx, &utm, &geo_back);
        if (ret == COORD_SUCCESS)
        {
            double lat_diff = fabs(geo_back.latitude - test_coord.latitude);
            double lon_diff = fabs(geo_back.longitude - test_coord.longitude);
            printf("  UTM -> Geographic: (%.6f, %.6f), error: Δlat=%.8f°, Δlon=%.8f°\n",
                   geo_back.latitude, geo_back.longitude, lat_diff, lon_diff);
        }
    }
    else
    {
        printf("  UTM conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test MGRS conversion
    printf("\nMGRS conversion:\n");
    MGRSPoint mgrs;
    ret = coord_to_mgrs(ctx, &test_coord, &mgrs);
    if (ret == COORD_SUCCESS)
    {
        printf("  Geographic -> MGRS: %d%c %s %05.0f %05.0f\n",
               mgrs.zone, mgrs.band, mgrs.square, mgrs.easting, mgrs.northing);
        // Convert back
        GeoCoord geo_back_mgrs;
        ret = coord_from_mgrs(ctx, &mgrs, &geo_back_mgrs);
        if (ret == COORD_SUCCESS)
        {
            double lat_diff = fabs(geo_back_mgrs.latitude - test_coord.latitude);
            double lon_diff = fabs(geo_back_mgrs.longitude - test_coord.longitude);
            printf("  MGRS -> Geographic: (%.6f, %.6f), error: Δlat=%.8f°, Δlon=%.8f°\n",
                   geo_back_mgrs.latitude, geo_back_mgrs.longitude, lat_diff, lon_diff);
        }
    }
    else
    {
        printf("  MGRS conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test British Grid conversion
    printf("\nBritish Grid conversion:\n");
    BritishGridPoint bg;
    ret = coord_to_british_grid(ctx, &test_coord, &bg);
    if (ret == COORD_SUCCESS)
    {
        printf("  Geographic -> British Grid: %s %.0f %.0f\n",
               bg.letters, bg.easting, bg.northing);
    }
    else
    {
        printf("  British Grid conversion failed: %s\n", coord_get_error_string(ret));
    }
    // Test Japan Grid conversion
    printf("\nJapan Grid conversion:\n");
    JapanGridPoint jg;
    ret = coord_to_japan_grid(ctx, &test_coord, &jg);
    if (ret == COORD_SUCCESS)
    {
        printf("  Geographic -> Japan Grid: Zone %d: %.3f, %.3f\n",
               jg.zone, jg.x, jg.y);
    }
    else
    {
        printf("  Japan Grid conversion failed: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// Test geodesic calculations
void test_geodesic_calculation()
{
    printf("=== Test geodesic calculations ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("Failed to create context; cannot run geodesic tests\n");
        return;
    }
    // Define two test points: Shanghai and Beijing
    GeoCoord shanghai = {31.230416, 121.473701, 0.0, DATUM_WGS84};
    GeoCoord beijing = {39.904211, 116.407394, 0.0, DATUM_WGS84};
    // Test distance calculation
    double distance, azi1, azi2;
    int ret = coord_distance(ctx, &shanghai, &beijing, &distance, &azi1, &azi2);
    if (ret == COORD_SUCCESS)
    {
        printf("Shanghai to Beijing:\n");
        printf("  Distance: %.2f m (approx %.2f km)\n", distance, distance / 1000.0);
        printf("  Forward azimuth: %.2f°\n", azi1);
        printf("  Reverse azimuth: %.2f°\n", azi2);
    }
    else
    {
        printf("Distance calculation failed: %s\n", coord_get_error_string(ret));
    }
    // Test direct calculation
    printf("\nDirect calculation:\n");
    double test_distance = 100000.0; // 100 km
    double test_azimuth = 45.0; // 45-degree direction
    GeoCoord end_point;
    ret = coord_direct(ctx, &shanghai, test_distance, test_azimuth, &end_point);
    if (ret == COORD_SUCCESS)
    {
        printf("  From Shanghai, heading %.1f° for %.0f m:\n", test_azimuth,
               test_distance);
        printf("  Reached: (%.6f, %.6f)\n", end_point.latitude, end_point.longitude);
    }
    else
    {
        printf("  Direct calculation failed: %s\n", coord_get_error_string(ret));
    }
    // Test inverse calculation
    printf("\nInverse calculation:\n");
    GeodesicResult result;
    ret = coord_inverse(ctx, &shanghai, &beijing, &result);
    if (ret == COORD_SUCCESS)
    {
        printf("  Inverse calculation Shanghai to Beijing:\n");
        printf("  Distance: %.2f m\n", result.distance);
        printf("  Forward azimuth: %.2f°\n", result.azimuth1);
        printf("  Reverse azimuth: %.2f°\n", result.azimuth2);
    }
    else
    {
        printf("  Inverse calculation failed: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// Test datum transform tools
void test_datum_tools()
{
    printf("=== Test datum transform tools ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("Failed to create context\n");
        return;
    }
    // Test getting ellipsoid parameters
    const Ellipsoid *ellipsoid = coord_get_ellipsoid(DATUM_WGS84);
    if (ellipsoid)
    {
        printf("WGS84 ellipsoid parameters:\n");
        printf("  Semi-major axis: %.3f m\n", ellipsoid->a);
        printf("  Flattening: 1/%.9f\n", 1.0 / ellipsoid->f);
        printf("  Semi-minor axis: %.3f m\n", ellipsoid->b);
        printf("  Name: %s\n", ellipsoid->name);
    }
    else
    {
        printf("Failed to get ellipsoid parameters\n");
    }
    // Test NAD83 ellipsoid parameters
    ellipsoid = coord_get_ellipsoid(DATUM_NAD83);
    if (ellipsoid)
    {
        printf("\nNAD83 ellipsoid parameters:\n");
        printf("  Semi-major axis: %.3f m\n", ellipsoid->a);
        printf("  Flattening: 1/%.9f\n", 1.0 / ellipsoid->f);
        printf("  Name: %s\n", ellipsoid->name);
    }
    // Test setting and getting transform parameters
    DatumTransform transform = {100.0, 200.0, 300.0, 1.0, 2.0, 3.0, 10.0};
    int ret = coord_set_transform_params(ctx, DATUM_WGS84, DATUM_TOKYO, &transform);
    if (ret == COORD_SUCCESS)
    {
        printf("\nSet transform parameters successfully\n");
        DatumTransform get_transform;
        ret = coord_get_transform_params(ctx, DATUM_WGS84, DATUM_TOKYO, &get_transform);
        if (ret == COORD_SUCCESS)
        {
            if (compare_double(transform.dx, get_transform.dx, 0.001) &&
                    compare_double(transform.dy, get_transform.dy, 0.001) &&
                    compare_double(transform.dz, get_transform.dz, 0.001))
            {
                printf("Got transform parameters successfully\n");
            }
            else
            {
                printf("Retrieved transform parameters do not match\n");
            }
        }
        else
        {
            printf("Failed to get transform parameters: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("Failed to set transform parameters: %s\n", coord_get_error_string(ret));
    }
    // Test custom ellipsoid
    ret = coord_set_custom_ellipsoid(ctx, 6371000.0, 1.0 / 298.3);
    if (ret == COORD_SUCCESS)
    {
        printf("Set custom ellipsoid successfully\n");
    }
    else
    {
        printf("Failed to set custom ellipsoid: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// Test error handling
void test_error_handling()
{
    printf("=== Test error handling ===\n");
    // Test invalid input
    CoordContext *ctx = NULL;
    GeoCoord invalid_coord = {100.0, 200.0, 0.0, DATUM_WGS84};
    char buffer[256];
    int ret = coord_convert(ctx, &invalid_coord, COORD_FORMAT_DD, DATUM_WGS84,
                            buffer, sizeof(buffer));
    printf("Null context conversion test: %s (expected: invalid input)\n",
           ret == COORD_ERROR_INVALID_INPUT ? "pass" : "fail");
    // Test error message retrieval
    printf("\nError message test:\n");
    for (int i = 0; i <= 10; i++)
    {
        const char *error_msg = coord_get_error_string(i);
        printf("  Error code %d: %s\n", i, error_msg);
    }
    printf("\n");
}

// Comprehensive tests
void test_comprehensive()
{
    printf("=== Comprehensive tests ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("Failed to create context\n");
        return;
    }
    // Define multiple test points
    struct
    {
        const char *name;
        double lat, lon;
    } test_points[] =
    {
        {"Shanghai", 31.230416, 121.473701},
        {"Beijing", 39.904211, 116.407394},
        {"New York", 40.712776, -74.005974},
        {"London", 51.507351, -0.127758},
        {"Sydney", -33.868820, 151.209290},
        {"Tokyo", 35.689487, 139.691711},
        {"Paris", 48.856614, 2.352222}
    };
    int num_points = sizeof(test_points) / sizeof(test_points[0]);
    // Test format conversion for each point
    for (int i = 0; i < num_points; i++)
    {
        printf("%s coordinate conversion:\n", test_points[i].name);
        GeoCoord coord = {test_points[i].lat, test_points[i].lon, 0.0, DATUM_WGS84};
        char buffer[256];
        // Test multiple formats
        CoordFormat formats[] = {COORD_FORMAT_DD, COORD_FORMAT_DMM, COORD_FORMAT_DMS};
        const char *format_names[] = {"DD", "DMM", "DMS"};
        for (int j = 0; j < 3; j++)
        {
            int ret = coord_convert(ctx, &coord, formats[j], DATUM_WGS84, buffer,
                                    sizeof(buffer));
            if (ret == COORD_SUCCESS)
            {
                printf("  %s: %s\n", format_names[j], buffer);
            }
            else
            {
                printf("  %s format failed: %s\n", format_names[j], coord_get_error_string(ret));
            }
        }
        // Test UTM conversion
        UTMPoint utm;
        int ret = coord_to_utm(ctx, &coord, &utm);
        if (ret == COORD_SUCCESS)
        {
            printf("  UTM: zone %d%c\n", utm.zone, utm.band);
        }
        else
        {
            printf("  UTM conversion failed: %s\n", coord_get_error_string(ret));
        }
        // Test MGRS conversion
        MGRSPoint mgrs;
        ret = coord_to_mgrs(ctx, &coord, &mgrs);
        if (ret == COORD_SUCCESS)
        {
            printf("  MGRS: zone %d%c\n", mgrs.zone, mgrs.band);
        }
        else
        {
            printf("  MGRS conversion failed: %s\n", coord_get_error_string(ret));
        }
        printf("\n");
    }
    // Test point-to-point distance
    printf("Point-to-point distance:\n");
    for (int i = 0; i < num_points; i++)
    {
        for (int j = i + 1; j < num_points; j++)
        {
            GeoCoord p1 = {test_points[i].lat, test_points[i].lon, 0.0, DATUM_WGS84};
            GeoCoord p2 = {test_points[j].lat, test_points[j].lon, 0.0, DATUM_WGS84};
            double distance;
            int ret = coord_distance(ctx, &p1, &p2, &distance, NULL, NULL);
            if (ret == COORD_SUCCESS)
            {
                printf("  %s -> %s: %.2f km\n",
                       test_points[i].name, test_points[j].name, distance / 1000.0);
            }
        }
    }
    printf("\nMGRS coordinate conversion test:\n");
    // Test MGRS coordinate conversion
    struct
    {
        const char *name;
        double lat, lon;
        const char *expected_mgrs;
    } mgrs_test_points[] =
    {
        {"Shanghai", 31.230416, 121.473701, "51R"},
        {"Beijing", 39.904211, 116.407394, "50S"},
        {"Sydney", -33.868820, 151.209290, "56H"}
    };
    for (int i = 0; i < 3; i++)
    {
        GeoCoord coord = {mgrs_test_points[i].lat, mgrs_test_points[i].lon, 0.0, DATUM_WGS84};
        MGRSPoint mgrs;
        int ret = coord_to_mgrs(ctx, &coord, &mgrs);
        if (ret == COORD_SUCCESS)
        {
            char buffer[256];
            ret = coord_format_mgrs(&mgrs, buffer, sizeof(buffer));
            if (ret == COORD_SUCCESS)
            {
                printf("  %s: %s (expected zone: %s)\n",
                       mgrs_test_points[i].name, buffer, mgrs_test_points[i].expected_mgrs);
                // Convert back
                GeoCoord back_coord;
                ret = coord_from_mgrs(ctx, &mgrs, &back_coord);
                if (ret == COORD_SUCCESS)
                {
                    double lat_diff = fabs(back_coord.latitude - mgrs_test_points[i].lat);
                    double lon_diff = fabs(back_coord.longitude - mgrs_test_points[i].lon);
                    printf("    Round-trip error: Δlat=%.6f°, Δlon=%.6f°\n", lat_diff, lon_diff);
                }
            }
        }
        else
        {
            printf("  %s MGRS conversion failed: %s\n", mgrs_test_points[i].name,
                   coord_get_error_string(ret));
        }
    }
    coord_destroy_context(ctx);
    printf("\n");
}

int main()
{
    printf("=== Coordinate Transformation System Enhanced Test Program ===\n\n");
    // Set error callback
    coord_set_error_callback(error_handler);
    // Run all tests
    test_context_creation();
    test_utility_functions();
    test_coord_parsing();
    test_coord_formatting();
    test_coord_conversion();
    test_geodesic_calculation();
    test_datum_tools();
    test_error_handling();
    test_comprehensive();
    printf("=== All tests completed ===\n");
    return 0;
}
