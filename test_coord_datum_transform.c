#include "coord_transform.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// 错误回调函数
static void error_handler(int code, const char *message)
{
    fprintf(stderr, "坐标转换错误 [%d]: %s\n", code, message);
}

// 工具函数：比较浮点数
static int compare_double(double a, double b, double epsilon)
{
    return fabs(a - b) < epsilon;
}

// 测试上下文创建和销毁
void test_context_creation()
{
    printf("=== 测试上下文创建和销毁 ===\n");
    // 测试创建WGS84上下文
    CoordContext *ctx1 = coord_create_context(DATUM_WGS84);
    if (ctx1)
    {
        printf("创建WGS84上下文成功\n");
        // 测试设置基准
        int ret = coord_set_datum(ctx1, DATUM_NAD83);
        if (ret == COORD_SUCCESS)
        {
            printf("设置基准到NAD83成功\n");
        }
        else
        {
            printf("设置基准到NAD83失败: %s\n", coord_get_error_string(ret));
        }
        coord_destroy_context(ctx1);
        printf("销毁WGS84上下文成功\n");
    }
    else
    {
        printf("创建WGS84上下文失败\n");
    }
    // 测试创建Tokyo基准上下文
    CoordContext *ctx2 = coord_create_context(DATUM_TOKYO);
    if (ctx2)
    {
        printf("创建Tokyo基准上下文成功\n");
        coord_destroy_context(ctx2);
        printf("销毁Tokyo基准上下文成功\n");
    }
    else
    {
        printf("创建Tokyo基准上下文失败\n");
    }
    printf("\n");
}

// 测试工具函数
void test_utility_functions()
{
    printf("=== 测试工具函数 ===\n");
    // 测试纬度验证
    printf("纬度验证测试:\n");
    printf("  90.0 -> %s\n", coord_is_valid_latitude(90.0) ? "有效" : "无效");
    printf("  -90.0 -> %s\n", coord_is_valid_latitude(-90.0) ? "有效" : "无效");
    printf("  91.0 -> %s\n", coord_is_valid_latitude(91.0) ? "有效" : "无效");
    printf("  -91.0 -> %s\n", coord_is_valid_latitude(-91.0) ? "有效" : "无效");
    // 测试经度验证
    printf("经度验证测试:\n");
    printf("  180.0 -> %s\n",
           coord_is_valid_longitude(180.0) ? "有效" : "无效");
    printf("  -180.0 -> %s\n",
           coord_is_valid_longitude(-180.0) ? "有效" : "无效");
    printf("  181.0 -> %s\n",
           coord_is_valid_longitude(181.0) ? "有效" : "无效");
    printf("  -181.0 -> %s\n",
           coord_is_valid_longitude(-181.0) ? "有效" : "无效");
    // 测试UTM区域计算
    printf("UTM区域计算测试:\n");
    printf("  上海(31.23, 121.47) -> 区域 %d\n", coord_get_utm_zone(121.47,
            31.23));
    printf("  纽约(40.71, -74.01) -> 区域 %d\n", coord_get_utm_zone(-74.01,
            40.71));
    printf("  伦敦(51.51, -0.13) -> 区域 %d\n", coord_get_utm_zone(-0.13,
            51.51));
    printf("  悉尼(-33.87, 151.21) -> 区域 %d\n", coord_get_utm_zone(151.21,
            -33.87));
    // 测试UTM纬度带计算
    printf("UTM纬度带测试:\n");
    printf("  31.23° -> 带 %c\n", coord_get_utm_band(31.23));
    printf("  40.71° -> 带 %c\n", coord_get_utm_band(40.71));
    printf("  -33.87° -> 带 %c\n", coord_get_utm_band(-33.87));
    printf("  51.51° -> 带 %c\n", coord_get_utm_band(51.51));
    // 测试坐标验证
    GeoCoord valid_coord = {31.23, 121.47, 0.0, DATUM_WGS84};
    GeoCoord invalid_coord = {100.0, 200.0, 0.0, DATUM_WGS84};
    printf("坐标验证测试:\n");
    printf("  有效坐标: %s\n",
           coord_validate_point(&valid_coord) ? "通过" : "失败");
    printf("  无效坐标: %s\n",
           coord_validate_point(&invalid_coord) ? "通过" : "失败");
    printf("\n");
}

// 测试坐标解析
void test_coord_parsing()
{
    printf("=== 测试坐标解析 ===\n");
    // 测试DD格式解析
    printf("DD格式解析:\n");
    ParseResult result1 = coord_parse_string("31.230416°N, 121.473701°E",
                          COORD_FORMAT_DD, DATUM_WGS84);
    if (result1.success)
    {
        printf("  解析成功: %.6f, %.6f (基准: %d)\n",
               result1.coord.latitude, result1.coord.longitude, result1.coord.datum);
    }
    else
    {
        printf("  解析失败: %s\n", result1.error_msg);
    }
    // 测试DMM格式解析
    printf("DMM格式解析:\n");
    ParseResult result2 = coord_parse_string("31°13.825'N, 121°28.422'E",
                          COORD_FORMAT_DMM, DATUM_WGS84);
    if (result2.success)
    {
        printf("  解析成功: %.6f, %.6f\n", result2.coord.latitude,
               result2.coord.longitude);
    }
    else
    {
        printf("  解析失败: %s\n", result2.error_msg);
    }
    // 测试DMS格式解析
    printf("DMS格式解析:\n");
    ParseResult result3 = coord_parse_string("31°13'49.50\"N, 121°28'25.32\"E",
                          COORD_FORMAT_DMS, DATUM_WGS84);
    if (result3.success)
    {
        printf("  解析成功: %.6f, %.6f\n", result3.coord.latitude,
               result3.coord.longitude);
    }
    else
    {
        printf("  解析失败: %s\n", result3.error_msg);
    }
    // 测试UTM格式解析
    printf("UTM格式解析 (Zone 50N, 447600E 4419300N):\n");
    ParseResult result4 = coord_parse_string("50N 447600E 4419300N",
                          COORD_FORMAT_UTM, DATUM_WGS84);
    if (result4.success)
    {
        printf("  解析成功: %.6f, %.6f\n", result4.coord.latitude,
               result4.coord.longitude);
    }
    else
    {
        printf("  解析失败: %s\n", result4.error_msg);
    }
    // 在测试程序中修改MGRS解析部分
    printf("MGRS格式解析 (51Q DQ 54634 56142):\n");
    ParseResult result5 = coord_parse_string("51Q DQ 54634 56142",
                          COORD_FORMAT_MGRS, DATUM_WGS84);
    if (result5.success)
    {
        printf("  解析成功: %.6f, %.6f\n", result5.coord.latitude,
               result5.coord.longitude);
    }
    else
    {
        printf("  解析失败: %s\n", result5.error_msg);
        // 添加调试：手动解析MGRS参数
        int zone;
        char band;
        char square[3];
        double easting, northing;
        int count = sscanf("51Q SB 54634 56142", "%d%c%2s %lf %lf",
                           &zone, &band, square, &easting, &northing);
        if (count == 5)
        {
            printf("  手动解析: zone=%d, band=%c, square=%s, easting=%.0f, northing=%.0f\n",
                   zone, band, square, easting, northing);
        }
    }
    // 测试自动解析
    printf("自动解析:\n");
    ParseResult result6 = coord_auto_parse("31.230416, 121.473701");
    if (result6.success)
    {
        printf("  自动解析成功: 格式=%d, 基准=%d, 坐标=(%.6f, %.6f)\n",
               result6.format, result6.datum, result6.coord.latitude, result6.coord.longitude);
    }
    else
    {
        printf("  自动解析失败: %s\n", result6.error_msg);
    }
    // 测试自动解析UTM
    printf("自动解析UTM (50N 447600 4419300):\n");
    ParseResult result7 = coord_auto_parse("50N 447600 4419300");
    if (result7.success)
    {
        printf("  自动解析成功: 格式=%d, 基准=%d, 坐标=(%.6f, %.6f)\n",
               result7.format, result7.datum, result7.coord.latitude, result7.coord.longitude);
    }
    else
    {
        printf("  自动解析失败: %s\n", result7.error_msg);
    }
    printf("\n");
}

// 测试坐标格式化
void test_coord_formatting()
{
    printf("=== 测试坐标格式化 ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("创建上下文失败，无法进行格式化测试\n");
        return;
    }
    GeoCoord test_coord = {31.230416, 121.473701, 0.0, DATUM_WGS84};
    char buffer[256];
    // 测试DD格式
    int ret = coord_format_dd(&test_coord, buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("DD格式: %s\n", buffer);
    }
    else
    {
        printf("DD格式化失败: %s\n", coord_get_error_string(ret));
    }
    // 测试DMM格式
    ret = coord_format_dmm(&test_coord, buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("DMM格式: %s\n", buffer);
    }
    else
    {
        printf("DMM格式化失败: %s\n", coord_get_error_string(ret));
    }
    // 测试DMS格式
    ret = coord_format_dms(&test_coord, buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("DMS格式: %s\n", buffer);
    }
    else
    {
        printf("DMS格式化失败: %s\n", coord_get_error_string(ret));
    }
    // 测试UTM格式
    UTMPoint utm;
    ret = coord_to_utm(ctx, &test_coord, &utm);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_utm(&utm, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("UTM格式: %s\n", buffer);
        }
        else
        {
            printf("UTM格式化失败: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("UTM转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试MGRS格式
    MGRSPoint mgrs;
    ret = coord_to_mgrs(ctx, &test_coord, &mgrs);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_mgrs(&mgrs, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("MGRS格式: %s\n", buffer);
        }
        else
        {
            printf("MGRS格式化失败: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("MGRS转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试British Grid格式
    BritishGridPoint bg;
    ret = coord_to_british_grid(ctx, &test_coord, &bg);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_british_grid(&bg, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("British Grid格式: %s\n", buffer);
        }
        else
        {
            printf("British Grid格式化失败: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("British Grid转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试Japan Grid格式
    JapanGridPoint jg;
    ret = coord_to_japan_grid(ctx, &test_coord, &jg);
    if (ret == COORD_SUCCESS)
    {
        ret = coord_format_japan_grid(&jg, buffer, sizeof(buffer));
        if (ret == COORD_SUCCESS)
        {
            printf("Japan Grid格式: %s\n", buffer);
        }
        else
        {
            printf("Japan Grid格式化失败: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("Japan Grid转换失败: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// 测试坐标转换
void test_coord_conversion()
{
    printf("=== 测试坐标转换 ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("创建上下文失败，无法进行转换测试\n");
        return;
    }
    GeoCoord test_coord = {31.230416, 121.473701, 0.0, DATUM_WGS84};
    char buffer[256];
    // 测试主转换函数
    printf("主转换函数测试:\n");
    // 测试DD格式输出
    int ret = coord_convert(ctx, &test_coord, COORD_FORMAT_DD, DATUM_WGS84, buffer,
                            sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  DD格式: %s\n", buffer);
    }
    else
    {
        printf("  DD转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试UTM格式输出
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_UTM, DATUM_UTM_GRID, buffer,
                        sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  UTM格式: %s\n", buffer);
    }
    else
    {
        printf("  UTM转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试MGRS格式输出
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_MGRS, DATUM_MGRS_GRID,
                        buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  MGRS格式: %s\n", buffer);
    }
    else
    {
        printf("  MGRS转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试British Grid格式输出
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_BRITISH_GRID, DATUM_ED50,
                        buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  British Grid格式: %s\n", buffer);
    }
    else
    {
        printf("  British Grid转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试Japan Grid格式输出
    ret = coord_convert(ctx, &test_coord, COORD_FORMAT_JAPAN_GRID, DATUM_TOKYO,
                        buffer, sizeof(buffer));
    if (ret == COORD_SUCCESS)
    {
        printf("  Japan Grid格式: %s\n", buffer);
    }
    else
    {
        printf("  Japan Grid转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试基准转换
    printf("\n基准转换测试:\n");
    // WGS84 -> NAD83
    GeoCoord converted_nad83;
    ret = coord_convert_datum(ctx, &test_coord, DATUM_NAD83, &converted_nad83);
    if (ret == COORD_SUCCESS)
    {
        printf("  WGS84 -> NAD83: (%.6f, %.6f)\n",
               converted_nad83.latitude, converted_nad83.longitude);
        // 转换回来验证
        GeoCoord back_coord;
        ret = coord_convert_datum(ctx, &converted_nad83, DATUM_WGS84, &back_coord);
        if (ret == COORD_SUCCESS)
        {
            double lat_diff = fabs(back_coord.latitude - test_coord.latitude);
            double lon_diff = fabs(back_coord.longitude - test_coord.longitude);
            printf("  回算误差: Δlat=%.8f°, Δlon=%.8f°\n", lat_diff, lon_diff);
        }
    }
    else
    {
        printf("  WGS84 -> NAD83转换失败: %s\n", coord_get_error_string(ret));
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
        printf("  WGS84 -> NAD27转换失败: %s\n", coord_get_error_string(ret));
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
        printf("  WGS84 -> ED50转换失败: %s\n", coord_get_error_string(ret));
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
        printf("  WGS84 -> Tokyo转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试UTM转换
    printf("\nUTM转换测试:\n");
    UTMPoint utm;
    ret = coord_to_utm(ctx, &test_coord, &utm);
    if (ret == COORD_SUCCESS)
    {
        printf("  地理->UTM: %d%c %.3fE %.3fN\n",
               utm.zone, utm.band, utm.easting, utm.northing);
        // 转换回来
        GeoCoord geo_back;
        ret = coord_from_utm(ctx, &utm, &geo_back);
        if (ret == COORD_SUCCESS)
        {
            double lat_diff = fabs(geo_back.latitude - test_coord.latitude);
            double lon_diff = fabs(geo_back.longitude - test_coord.longitude);
            printf("  UTM->地理: (%.6f, %.6f), 误差: Δlat=%.8f°, Δlon=%.8f°\n",
                   geo_back.latitude, geo_back.longitude, lat_diff, lon_diff);
        }
    }
    else
    {
        printf("  UTM转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试MGRS转换
    printf("\nMGRS转换测试:\n");
    MGRSPoint mgrs;
    ret = coord_to_mgrs(ctx, &test_coord, &mgrs);
    if (ret == COORD_SUCCESS)
    {
        printf("  地理->MGRS: %d%c %s %05.0f %05.0f\n",
               mgrs.zone, mgrs.band, mgrs.square, mgrs.easting, mgrs.northing);
        // 转换回来
        GeoCoord geo_back_mgrs;
        ret = coord_from_mgrs(ctx, &mgrs, &geo_back_mgrs);
        if (ret == COORD_SUCCESS)
        {
            double lat_diff = fabs(geo_back_mgrs.latitude - test_coord.latitude);
            double lon_diff = fabs(geo_back_mgrs.longitude - test_coord.longitude);
            printf("  MGRS->地理: (%.6f, %.6f), 误差: Δlat=%.8f°, Δlon=%.8f°\n",
                   geo_back_mgrs.latitude, geo_back_mgrs.longitude, lat_diff, lon_diff);
        }
    }
    else
    {
        printf("  MGRS转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试British Grid转换
    printf("\nBritish Grid转换测试:\n");
    BritishGridPoint bg;
    ret = coord_to_british_grid(ctx, &test_coord, &bg);
    if (ret == COORD_SUCCESS)
    {
        printf("  地理->British Grid: %s %.0f %.0f\n",
               bg.letters, bg.easting, bg.northing);
    }
    else
    {
        printf("  British Grid转换失败: %s\n", coord_get_error_string(ret));
    }
    // 测试Japan Grid转换
    printf("\nJapan Grid转换测试:\n");
    JapanGridPoint jg;
    ret = coord_to_japan_grid(ctx, &test_coord, &jg);
    if (ret == COORD_SUCCESS)
    {
        printf("  地理->Japan Grid: Zone %d: %.3f, %.3f\n",
               jg.zone, jg.x, jg.y);
    }
    else
    {
        printf("  Japan Grid转换失败: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// 测试测地线计算
void test_geodesic_calculation()
{
    printf("=== 测试测地线计算 ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("创建上下文失败，无法进行测地线计算测试\n");
        return;
    }
    // 定义两个测试点：上海和北京
    GeoCoord shanghai = {31.230416, 121.473701, 0.0, DATUM_WGS84};
    GeoCoord beijing = {39.904211, 116.407394, 0.0, DATUM_WGS84};
    // 测试距离计算
    double distance, azi1, azi2;
    int ret = coord_distance(ctx, &shanghai, &beijing, &distance, &azi1, &azi2);
    if (ret == COORD_SUCCESS)
    {
        printf("上海到北京:\n");
        printf("  距离: %.2f 米 (约 %.2f 公里)\n", distance, distance / 1000.0);
        printf("  正向方位角: %.2f°\n", azi1);
        printf("  反向方位角: %.2f°\n", azi2);
    }
    else
    {
        printf("距离计算失败: %s\n", coord_get_error_string(ret));
    }
    // 测试正向计算
    printf("\n正向计算测试:\n");
    double test_distance = 100000.0; // 100公里
    double test_azimuth = 45.0; // 45度方向
    GeoCoord end_point;
    ret = coord_direct(ctx, &shanghai, test_distance, test_azimuth, &end_point);
    if (ret == COORD_SUCCESS)
    {
        printf("  从上海出发，沿%.1f°方向走%.0f米:\n", test_azimuth,
               test_distance);
        printf("  到达: (%.6f, %.6f)\n", end_point.latitude, end_point.longitude);
    }
    else
    {
        printf("  正向计算失败: %s\n", coord_get_error_string(ret));
    }
    // 测试反向计算
    printf("\n反向计算测试:\n");
    GeodesicResult result;
    ret = coord_inverse(ctx, &shanghai, &beijing, &result);
    if (ret == COORD_SUCCESS)
    {
        printf("  上海到北京反向计算:\n");
        printf("  距离: %.2f 米\n", result.distance);
        printf("  正向方位角: %.2f°\n", result.azimuth1);
        printf("  反向方位角: %.2f°\n", result.azimuth2);
    }
    else
    {
        printf("  反向计算失败: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// 测试基准转换工具
void test_datum_tools()
{
    printf("=== 测试基准转换工具 ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("创建上下文失败\n");
        return;
    }
    // 测试获取椭球体参数
    const Ellipsoid *ellipsoid = coord_get_ellipsoid(DATUM_WGS84);
    if (ellipsoid)
    {
        printf("WGS84椭球体参数:\n");
        printf("  长半轴: %.3f 米\n", ellipsoid->a);
        printf("  扁率: 1/%.9f\n", 1.0 / ellipsoid->f);
        printf("  短半轴: %.3f 米\n", ellipsoid->b);
        printf("  名称: %s\n", ellipsoid->name);
    }
    else
    {
        printf("获取椭球体参数失败\n");
    }
    // 测试NAD83椭球体参数
    ellipsoid = coord_get_ellipsoid(DATUM_NAD83);
    if (ellipsoid)
    {
        printf("\nNAD83椭球体参数:\n");
        printf("  长半轴: %.3f 米\n", ellipsoid->a);
        printf("  扁率: 1/%.9f\n", 1.0 / ellipsoid->f);
        printf("  名称: %s\n", ellipsoid->name);
    }
    // 测试设置和获取转换参数
    DatumTransform transform = {100.0, 200.0, 300.0, 1.0, 2.0, 3.0, 10.0};
    int ret = coord_set_transform_params(ctx, DATUM_WGS84, DATUM_TOKYO, &transform);
    if (ret == COORD_SUCCESS)
    {
        printf("\n设置转换参数成功\n");
        DatumTransform get_transform;
        ret = coord_get_transform_params(ctx, DATUM_WGS84, DATUM_TOKYO, &get_transform);
        if (ret == COORD_SUCCESS)
        {
            if (compare_double(transform.dx, get_transform.dx, 0.001) &&
                    compare_double(transform.dy, get_transform.dy, 0.001) &&
                    compare_double(transform.dz, get_transform.dz, 0.001))
            {
                printf("获取转换参数成功\n");
            }
            else
            {
                printf("获取的转换参数不匹配\n");
            }
        }
        else
        {
            printf("获取转换参数失败: %s\n", coord_get_error_string(ret));
        }
    }
    else
    {
        printf("设置转换参数失败: %s\n", coord_get_error_string(ret));
    }
    // 测试自定义椭球体
    ret = coord_set_custom_ellipsoid(ctx, 6371000.0, 1.0 / 298.3);
    if (ret == COORD_SUCCESS)
    {
        printf("设置自定义椭球体成功\n");
    }
    else
    {
        printf("设置自定义椭球体失败: %s\n", coord_get_error_string(ret));
    }
    coord_destroy_context(ctx);
    printf("\n");
}

// 测试错误处理
void test_error_handling()
{
    printf("=== 测试错误处理 ===\n");
    // 测试无效输入
    CoordContext *ctx = NULL;
    GeoCoord invalid_coord = {100.0, 200.0, 0.0, DATUM_WGS84};
    char buffer[256];
    int ret = coord_convert(ctx, &invalid_coord, COORD_FORMAT_DD, DATUM_WGS84,
                            buffer, sizeof(buffer));
    printf("空上下文转换测试: %s (期望: 无效输入)\n",
           ret == COORD_ERROR_INVALID_INPUT ? "通过" : "失败");
    // 测试错误信息获取
    printf("\n错误信息测试:\n");
    for (int i = 0; i <= 10; i++)
    {
        const char *error_msg = coord_get_error_string(i);
        printf("  错误码 %d: %s\n", i, error_msg);
    }
    printf("\n");
}

// 综合测试
void test_comprehensive()
{
    printf("=== 综合测试 ===\n");
    CoordContext *ctx = coord_create_context(DATUM_WGS84);
    if (!ctx)
    {
        printf("创建上下文失败\n");
        return;
    }
    // 定义多个测试点
    struct
    {
        const char *name;
        double lat, lon;
    } test_points[] =
    {
        {"上海", 31.230416, 121.473701},
        {"北京", 39.904211, 116.407394},
        {"纽约", 40.712776, -74.005974},
        {"伦敦", 51.507351, -0.127758},
        {"悉尼", -33.868820, 151.209290},
        {"东京", 35.689487, 139.691711},
        {"巴黎", 48.856614, 2.352222}
    };
    int num_points = sizeof(test_points) / sizeof(test_points[0]);
    // 测试每个点的格式转换
    for (int i = 0; i < num_points; i++)
    {
        printf("%s坐标转换测试:\n", test_points[i].name);
        GeoCoord coord = {test_points[i].lat, test_points[i].lon, 0.0, DATUM_WGS84};
        char buffer[256];
        // 测试各种格式
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
                printf("  %s格式失败: %s\n", format_names[j], coord_get_error_string(ret));
            }
        }
        // 测试UTM转换
        UTMPoint utm;
        int ret = coord_to_utm(ctx, &coord, &utm);
        if (ret == COORD_SUCCESS)
        {
            printf("  UTM: 区域 %d%c\n", utm.zone, utm.band);
        }
        else
        {
            printf("  UTM转换失败: %s\n", coord_get_error_string(ret));
        }
        // 测试MGRS转换
        MGRSPoint mgrs;
        ret = coord_to_mgrs(ctx, &coord, &mgrs);
        if (ret == COORD_SUCCESS)
        {
            printf("  MGRS: 区域 %d%c\n", mgrs.zone, mgrs.band);
        }
        else
        {
            printf("  MGRS转换失败: %s\n", coord_get_error_string(ret));
        }
        printf("\n");
    }
    // 测试点对点距离计算
    printf("点对点距离计算:\n");
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
                printf("  %s -> %s: %.2f 公里\n",
                       test_points[i].name, test_points[j].name, distance / 1000.0);
            }
        }
    }
    printf("\nMGRS坐标转换测试:\n");
    // 测试MGRS坐标转换
    struct
    {
        const char *name;
        double lat, lon;
        const char *expected_mgrs;
    } mgrs_test_points[] =
    {
        {"上海", 31.230416, 121.473701, "51R"},
        {"北京", 39.904211, 116.407394, "50S"},
        {"悉尼", -33.868820, 151.209290, "56H"}
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
                printf("  %s: %s (期望区域: %s)\n",
                       mgrs_test_points[i].name, buffer, mgrs_test_points[i].expected_mgrs);
                // 转换回来
                GeoCoord back_coord;
                ret = coord_from_mgrs(ctx, &mgrs, &back_coord);
                if (ret == COORD_SUCCESS)
                {
                    double lat_diff = fabs(back_coord.latitude - mgrs_test_points[i].lat);
                    double lon_diff = fabs(back_coord.longitude - mgrs_test_points[i].lon);
                    printf("    回转误差: Δlat=%.6f°, Δlon=%.6f°\n", lat_diff, lon_diff);
                }
            }
        }
        else
        {
            printf("  %s MGRS转换失败: %s\n", mgrs_test_points[i].name,
                   coord_get_error_string(ret));
        }
    }
    coord_destroy_context(ctx);
    printf("\n");
}

int main()
{
    printf("=== 坐标转换系统增强测试程序 ===\n\n");
    // 设置错误回调
    coord_set_error_callback(error_handler);
    // 运行所有测试
    test_context_creation();
    test_utility_functions();
    test_coord_parsing();
    test_coord_formatting();
    test_coord_conversion();
    test_geodesic_calculation();
    test_datum_tools();
    test_error_handling();
    test_comprehensive();
    printf("=== 所有测试完成 ===\n");
    return 0;
}