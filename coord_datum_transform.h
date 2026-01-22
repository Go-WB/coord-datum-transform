#ifndef COORD_TRANSFORM_H
#define COORD_TRANSFORM_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

// 前向声明GeographicLib结构体
struct geod_geodesic;

// 坐标格式枚举
typedef enum
{
    COORD_FORMAT_DD = 0,        // 十进制度 DD.ddddd°
    COORD_FORMAT_DMM,           // 度分 DD°MM.mmm'
    COORD_FORMAT_DMS,           // 度分秒 DD°MM'SS"
    COORD_FORMAT_UTM,           // UTM坐标
    COORD_FORMAT_MGRS,          // MGRS坐标（默认）
    COORD_FORMAT_BRITISH_GRID,  // 英国网格
    COORD_FORMAT_JAPAN_GRID,    // 日本网格
    COORD_FORMAT_MAX
} CoordFormat;

// 地图基准枚举
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

// 椭球体参数
typedef struct
{
    double a;           // 长半轴 (meters)
    double f;           // 扁率
    double b;           // 短半轴 (meters)
    double e2;          // 第一偏心率平方
    double ep2;         // 第二偏心率平方
    const char *name;   // 椭球体名称
} Ellipsoid;

// 基准转换七参数
typedef struct
{
    double dx, dy, dz;          // 平移参数 (meters)
    double rx, ry, rz;          // 旋转参数 (arc-seconds)
    double scale;               // 尺度因子 (ppm)
} DatumTransform;

// 地理坐标
typedef struct
{
    double latitude;            // 纬度 (degrees)
    double longitude;           // 经度 (degrees)
    double altitude;            // 海拔高度 (meters)
    MapDatum datum;             // 坐标基准
} GeoCoord;

// UTM坐标
typedef struct
{
    int zone;                   // UTM区域 (1-60)
    char band;                  // 纬度带 (C-X)
    double easting;             // 东距 (meters)
    double northing;            // 北距 (meters)
    double convergence;         // 子午线收敛角 (degrees)
    double scale_factor;        // 比例因子
    MapDatum datum;             // 基准
} UTMPoint;

// MGRS坐标
typedef struct
{
    int zone;                   // UTM区域 (1-60)
    char band;                  // 纬度带 (C-X)
    char square[3];             // 100km网格方 (2字符)
    double easting;             // 东距 (meters, 在网格方内)
    double northing;            // 北距 (meters, 在网格方内)
    MapDatum datum;             // 基准
} MGRSPoint;

// 英国国家网格坐标
typedef struct
{
    char letters[3];           // 两字母编码
    double easting;             // 东距 (meters)
    double northing;            // 北距 (meters)
    MapDatum datum;             // 基准
} BritishGridPoint;

// 日本网格坐标
typedef struct
{
    int zone;                   // 区域号
    double x;                   // X坐标
    double y;                   // Y坐标
    MapDatum datum;             // 基准
} JapanGridPoint;

// 解析结果
typedef struct
{
    int success;                // 是否成功
    GeoCoord coord;             // 解析出的坐标
    CoordFormat format;         // 检测到的格式
    MapDatum datum;             // 检测到的基准
    char error_msg[256];        // 错误信息
} ParseResult;

// 测地线结果
typedef struct
{
    double distance;            // 距离 (meters)
    double azimuth1;            // 正向方位角 (degrees)
    double azimuth2;            // 反向方位角 (degrees)
} GeodesicResult;

// 坐标转换上下文
typedef struct
{
    struct geod_geodesic *geod;  // GeographicLib测地线对象指针
    Ellipsoid ellipsoid;        // 当前椭球体
    DatumTransform transforms[DATUM_MAX][DATUM_MAX]; // 转换参数表
} CoordContext;

// ============================ 对外接口列表 ============================

// 错误码
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

// ==================== 初始化与清理函数 ====================
CoordContext *coord_create_context(MapDatum datum);
void coord_destroy_context(CoordContext *ctx);
int coord_set_datum(CoordContext *ctx, MapDatum datum);

// ==================== 坐标解析函数 ====================
ParseResult coord_parse_string(const char *str, CoordFormat format,
                               MapDatum datum);
ParseResult coord_auto_parse(const char *str);

// ==================== 坐标格式化函数 ====================
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

// ==================== 坐标转换函数 ====================
// 地理坐标转其他格式
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

// 基准转换
int coord_convert_datum(CoordContext *ctx, const GeoCoord *src,
                        MapDatum target_datum, GeoCoord *dst);

// ==================== 测地线计算 ====================
int coord_distance(CoordContext *ctx, const GeoCoord *p1, const GeoCoord *p2,
                   double *distance, double *azi1, double *azi2);
int coord_direct(CoordContext *ctx, const GeoCoord *start,
                 double distance, double azimuth, GeoCoord *end);
int coord_inverse(CoordContext *ctx, const GeoCoord *p1, const GeoCoord *p2,
                  GeodesicResult *result);

// ==================== 工具函数 ====================
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

// ==================== 基准转换工具 ====================
int coord_set_transform_params(CoordContext *ctx, MapDatum from, MapDatum to,
                               const DatumTransform *params);
int coord_get_transform_params(CoordContext *ctx, MapDatum from, MapDatum to,
                               DatumTransform *params);

// ==================== 椭球体工具 ====================
const Ellipsoid *coord_get_ellipsoid(MapDatum datum);
int coord_set_custom_ellipsoid(CoordContext *ctx, double a, double f);

// ==================== 错误处理 ====================
const char *coord_get_error_string(int error_code);
void coord_set_error_callback(void (*callback)(int, const char *));

// ==================== 格式转换主函数 ====================
int coord_convert(CoordContext *ctx, const GeoCoord *src,
                  CoordFormat target_format, MapDatum target_datum,
                  char *result_buffer, size_t buffer_size);

#endif // COORD_TRANSFORM_H