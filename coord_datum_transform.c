#include "coord_transform.h"
#include "geodesic.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

// 常量定义
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 / M_PI)
#define ARC_SEC_TO_RAD (M_PI / (180.0 * 3600.0))
#define PPM_TO_SCALE 1e-6
#define METERS_TO_FEET 3.280839895
#define FEET_TO_METERS 0.3048

// 椭球体定义
static const Ellipsoid ELLIPSOIDS[] =
{
    // WGS84
    {
        6378137.0, 1.0 / 298.257223563, 6356752.314245,
        0.0066943799901413165, 0.0067394967422764341, "WGS84"
    },
    // MGRS Grid (使用WGS84椭球)
    {
        6378137.0, 1.0 / 298.257223563, 6356752.314245,
        0.0066943799901413165, 0.0067394967422764341, "WGS84"
    },
    // UTM Grid (使用WGS84椭球)
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
    // Airy 1830 (OSGB36 - 英国国家网格)
    {
        6377563.396, 1.0 / 299.3249646, 6356256.909,
        0.0066705397616, 0.006715826523, "Airy1830"
    }
};

// 英国国家网格参数
static const double OSGB36_A = 6377563.396;    // Airy 1830椭球长半轴
static const double OSGB36_F = 1.0 / 299.3249646; // Airy 1830扁率
static const double OSGB36_N0 = -100000.0;     // 北偏移
static const double OSGB36_E0 = 400000.0;     // 东偏移
static const double OSGB36_F0 = 0.9996012717; // 中央子午线比例因子
static const double OSGB36_LAT0 = 49.0 * DEG_TO_RAD; // 真原点纬度
static const double OSGB36_LON0 = -2.0 * DEG_TO_RAD; // 真原点经度

// 日本网格参数 (Tokyo Datum, Bessel 1841椭球)
static const double JAPAN_GRID_A = 6377397.155;
static const double JAPAN_GRID_F = 1.0 / 299.1528128;

// 错误信息
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

// 全局错误回调
static void (*error_callback)(int, const char *) = NULL;

// 设置错误
static void set_error(int code, const char *message)
{
    if (error_callback)
    {
        error_callback(code, message);
    }
}

// ==================== 基础工具函数 ====================
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

// ==================== 上下文管理 ====================
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
    // 设置椭球体
    ctx->ellipsoid = ELLIPSOIDS[datum];
    // 初始化GeographicLib测地线对象
    ctx->geod = (struct geod_geodesic *)malloc(sizeof(struct geod_geodesic));
    if (!ctx->geod)
    {
        free(ctx);
        set_error(COORD_ERROR_MEMORY, "Failed to create geodesic object");
        return NULL;
    }
    geod_init(ctx->geod, ctx->ellipsoid.a, ctx->ellipsoid.f);
    // 初始化转换参数表
    memset(ctx->transforms, 0, sizeof(ctx->transforms));
    // 设置默认转换参数
    // WGS84 <-> NAD83 (基本相同)
    ctx->transforms[DATUM_WGS84][DATUM_NAD83].dx = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD83].dy = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD83].dz = 0.0;
    ctx->transforms[DATUM_NAD83][DATUM_WGS84] =
        ctx->transforms[DATUM_WGS84][DATUM_NAD83];
    // WGS84 <-> MGRS Grid (相同)
    ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID].dx = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID].dy = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID].dz = 0.0;
    ctx->transforms[DATUM_MGRS_GRID][DATUM_WGS84] =
        ctx->transforms[DATUM_WGS84][DATUM_MGRS_GRID];
    // WGS84 <-> UTM Grid (相同)
    ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID].dx = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID].dy = 0.0;
    ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID].dz = 0.0;
    ctx->transforms[DATUM_UTM_GRID][DATUM_WGS84] =
        ctx->transforms[DATUM_WGS84][DATUM_UTM_GRID];
    // WGS84 -> NAD27 (NADCON 参数, CONUS)
    // 来源: National Geodetic Survey
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].dx = -8.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].dy = 160.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].dz = 176.0;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].rx = -0.25;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].ry = 0.75;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].rz = -0.06;
    ctx->transforms[DATUM_WGS84][DATUM_NAD27].scale = -0.34;
    // WGS84 -> ED50 (EPSG 参数)
    // 来源: EPSG Dataset
    ctx->transforms[DATUM_WGS84][DATUM_ED50].dx = -87.0;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].dy = -98.0;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].dz = -121.0;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].rx = -0.59;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].ry = -0.32;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].rz = -1.12;
    ctx->transforms[DATUM_WGS84][DATUM_ED50].scale = -3.72;
    // WGS84 -> Tokyo (近似参数)
    ctx->transforms[DATUM_WGS84][DATUM_TOKYO].dx = -148.0;
    ctx->transforms[DATUM_WGS84][DATUM_TOKYO].dy = 507.0;
    ctx->transforms[DATUM_WGS84][DATUM_TOKYO].dz = 685.0;
    // WGS84 -> OSGB36 (OSTN15 参数)
    // 来源: Ordnance Survey National Grid (OSTN15)
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

// ==================== UTM区域计算 ====================
int coord_get_utm_zone(double longitude, double latitude)
{
    if (!coord_is_valid_longitude(longitude) || !coord_is_valid_latitude(latitude))
    {
        return 0;
    }
    // 标准化经度
    double lon_norm = longitude;
    while (lon_norm < -180.0)
    {
        lon_norm += 360.0;
    }
    while (lon_norm >= 180.0)
    {
        lon_norm -= 360.0;
    }
    // 特殊区域处理
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
    // 标准UTM区域
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
    // UTM纬度带定义表（8度一个带，跳过I和O）
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
    return 'Z'; // 不应该到达这里
}

// ==================== 坐标验证 ====================
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
    // 放宽东距检查范围
    if (utm->easting < 100000.0 || utm->easting > 900000.0)
    {
        return 0;
    }
    // 放宽北距检查范围，考虑南北半球
    if (utm->band >= 'N' && utm->band <= 'X')
    {
        // 北半球：0-10,000,000米
        if (utm->northing < 0.0 || utm->northing > 10000000.0)
        {
            return 0;
        }
    }
    else
    {
        // 南半球：10,000,000-20,000,000米（假北距）
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
    // 验证UTM区域
    if (mgrs->zone < 1 || mgrs->zone > 60)
    {
        return 0;
    }
    // 验证纬度带
    if (mgrs->band < 'C' || mgrs->band > 'X' ||
            mgrs->band == 'I' || mgrs->band == 'O')
    {
        return 0;
    }
    // 验证网格方标识符
    if (strlen(mgrs->square) != 2)
    {
        return 0;
    }
    // 验证网格方字母（跳过I和O）
    if (mgrs->square[0] < 'A' || mgrs->square[0] > 'Z' ||
            mgrs->square[0] == 'I' || mgrs->square[0] == 'O' ||
            mgrs->square[1] < 'A' || mgrs->square[1] > 'Z' ||
            mgrs->square[1] == 'I' || mgrs->square[1] == 'O')
    {
        return 0;
    }
    // 验证东距和北距
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

// ==================== 坐标解析 ====================
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
    // 跳过前导空格
    while (isspace((unsigned char)*str))
    {
        str++;
    }
    switch (format)
    {
        case COORD_FORMAT_DD:
        {
            // 格式: "31.230416°N, 121.473701°E" 或 "31.230416, 121.473701"
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
                // 尝试不带方向字符的格式
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
            // 格式: "31°13'49.5\"N, 121°28'25.32\"E"
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
            // 格式: "31°13.825'N, 121°28.422'E"
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
            // 格式: "50N 447600E 4419300N" 或 "50N 447600 4419300"
            int zone;
            char band;
            double easting, northing;
            char east_dir = 'E', north_dir = 'N';
            int count = sscanf(str, "%d%c %lf%c %lf%c",
                               &zone, &band, &easting, &east_dir, &northing, &north_dir);
            if (count != 6)
            {
                // 尝试不带方向字符的格式
                count = sscanf(str, "%d%c %lf %lf", &zone, &band, &easting, &northing);
                if (count != 4)
                {
                    strcpy(result.error_msg, "Failed to parse UTM format");
                    return result;
                }
            }
            // 创建UTM点
            UTMPoint utm = {zone, band, easting, northing, 0.0, 0.9996, datum};
            if (!coord_validate_utm(&utm))
            {
                strcpy(result.error_msg, "Invalid UTM coordinate");
                return result;
            }
            // 转换为地理坐标
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
            // 格式: "51Q SB 54634 56142" 或 "51QSB 54634 56142"
            int zone;
            char band;
            char square[3] = {0};  // 网格方 (2个字母 + null终止符)
            double easting, northing;
            // 尝试多种格式
            int count = sscanf(str, "%d%c%2s %lf %lf",
                               &zone, &band, square, &easting, &northing);
            if (count != 5)
            {
                // 尝试带空格的格式: "51Q SB 54634 56142"
                count = sscanf(str, "%d%c %2s %lf %lf",
                               &zone, &band, square, &easting, &northing);
                if (count != 5)
                {
                    strcpy(result.error_msg, "Failed to parse MGRS format");
                    return result;
                }
            }
            // 验证MGRS参数
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
            // 验证网格方字母
            if (square[0] < 'A' || square[0] > 'Z' || square[0] == 'I' || square[0] == 'O'
                    ||
                    square[1] < 'A' || square[1] > 'Z' || square[1] == 'I' || square[1] == 'O')
            {
                strcpy(result.error_msg, "Invalid MGRS square letters");
                return result;
            }
            // 验证东距和北距
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
            // 创建MGRS点
            MGRSPoint mgrs;
            mgrs.zone = zone;
            mgrs.band = band;
            mgrs.square[0] = square[0];
            mgrs.square[1] = square[1];
            mgrs.square[2] = '\0';
            mgrs.easting = easting;
            mgrs.northing = northing;
            mgrs.datum = datum;
            // 使用coord_validate_mgrs验证
            if (!coord_validate_mgrs(&mgrs))
            {
                strcpy(result.error_msg, "Invalid MGRS coordinate");
                return result;
            }
            // 转换为地理坐标
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
            // 格式: "TQ 12345 67890" 或 "TQ1234567890"
            char letters[3] = {0};
            double easting, northing;
            int count = sscanf(str, "%2s %lf %lf", letters, &easting, &northing);
            if (count != 3)
            {
                // 尝试不带空格的格式: "TQ1234567890"
                char buffer[32];
                strncpy(buffer, str, sizeof(buffer) - 1);
                buffer[sizeof(buffer) - 1] = '\0';
                if (strlen(buffer) >= 2)
                {
                    letters[0] = buffer[0];
                    letters[1] = buffer[1];
                    letters[2] = '\0';
                    // 解析数字部分
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
            // 创建British Grid点
            BritishGridPoint bg;
            strncpy(bg.letters, letters, 3);
            bg.easting = easting;
            bg.northing = northing;
            bg.datum = datum;
            // 转换为地理坐标
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
            // 格式: "Zone 3: 12345.6, 67890.1" 或 "3 12345.6 67890.1"
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
            // 创建Japan Grid点
            JapanGridPoint jg;
            jg.zone = zone;
            jg.x = x;
            jg.y = y;
            jg.datum = datum;
            // 转换为地理坐标
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

// 在coord_auto_parse函数中，添加对MGRS格式的自动检测
ParseResult coord_auto_parse(const char *str)
{
    ParseResult result = {0};
    if (!str)
    {
        strcpy(result.error_msg, "Input string is NULL");
        return result;
    }
    // 跳过前导空格
    const char *s = str;
    while (isspace((unsigned char)*s))
    {
        s++;
    }
    // 检查是否是MGRS格式
    int zone;
    char band;
    char square[3];
    double easting, northing;
    // 尝试解析MGRS格式
    int count = sscanf(s, "%d%c%2s %lf %lf", &zone, &band, square, &easting,
                       &northing);
    if (count == 5)
    {
        // 验证MGRS参数
        if (zone >= 1 && zone <= 60 &&
                band >= 'C' && band <= 'X' && band != 'I' && band != 'O' &&
                strlen(square) == 2 &&
                square[0] >= 'A' && square[0] <= 'Z' && square[0] != 'I' && square[0] != 'O' &&
                square[1] >= 'A' && square[1] <= 'Z' && square[1] != 'I' && square[1] != 'O' &&
                easting >= 0.0 && easting <= 100000.0 &&
                northing >= 0.0 && northing <= 100000.0)
        {
            // 看起来像MGRS格式
            result = coord_parse_string(str, COORD_FORMAT_MGRS, DATUM_WGS84);
            if (result.success)
            {
                return result;
            }
        }
    }
    // 检查是否是UTM格式
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
    // 检查是否是British Grid格式
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
    // 检查是否是Japan Grid格式
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
    // 尝试其他格式
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


// ==================== 坐标格式化函数 ====================
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

// ==================== 坐标转换函数 ====================
// 地理坐标转UTM
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
    // 计算UTM区域
    int zone = coord_get_utm_zone(geo->longitude, geo->latitude);
    if (zone < 1 || zone > 60)
    {
        return COORD_ERROR_INVALID_UTM_ZONE;
    }
    // 计算中央子午线
    double lon_center = (zone - 1) * 6.0 - 180.0 + 3.0;
    // 转换为弧度
    double lat_rad = coord_deg_to_rad(geo->latitude);
    double lon_rad = coord_deg_to_rad(geo->longitude);
    double lon_center_rad = coord_deg_to_rad(lon_center);
    // UTM转换参数
    double k0 = 0.9996;  // UTM比例因子
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
    // 计算M（子午线弧长）
    double M = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 /
                     256.0) * lat_rad
                    - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(
                        2.0 * lat_rad)
                    + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(4.0 * lat_rad)
                    - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * lat_rad));
    // 计算UTM坐标
    double A2 = A * A;
    double A3 = A2 * A;
    double A4 = A3 * A;
    double A5 = A4 * A;
    double A6 = A5 * A;
    // 东距
    utm->easting = k0 * N * (A + (1.0 - T + C) * A3 / 6.0
                             + (5.0 - 18.0 * T + T * T + 72.0 * C - 58.0 * e2) * A5 / 120.0)
                   + 500000.0;  // 假东距
    // 北距
    utm->northing = k0 * (M + N * tan_lat *
                          (A2 / 2.0 + (5.0 - T + 9.0 * C + 4.0 * C * C) * A4 / 24.0
                           + (61.0 - 58.0 * T + T * T + 600.0 * C - 330.0 * e2) * A6 / 720.0));
    // 如果是南半球，加上假北距
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
    // 计算中央子午线
    double lon_center = (utm->zone - 1) * 6.0 - 180.0 + 3.0;
    double k0 = 0.9996;
    double a = ctx->ellipsoid.a;
    double f = ctx->ellipsoid.f;
    double e2 = 2 * f - f * f;
    // 移除假东距
    double x = utm->easting - 500000.0;
    double y = utm->northing;
    // 如果是南半球，移除假北距
    if (utm->band < 'N')
    {
        y -= 10000000.0;
    }
    // 计算脚点纬度
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

// 辅助函数：获取字母在MGRS字母表中的索引（跳过I和O）
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

// 辅助函数：从索引获取MGRS字母（跳过I和O）
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

// ==================== 修复MGRS转换的关键函数 ====================
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

    // MGRS 100k 网格列字母反向计算 (WGS84/现代基准)
    // 使用与 coord_to_mgrs 相同的 6-set 系统进行反向计算
    //
    // 从列字母反推 col_100k 的逻辑:
    // 1. 确定 set (zone % 6)
    // 2. 获取该 set 的起始字母
    // 3. 从起始字母开始计数，直到达到目标列字母

    // 确定 set (1-6)
    int set = zone % 6;
    if (set == 0)
        set = 6;

    // 获取该 set 的起始列字母
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

    // 从起始字母计数到目标列字母，跳过 I 和 O
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

    // 行字母反向计算
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

    // 处理南半球北距
    // 南半球的UTM北距 = 10000000 + 真实北距（从赤道向南为正）
    // 但 MGRS 需要使用真实北距来计算 100km 网格
    double utm_northing_for_mgrs = utm.northing;
    if (band < 'N')
    {
        // 南半球：移除假北距，得到真实北距（相对于南半球起点）
        // 注意：在南半球，真实北距可能是负数（相对于赤道）
        // 但 MGRS 在南半球使用相对于南极方向计算
        utm_northing_for_mgrs -= 10000000.0;
    }

    // 计算 100km 网格索引
    int row_100k = (int)(utm_northing_for_mgrs / 100000.0);

    // 南半球特殊处理：row_100k 可能是负数
    // MGRS 在南半球的 row_100k 从 0 开始递增（向南）
    // 如果是负数，说明坐标在南半球高纬度区域
    if (row_100k < 0)
    {
        // 南半球极地区域，需要特殊处理
        // 将负数索引转换为正数索引
        row_100k = 20 + (row_100k % 20);
    }

    // MGRS 100k 网格列字母计算 (WGS84/现代基准)
    // 参考: GeoTrellis MGRS 实现 (基于 proj4js/mgrs)
    // 使用 6-set 循环系统，而不是简单的 (zone-1)%3
    //
    // SET_ORIGIN_COLUMN_LETTERS = "AJSAJS" (6 个 set 循环)
    // - Zone % 6 = 1: A (ASCII 65)
    // - Zone % 6 = 2: J (ASCII 74)
    // - Zone % 6 = 3: S (ASCII 83)
    // - Zone % 6 = 4: A (ASCII 65)
    // - Zone % 6 = 5: J (ASCII 74)
    // - Zone % 6 = 0: S (ASCII 65)
    //
    // 计算公式: col_letter = origin + col_100k - 1
    // 需要跳过 I (73) 和 O (79)
    //
    // 示例: Zone 50, col_100k=5
    // - set = 50 % 6 = 2
    // - origin = 'J' (ASCII 74)
    // - col = 'J' + 5 - 1 = 78 = 'N' ✓

    // 确定 set (1-6)
    int set = zone % 6;
    if (set == 0)
        set = 6;

    // 获取该 set 的起始列字母
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

    // 计算列字母，跳过 I 和 O
    int col_letter_code = col_origin + col_100k - 1;
    char col_letter;

    // 处理 I 和 O 的跳过逻辑
    int temp_code = col_origin;
    for (int i = 1; i < col_100k; i++)
    {
        temp_code++;
        if (temp_code == 'I' || temp_code == 'O')
            temp_code++;
    }
    col_letter_code = temp_code;

    // 处理循环 (超过 Z 时回到 A)
    if (col_letter_code > 'Z')
    {
        col_letter_code = col_letter_code - 'Z' + 'A' - 1;
        // 再次检查 I 和 O
        if (col_letter_code == 'I' || col_letter_code == 'O')
            col_letter_code++;
    }

    col_letter = (char)col_letter_code;

    // 计算行字母 (使用原有的奇偶 zone 逻辑)
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
    // 确保 northing 为正数
    if (mgrs->northing < 0)
    {
        mgrs->northing += 100000.0;
    }
    mgrs->datum = utm.datum;
    return COORD_SUCCESS;
}

// 地理坐标转英国网格
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

    // British National Grid 必须使用 OSGB36 基准和 Airy 1830 椭球
    // 如果输入不是 OSGB36，需要先进行基准转换
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

    // 使用 OSGB36/Airy 1830 椭球参数
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
    // 计算M
    double M = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 /
                     256.0) * lat_rad
                    - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(
                        2.0 * lat_rad)
                    + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(4.0 * lat_rad)
                    - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * lat_rad));
    // 计算M0
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

    // 计算英国国家网格字母
    // 英国网格使用特殊的 500km 方字母系统
    // 注意：对于超出英国范围的坐标，字母没有标准定义
    // 这里使用扩展的循环计算方式

    // 计算 500km 方索引
    int e500k = (int)(bg->easting / 500000.0);
    int n500k = (int)(bg->northing / 500000.0);

    // 处理负值索引
    while (e500k < 0) e500k += 25;  // 确保正数
    while (n500k < 0) n500k += 25;

    // 100km 方内的字母
    int e100k = (int)(fmod(fabs(bg->easting), 500000.0) / 100000.0);
    int n100k = (int)(fmod(fabs(bg->northing), 500000.0) / 100000.0);

    // 使用英国网格字母表（跳过 I）
    const char *bg_letters = "ABCDEFGHJKLMNPQRSTUVWXYZ";

    // 计算东向字母（循环使用字母表）
    int e_idx = (e500k * 5 + e100k) % 25;
    bg->letters[0] = bg_letters[e_idx];

    // 计算北向字母（循环使用字母表）
    int n_idx = (n500k * 5 + n100k) % 25;
    bg->letters[1] = bg_letters[n_idx];

    bg->letters[2] = '\0';
    bg->datum = DATUM_OSGB36;  // British Grid 始终是 OSGB36 基准
    return COORD_SUCCESS;
}

int coord_from_british_grid(CoordContext *ctx, const BritishGridPoint *bg,
                            GeoCoord *geo)
{
    if (!ctx || !bg || !geo)
    {
        return COORD_ERROR_INVALID_INPUT;
    }

    // 英国国家网格 (OSGB36) 使用 Airy 1830 椭球
    // 首先将英国网格坐标转换为 OSGB36 经纬度
    // 然后转换为 WGS84

    // 移除东偏移和北偏移
    double E = bg->easting;
    double N = bg->northing;

    // 临时使用 WGS84 椭球进行简化计算
    // 注意：完整的 OSGB36 转换需要使用 Airy 1830 椭球参数
    double a = OSGB36_A;
    double f = OSGB36_F;
    double e2 = 2 * f - f * f;

    // 计算中央经线和纬度原点
    double lon0 = OSGB36_LON0;
    double E0 = OSGB36_E0;
    double N0 = OSGB36_N0;
    double F0 = OSGB36_F0;

    // 初始估计纬度
    double M_prime = (N - N0) / F0;
    double mu_prime = M_prime / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0
                                   - 5.0 * e2 * e2 * e2 / 256.0));

    // 计算初始纬度
    double e1 = (1.0 - sqrt(1.0 - e2)) / (1.0 + sqrt(1.0 - e2));
    double J1 = 3.0 * e1 / 2.0 - 27.0 * e1 * e1 * e1 / 32.0;
    double J2 = 21.0 * e1 * e1 / 16.0 - 55.0 * e1 * e1 * e1 * e1 / 32.0;
    double J3 = 151.0 * e1 * e1 * e1 / 96.0;
    double J4 = 1097.0 * e1 * e1 * e1 * e1 / 512.0;

    double lat_prime = mu_prime + J1 * sin(2.0 * mu_prime) + J2 * sin(4.0 * mu_prime)
                      + J3 * sin(6.0 * mu_prime) + J4 * sin(8.0 * mu_prime);

    // 迭代计算精确纬度和经度
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

        // 检查收敛
        if (fabs(lat_delta) < 1e-12)
        {
            break;
        }
    }

    // 转换为度
    geo->latitude = coord_rad_to_deg(lat_rad);
    geo->longitude = coord_rad_to_deg(lon_rad);
    geo->altitude = 0.0;

    // 从 OSGB36 转换到 WGS84（简化近似）
    // 英国区域的近似 Helmert 转换参数
    double tx = -446.448;  // 米
    double ty = 125.157;   // 米
    double tz = -542.060;  // 米
    double rx = -0.1502;   // 弧秒
    double ry = -0.2470;   // 弧秒
    double rz = -0.8421;   // 弧秒
    double s = 20.4894;    // ppm

    // 使用七参数转换
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

    // 应用七参数转换
    double rx_rad = rx * ARC_SEC_TO_RAD;
    double ry_rad = ry * ARC_SEC_TO_RAD;
    double rz_rad = rz * ARC_SEC_TO_RAD;
    double scale_factor = 1.0 + s * PPM_TO_SCALE;

    double X2 = tx + X * scale_factor + Y * rz_rad - Z * ry_rad;
    double Y2 = ty - X * rz_rad + Y * scale_factor + Z * rx_rad;
    double Z2 = tz + X * ry_rad - Y * rx_rad + Z * scale_factor;

    // 转换回 WGS84 大地坐标
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

// 日本平面直角坐标系区域参数
static const struct
{
    int zone;
    double lat0;
    double lon0;
    double false_e;   // 东偏移 (假东距)
    double false_n;   // 北偏移 (假北距)
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

// 地理坐标转日本网格
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
    // 转换为Tokyo基准
    GeoCoord tokyo_geo;
    int ret = coord_convert_datum(ctx, geo, DATUM_TOKYO, &tokyo_geo);
    if (ret != COORD_SUCCESS)
    {
        return ret;
    }

    double lat = tokyo_geo.latitude;
    double lon = tokyo_geo.longitude;

    // 根据经纬度查找合适的区域（选择最近的中央经线）
    // 不限制地理范围，支持任意坐标转换
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

    // 获取区域参数
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

    // Bessel 1841 椭球参数
    double a = JAPAN_GRID_A;
    double f = JAPAN_GRID_F;
    double e2 = 2 * f - f * f;

    // 计算子午线弧长 M
    double M = a * ((1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 / 256.0) * lat_rad
                   - (3.0 * e2 / 8.0 + 3.0 * e2 * e2 / 32.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(2.0 * lat_rad)
                   + (15.0 * e2 * e2 / 256.0 + 45.0 * e2 * e2 * e2 / 1024.0) * sin(4.0 * lat_rad)
                   - (35.0 * e2 * e2 * e2 / 3072.0) * sin(6.0 * lat_rad));

    // 计算辅助参数
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double T = tan_lat * tan_lat;
    double C = e2 * cos_lat * cos_lat / (1.0 - e2);
    double A = (lon_rad - lon0) * cos_lat;
    double A2 = A * A;
    double A3 = A2 * A;
    double A4 = A3 * A;
    double A5 = A4 * A;
    double A6 = A5 * A;

    // 计算 X（北坐标）
    jg->x = k0 * (M + N * tan_lat * (A2 / 2.0 + (5.0 - T + 9.0 * C + 4.0 * C * C) * A4 / 24.0
               + (61.0 - 58.0 * T + T * T + 600.0 * C - 330.0 * e2) * A6 / 720.0)) + false_n;

    // 计算 Y（东坐标）
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

    // 获取 Bessel 1841 椭球参数
    double a = JAPAN_GRID_A;
    double f = JAPAN_GRID_F;
    double e2 = 2 * f - f * f;

    // 查找区域参数（使用外部定义的 japan_zones 数组）
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

    // 高斯-克吕格投影反向转换
    // 移除假偏移
    // 注意：jg->x 是北坐标，jg->y 是东坐标
    double northing = jg->x - false_n;   // 北坐标
    double easting = jg->y - false_e;     // 东坐标

    // 计算辅助参数
    // M 是子午线弧长，从北坐标计算
    double M = northing / k0;
    double mu = M / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0
                        - 5.0 * e2 * e2 * e2 / 256.0));

    // 计算脚点纬度
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

    // D 从东坐标计算
    double D = easting / (N1 * k0);
    double D2 = D * D;
    double D3 = D2 * D;
    double D4 = D3 * D;
    double D5 = D4 * D;
    double D6 = D5 * D;

    // 计算纬度修正
    double Q1 = N1 * tan_fp / R1;
    double Q2 = 0.5 * D2;
    double Q3 = (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C1 * C1
                - 9.0 * e2) * D4 / 24.0;
    double Q4 = (61.0 + 90.0 * T1 + 298.0 * C1 + 45.0 * T1 * T1
                - 252.0 * e2 - 3.0 * C1 * C1) * D6 / 720.0;

    double lat_rad = fp - Q1 * (Q2 - Q3 + Q4);

    // 计算经度修正
    double Q5 = D;
    double Q6 = (1.0 + 2.0 * T1 + C1) * D3 / 6.0;
    double Q7 = (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C1 * C1
                + 8.0 * e2 + 24.0 * T1 * T1) * D5 / 120.0;

    double lon_rad = lon0 + (Q5 - Q6 + Q7) / cos_fp;

    // 转换为度
    geo->latitude = coord_rad_to_deg(lat_rad);
    geo->longitude = coord_rad_to_deg(lon_rad);
    geo->altitude = 0.0;
    geo->datum = DATUM_TOKYO;

    // 如果需要 WGS84 坐标，进行基准转换
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

// ==================== 基准转换函数 ====================
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
    // 获取转换参数
    DatumTransform *params = &ctx->transforms[src->datum][target_datum];
    if (params->dx == 0.0 && params->dy == 0.0 && params->dz == 0.0 &&
            params->rx == 0.0 && params->ry == 0.0 && params->rz == 0.0 &&
            params->scale == 0.0)
    {
        // 没有转换参数，直接返回
        *dst = *src;
        dst->datum = target_datum;
        return COORD_SUCCESS;
    }
    // 获取源椭球体参数
    const Ellipsoid *src_ell = &ELLIPSOIDS[src->datum];
    const Ellipsoid *dst_ell = &ELLIPSOIDS[target_datum];
    // 将经纬度转换为地心直角坐标
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
    // 应用七参数转换
    double rx_rad = params->rx * ARC_SEC_TO_RAD;
    double ry_rad = params->ry * ARC_SEC_TO_RAD;
    double rz_rad = params->rz * ARC_SEC_TO_RAD;
    double scale_factor = 1.0 + params->scale * PPM_TO_SCALE;
    double X2 = params->dx + X * scale_factor + Y * rz_rad - Z * ry_rad;
    double Y2 = params->dy - X * rz_rad + Y * scale_factor + Z * rx_rad;
    double Z2 = params->dz + X * ry_rad - Y * rx_rad + Z * scale_factor;
    // 转换回大地坐标
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

// ==================== 测地线计算函数 ====================
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
    // 如果基准不同，需要转换
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
    // 如果基准不同，需要转换
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

// ==================== 基准转换工具函数 ====================
int coord_set_transform_params(CoordContext *ctx, MapDatum from, MapDatum to,
                               const DatumTransform *params)
{
    if (!ctx || from >= DATUM_MAX || to >= DATUM_MAX || !params)
    {
        return COORD_ERROR_INVALID_INPUT;
    }
    ctx->transforms[from][to] = *params;
    // 设置反向转换参数（正确计算七参数反向转换）
    if (from != to)
    {
        // 七参数反向转换的正确公式：
        // 反向转换需要通过旋转矩阵求逆来计算
        // 如果正向: X2 = T + s*R*X1
        // 则反向: X1 = (-1/s)*R^T*(X2-T) = (-1/s)*R^T*X2 + (1/s)*R^T*T
        //
        // 对于小角度近似（弧秒单位），反向参数约为：
        // 反向平移 = -R^T * T / (1+s)
        // 反向旋转 ≈ -旋转（但需要更精确计算）

        double s = params->scale * PPM_TO_SCALE;

        // 构建旋转矩阵 R
        // R = Rz(θz) * Ry(θy) * Rx(θx)
        // 对于小角度，可以近似为：
        // R ≈ [ 1    -rz    ry
        //       rz     1   -rx
        //      -ry    rx     1 ]

        // 计算反向转换参数
        // 反向尺度因子
        ctx->transforms[to][from].scale = -params->scale;

        // 反向旋转参数（近似，对于小角度）
        ctx->transforms[to][from].rx = -params->rx;
        ctx->transforms[to][from].ry = -params->ry;
        ctx->transforms[to][from].rz = -params->rz;

        // 反向平移参数: T_back = -(dx,dy,dz) / (1+s)
        // 更精确的计算需要考虑旋转矩阵的转置
        double factor = 1.0 / (1.0 + s);
        ctx->transforms[to][from].dx = -params->dx * factor;
        ctx->transforms[to][from].dy = -params->dy * factor;
        ctx->transforms[to][from].dz = -params->dz * factor;

        // 对于小角度，可以加入旋转修正项
        // 修正: T_back ≈ -(T + R×T) / (1+s)
        // 这里使用一阶近似
        double dx_corr = params->ry * params->dz - params->rz * params->dy;
        double dy_corr = params->rz * params->dx - params->rx * params->dz;
        double dz_corr = params->rx * params->dy - params->ry * params->dx;
        // 将弧秒转换为弧度的修正
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

// ==================== 椭球体工具函数 ====================
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
    // 重新初始化GeographicLib测地线对象
    geod_init(ctx->geod, a, f);
    return COORD_SUCCESS;
}

// ==================== 错误处理函数 ====================
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

// ==================== 格式转换主函数 ====================
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
    // 转换为目标基准
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
    // 根据目标格式进行格式化
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