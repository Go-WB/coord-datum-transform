# 坐标转换系统

基于GeographicLib的C语言坐标转换系统，支持多种坐标格式和地图基准的两两组合转换。

## 功能特性

### 支持的坐标格式
- DD.ddddd° (十进制度)
- DD°MM.mmm' (度分)
- DD°MM'SS" (度分秒)
- UTM坐标
- MGRS坐标
- British Grid (英国国家网格)
- Japan Grid (日本网格)

### 支持的地图基准
- WGS84 (World Geodetic System 1984)
- MGRS Grid (基于WGS84)
- UTM Grid (基于WGS84)
- NAD83 (North American Datum 1983)
- NAD27 (North American Datum 1927)
- ED50 (European Datum 1950)

### 核心功能
1. **统一转换接口**：输入为WGS84 DD.ddddd°格式，输出为指定格式和基准
2. **基准转换**：使用七参数法进行基准转换
3. **坐标格式转换**：支持所有格式之间的相互转换
4. **线程安全**：每个上下文独立，支持