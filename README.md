# Coordinate Datum Transform Library

A high-precision coordinate transformation library supporting multiple coordinate formats and map datums.

## Features

- **Multiple Coordinate Formats**: DD, DMM, DMS, UTM, MGRS, British Grid, Japan Grid
- **Map Datums**: WGS84, NAD83, NAD27, ED50, Tokyo, OSGB36
- **Seven-Parameter Helmert Transformations**
- **Geodesic Calculations** using GeographicLib
- **Correct MGRS 6-set cycle algorithm**

## Files

- `coord_datum_transform.h` - Header file
- `coord_datum_transform.c` - Core implementation
- `geodesic.h` / `geodesic.c` - GeographicLib
- `test_coord_datum_transform.c` - Test suite
- `markdown` - Documentation notes

## Compilation (Windows MSYS2)

```bash
export PATH="D:\msys64\ucrt64\bin:$PATH"
gcc -c coord_datum_transform.c -o coord_datum_transform.o
gcc -c geodesic.c -o geodesic.o
gcc test_coord_datum_transform.c coord_datum_transform.o geodesic.o -o test_transform.exe -lm
```

## Author

**Wenbing Wang (Go-WB)**  
Email: wbwang@mail.ustc.edu.cn  
GitHub: https://github.com/Go-WB

## License

MIT License
