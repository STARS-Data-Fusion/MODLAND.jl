module MODLAND

import ArchGDAL as AG
import ArchGDAL
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
using DataFrames
using Rasters
using DimensionalData.Dimensions.LookupArrays
import JSON

# boundaries of sinusodial projection
UPPER_LEFT_X_METERS = -20015109.355798
UPPER_LEFT_Y_METERS = 10007554.677899
LOWER_RIGHT_X_METERS = 20015109.355798
LOWER_RIGHT_Y_METERS = -10007554.677899

# size across (width or height) of any equal-area sinusoidal target
TILE_SIZE_METERS = 1111950.5197665554

# boundaries of MODIS land grid
TOTAL_ROWS = 18
TOTAL_COLUMNS = 36

WGS84 = ProjString("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
SINUSOIDAL = ProjString("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

function latlon_to_sinusoidal(lat::Float64, lon::Float64)::Tuple{Float64,Float64}
    if lat < -90 || lat > 90
        error("latitude ($(lat)) out of bounds")
    end

    if lon < -180 || lon > 180
        error("longitude ($(lon))) out of bounds")
    end

    point_latlon = AG.createpoint(lon, lat)
    point_sinusoidal = AG.reproject(point_latlon, WGS84, SINUSOIDAL)
    x = AG.getx(point_sinusoidal, 0)
    y = AG.gety(point_sinusoidal, 0)

    return x, y
end

export latlon_to_sinusoidal

# MODIS land target indices for target containing sinusoidal coordinate
function sinusoidal_to_MODLAND(x::Float64, y::Float64)::String
    if x < UPPER_LEFT_X_METERS || x > LOWER_RIGHT_X_METERS
        error("sinusoidal x coordinate ($(x)) out of bounds")
    end

    if y < LOWER_RIGHT_Y_METERS || y > UPPER_LEFT_Y_METERS
        error("sinusoidal y ($(y))) coordinate out of bounds")
    end

    horizontal_index = Int(floor((x - UPPER_LEFT_X_METERS) / TILE_SIZE_METERS))
    vertical_index = Int(floor((-1 * (y + LOWER_RIGHT_Y_METERS)) / TILE_SIZE_METERS))

    if horizontal_index == TOTAL_COLUMNS
        horizontal_index -= 1
    end

    if vertical_index == TOTAL_ROWS
        vertical_index -= 1
    end

    tile = "h$(lpad(horizontal_index, 2, '0'))v$(lpad(vertical_index, 2, '0'))"

    return tile
end

export sinusoidal_to_MODLAND

function latlon_to_MODLAND(lat, lon)
    x, y = latlon_to_sinusoidal(lat, lon)
    tile = sinusoidal_to_MODLAND(x, y)

    return tile
end

export latlon_to_MODLAND

function MODLAND_tiles_in_polygon(poly::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon25D})::Set{String}
    Set([latlon_to_MODLAND(lat, lon) for (lon, lat) in JSON.parse(AG.toJSON(poly))["coordinates"][1]])
end

function MODLAND_tiles_in_polygon(poly::AG.IGeometry{AG.wkbPolygon})::Set{String}
    Set([latlon_to_MODLAND(lat, lon) for (lon, lat) in JSON.parse(AG.toJSON(poly))["coordinates"][1]])
end

export MODLAND_tiles_in_polygon

function sinusoidal_tile_dims(h::Int, v::Int, tile_width_cells::Int)::Tuple{X,Y}
    cell_size = TILE_SIZE / tile_width_cells
    tile_left_x = GLOBE_UPPER_LEFT_X + h * TILE_SIZE
    tile_right_x = GLOBE_UPPER_LEFT_X + (h + 1) * TILE_SIZE - cell_size
    tile_upper_y = GLOBE_LOWER_RIGHT_Y + (TOTAL_ROWS - v) * TILE_SIZE - cell_size
    tile_lower_y = GLOBE_LOWER_RIGHT_Y + (TOTAL_ROWS - 1 - v) * TILE_SIZE
    sampling = Intervals(Start())
    x_dim = X(Projected(LinRange(tile_left_x, tile_right_x, tile_width_cells), order=ForwardOrdered(), span=Regular(cell_size), sampling=sampling, crs=SINUSOIDAL_CRS))
    y_dim = Y(Projected(LinRange(tile_upper_y, tile_lower_y, tile_width_cells), order=ReverseOrdered(), span=Regular(-cell_size), sampling=sampling, crs=SINUSOIDAL_CRS))
    dims = (x_dim, y_dim)

    return dims
end

export sinusoidal_tile_dims

end
