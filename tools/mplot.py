import sys
import numpy as np

from Magics import macro as magics

def main():

    data =  magics.mgrib(grib_input_file_name = sys.argv[1],)

    projection = magics.mmap(
        subpage_map_library_area = "on",
        subpage_map_area_name    = "global",
        page_id_line             = "off"
    )

    coast = magics.mcoast()

    legend = magics.mlegend(legend_display_type   = 'continuous')

    # contour = magics.mcont(
    #     contour_automatic_setting = "ecmwf",
    #     legend = "on")

    output = magics.output(output_formats = ['png'],
        output_width = 1200,
        output_cairo_transparent_background = "on",
        output_name_first_page_number = "off",
        output_name = "magics")
 
    contour = magics.mcont( 
        contour                        = "on",
        contour_level_selection_type   = "interval",
        contour_shade                  = "on",
        contour_shade_method           = "area_fill",
        # contour_interval               = 0.1,
        # contour_shade_min_level        = -1.0,
        # contour_shade_max_level        = 1.0,
        contour_shade_colour_method    = "calculate",
        contour_shade_colour_direction = "clockwise",
        contour_highlight              = "off",
        contour_label                  = "off",
        contour_shade_max_level_colour = "red",
        contour_shade_min_level_colour = "blue",
        # legend                         = "on"
    )

    magics.plot(output, projection, data, contour, coast, legend)

if __name__ == "__main__":
    sys.exit(main())
