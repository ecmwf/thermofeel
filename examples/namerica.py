# (C) Copyright 1996- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import sys

from Magics import macro as magics


def main():

    data = magics.mgrib(
        grib_input_file_name=sys.argv[1],
    )

    namerica = magics.mmap(
        subpage_lower_left_latitude=11.00,
        subpage_lower_left_longitude=-120.00,
        subpage_upper_right_latitude=47.00,
        subpage_upper_right_longitude=-32.00,
        subpage_map_vertical_longitude=-100.00,
        subpage_map_projection="polar_stereographic",
        page_x_length=20.0,
        page_y_length=20.0,
        subpage_x_length=20.0,
        subpage_y_length=20.0,
        page_x_position=0.5,
        page_y_position=1,
        layout="positional",
        super_page_x_length=20.00,
        super_page_y_length=20.00,
        page_id_line="off",
    )

    line_text_colour = "charcoal"
    text_height = 0.5
    coast = magics.mcoast(
        map_coastline_resolution="medium",
        map_grid="on",
        map_grid_colour=line_text_colour,
        map_grid_latitude_increment=20.00,
        map_grid_longitude_increment=20.00,
        map_label_height=text_height,
        map_label_right="off",
        map_label_left="off",
        map_label_latitude_frequency=1,
        map_label_longitude_frequency=2,
        map_label_colour=line_text_colour,
        map_label_blanking="off",
    )

    legend = magics.mlegend(
        legend_display_type="continuous",
        legend_box_mode="positional",
        legend_box_x_position=17.00,
        legend_box_y_position=0.15,
        legend_box_x_length=2.0,
        legend_box_y_length=16.70,
        legend_automatic_position="right",
        legend_title="off",
        legend_text_font_size=text_height,
        legend_text_colour=line_text_colour,
        legend_label_frequency=10,
        legend_user_minimum="on",
        legend_user_minimum_text="< -50",
        legend_user_maximum="on",
        legend_user_maximum_text="> 50",
        legend_entry_border="off",
    )

    # utci_text = "Universal Thermal Climate Index"
    wbgt_text = "Wet Bulb Globe Temperature"

    title = magics.mtext(
        text_lines=[wbgt_text, "Valid date: 2 July 2021"],
        text_justification="centre",
        text_font_size=text_height + 0.2,
        text_colour=line_text_colour,
        text_mode="automatic",
    )

    output_namerica = magics.output(
        output_formats=["png"],
        output_width=800,
        output_cairo_transparent_background="on",
        output_name_first_page_number="off",
        output_name="/home/ma/mamv/2021/wbgt_usa_20210702",
    )

    contour = magics.mcont(
        contour="off",
        contour_level_selection_type="interval",
        contour_shade="on",
        contour_shade_method="area_fill",
        contour_interval=1.0,
        contour_shade_min_level=-50.0,
        contour_shade_max_level=50.0,
        contour_interpolation_ceiling=49.99,
        contour_interpolation_floor=-49.99,
        contour_label="off",
        contour_shade_colour_method="list",
        contour_shade_colour_list=[
            "#053061",
            "#073466",
            "#09376b",
            "#0b3b70",
            "#0e3f75",
            "#10427a",
            "#12467e",
            "#154a83",
            "#174e88",
            "#1a528c",
            "#1c5691",
            "#1f5a95",
            "#225e99",
            "#24629e",
            "#2766a2",
            "#2a6aa6",
            "#2c6faa",
            "#2f73ad",
            "#3277b1",
            "#357cb4",
            "#3880b8",
            "#3a85bb",
            "#3d89be",
            "#408ec0",
            "#4392c3",
            "#4a97c5",
            "#529bc7",
            "#5a9fc9",
            "#61a3cc",
            "#68a7ce",
            "#70abd0",
            "#77afd2",
            "#7eb3d4",
            "#85b7d7",
            "#8cbbd9",
            "#93bfdb",
            "#9ac3dd",
            "#a1c7df",
            "#a8cbe1",
            "#afcfe3",
            "#b6d3e5",
            "#bdd7e7",
            "#c4dbe9",
            "#cbdfeb",
            "#d2e3ed",
            "#dae7ef",
            "#e1ebf1",
            "#e8eff3",
            "#f0f3f5",
            "#f7f7f7",
            "#f7f7f7",
            "#f8f1ef",
            "#f8ebe6",
            "#f9e6de",
            "#f9e0d6",
            "#f9dace",
            "#f9d5c7",
            "#f8cfbf",
            "#f8c9b8",
            "#f7c4b1",
            "#f6beaa",
            "#f5b8a3",
            "#f4b39c",
            "#f3ad95",
            "#f1a78f",
            "#f0a289",
            "#ee9c83",
            "#ec977d",
            "#ea9177",
            "#e88b71",
            "#e6866c",
            "#e48067",
            "#e27b62",
            "#df755d",
            "#dd6f58",
            "#da6954",
            "#d86450",
            "#d55e4c",
            "#d15848",
            "#ce5345",
            "#ca4e42",
            "#c64940",
            "#c2443d",
            "#bd3f3a",
            "#b93a38",
            "#b43536",
            "#af3033",
            "#ab2c31",
            "#a5272f",
            "#a0232d",
            "#9b1e2c",
            "#961a2a",
            "#901628",
            "#8a1227",
            "#850e25",
            "#7f0a24",
            "#790623",
            "#730321",
            "#6d0120",
            "#67001f",
        ],
        legend="on",
    )

    magics.plot(output_namerica, namerica, data, contour, coast, legend, title)


if __name__ == "__main__":
    sys.exit(main())
