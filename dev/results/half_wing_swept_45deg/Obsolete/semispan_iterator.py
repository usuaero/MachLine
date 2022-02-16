import fileinput
import json
import pressure_plotting_test as ppt

#This script is to test iterating over the pressure_plotting script for all locations along the wing

#import pressure_plotting

# for line in fileinput.input():
#     process(line)


inputfile = "Swept_half_wing_conditions_input.json"

json_string = open(inputfile).read()
json_vals = json.loads(json_string)

if json_vals["plots"]["show_all"]:
    for i in range(len(json_vals["geometry"]["Percent Semispan Locations"])):
        ppt.Swept_Plotting(json_vals["geometry"]["Percent Semispan Locations"][i]).get_data()
else:
    x = input("Enter percent semispan for analysis results. Options are 4.1, 8.2, 16.3, 24.5, 36.7, 51.0, 65.3, 89.8, 94.9:   ",)
    ppt.Swept_Plotting(x).get_data()



