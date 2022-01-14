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

for i in range(len(json_vals["geometry"]["Percent Semispan Locations"])):
    ppt.Swept_Plotting.get_data((json_vals["geometry"]["Percent Semispan Locations"][i]))