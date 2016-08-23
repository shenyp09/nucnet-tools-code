#!/usr/bin/python
# Imports
from xml.dom.minidom import parse
import xml.dom.minidom
import sys
import matplotlib.pyplot as plt

# USAGE
if len(sys.argv) < 4:
    print "\nUsage: python plot_mass.py INPUT_XML_FILE property nuclide_list",\
          "\nFor the property, it can have up to two tags.",\
          "\nExample1: python plot_mass.py my_output.xml time h1 he4",\
          "\nExample2: python plot_mass.py my_new_output.xml \"exposure,n\" fe56 fe57 fe58\n"
    exit()

# Set if plot abundance instead of mass fractions
ifabundance = "yes"

# Plot settings
# Output file name
output_figure_name = "plot_mass.png"
# Plot title on top of the figure
title = ""
# X label
xlabel = r"$\tau_n$(mb$^{-1}$)"
# Y label
ylabel = "Abundance per nucleon"
# legend location (0 is auto best)
# reference: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.legend
legend_location = 0
# Ranges. ifaxis == "no" applies default data range. [xbegin,xend,ybegin,yend]
ifaxis = "yes"
axis = [0, 0.8, 1.e-8, 1.e-1]
# Log. Toggle to "yes" if logarithmic axis is needed.
ifxlog = "no"
ifylog = "yes"
# symlog is used which will handle 0 values as linear scale locally.
ifxsymlog = "no"
ifysymlog = "no"

# Open XML document using minidom parser
DOMTree = xml.dom.minidom.parse(sys.argv[1])
collection = DOMTree.documentElement

# Read inputs
# Read property to plot on x axis
plot_property = sys.argv[2]
plot_property_list = []
plot_property_list.append(plot_property)
# For property with tags:
if not plot_property.isalnum():
    plot_property_list = plot_property.split(",")
if len(plot_property_list) > 3:
    print "\nToo many property tags (at most 2)!\n"
    exit()
# Read nuclide list
nuclide_list = []
for i in range(3, len(sys.argv)):
    nuclide_list.append(sys.argv[i])

# Create nuclide names for plotting
print_names = []
mass_numbers = []
for nuclide in nuclide_list:
    letters = nuclide[0].upper()
    numbers = ""
    for i in range(1, len(nuclide)):
        if nuclide[i].isalpha():
            letters += nuclide[i]
        if nuclide[i].isdigit():
            numbers += nuclide[i]
    name = r"^{%s}\rm{%s}" % (numbers, letters)
    print_names.append(name)
    mass_numbers.append(int(numbers))

# Create a dictionary
dict = {}
plot_properties = []
dict[plot_property] = plot_properties
for i in range(len(nuclide_list)):
    mass_list = []
    dict[nuclide_list[i]] = mass_list

# Get all the zones in the collection
zones = collection.getElementsByTagName("zone")

# Assign data for each zone
for zone in zones:

    zone_number = int(zone.getAttribute("label1"))

    # Get property lists
    props = zone.getElementsByTagName("property")
    for prop in props:
        if len(plot_property_list) == 1:
            if prop.getAttribute("name") == plot_property_list[0]:
                dict[plot_property].append(prop.firstChild.data)
        if len(plot_property_list) == 2:
            if prop.getAttribute("name") == plot_property_list[0] and\
                    prop.getAttribute("tag1") == plot_property_list[1]:
                dict[plot_property].append(prop.firstChild.data)

    # Get mass fractions
    zone_nuclides = zone.getElementsByTagName("nuclide")
    for nuclide in zone_nuclides:
        nuclide_name = nuclide.getAttribute("name")
        if nuclide_name in nuclide_list:
            nuclide_mass = nuclide.getElementsByTagName("x")[0].firstChild.data
            dict[nuclide_name].append(nuclide_mass)

    # Append mass list with 0 if not existing
    for nuclide in nuclide_list:
        if len(dict[nuclide]) < zone_number:
            dict[nuclide].append(0)

# Change from mass fraction to abundance per nucleon
if ifabundance == "yes":
    for i in range(len(nuclide_list)):
        for j in range(len(dict[nuclide_list[i]])):
            dict[nuclide_list[i]][j] = float(dict[nuclide_list[i]][j])\
                / mass_numbers[i]

# Plot
for i in range(len(nuclide_list)):
    plt.plot(dict[plot_property], dict[nuclide_list[i]],
             label=r"$%s$" % print_names[i])
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc=legend_location)
if ifaxis == "yes":
    plt.axis(axis)
if ifxlog == "yes":
    plt.xscale("log")
if ifylog == "yes":
    plt.yscale("log")
if ifxsymlog == "yes":
    plt.xscale("symlog")
if ifysymlog == "yes":
    plt.yscale("symlog")
plt.savefig(output_figure_name)
