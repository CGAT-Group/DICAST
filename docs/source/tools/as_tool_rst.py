import argparse
import pandas as pd
from rst_render import *

parser = argparse.ArgumentParser(description='Document splicing tools for DICAST in .rst format.')
parser.add_argument('--tools', metavar='-t', type=str, nargs='*',
                    help='Tools to document separated by space. None to document all tools')
parser.add_argument('--inputdir', metavar='-i', type=str, nargs='?', default="./",
                    help='Input directory with documentation tables. Default=directory where the script is executed')
parser.add_argument('--outdir', metavar='-o', type=str, nargs='?', default="./",
                    help='Output directory for .rst files Default=directory where the script is executed')
args = parser.parse_args()

tools = args.tools
inputdir = args.inputdir

as_tools_table = pd.read_csv(inputdir + "/as_tools.csv", sep="\t", header=0, index_col=0, dtype=str)
compatibility_table = pd.read_csv(inputdir + "/compatibility_table.csv", sep="\t", header=0, index_col=0, dtype=str)
required_files = pd.read_csv(inputdir + "/input_files.csv", sep="\t", header=0, index_col=0, dtype=str)
parameter_definitions = pd.read_csv(inputdir + "/parameter_description.csv", sep="\t", header=0, index_col=0, dtype=str)



def rendertool(toolname):
    print(f"Starting documentation for {toolname} ...")
    tool_parameters = pd.read_csv(inputdir + "/tool_specific/" + toolname + ".csv", sep="\t", header=0, dtype=str)
    tool_info = as_tools_table.loc[toolname]

    print(tool_parameters)

    version = tool_info["version"]
    manual = tool_info["manual_link"]
    licensename = tool_info["license"]
    license_link = tool_info["license_link"]
    officialname = tool_info["officialname"]

    file = open(args.outdir + "/" + toolname + ".rst", "w")

    file.write(comment("Links"))
    file.write(reference("manual", manual))
    file.write(replace("tool", officialname))
    file.write(title(officialname))

    # Opional: Warning
    if not pd.isna(tool_info["warning"]):
        file.write(warning(tool_info["warning"]))

    if not pd.isna(tool_info["short_description"]):
        file.write(tool_info["short_description"]+"\n\n")

    # Differential
    if tool_info["differential"]=="both":
        file.write(note("|tool| can be used to calculate differential splicing as well as only alternative-splicing events.\nIf you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.Otherwise set ``differential=0``."))

    # Sidebar

    sidebar_content = ""

    table = small_table([[bold("Toolname"), italic(toolname)],
                         [bold("Version"), italic(version)],
                         [bold("License"), inline_link(licensename, license_link)]])

    sidebar_content += table + "\n"
    sidebar_content += bold("Required Files") + "\n"


    files_list = [ref_link(file, f"{file}Splicing")  for file in required_files.columns.values if
                  required_files[file].loc[toolname] == "Yes"]

    sidebar_content += indent_content("* " +" , ".join(files_list), 1)

    if ("bam" in files_list):
        compatible_list = [doc_link("../splicing/"+tool) for tool in compatibility_table.index.values if
                       compatibility_table[tool].loc[toolname] == "Yes"]
        if len(compatible_list)>0:
            sidebar_content += bold("Compatible mapping tools") + "\n"
            sidebar_content += li(compatible_list, indent=1)

    sidebar_content += bold("Links") + "\n"
    sidebar_content += li(["|tool| `manual`_", "|tool| publication: " + inline_link(tool_info["publication_name"],
                                                                                    tool_info["publication_link"])],
                          indent=1)

    file.write(directive("sidebar", "|tool| Factsheet", sidebar_content, False))

    # Opional: Note
    if not pd.isna(tool_info["note"]):
        file.write(note(tool_info["note"]))


    # Parameters
    file.write(title("Parameters", symbol="^"))
    file.write(
        "These are the default parameters set in the :guilabel:`src/{0}/ENTRYPOINT.sh` script. If you want to change it "
        "you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.\n\n".format(
            toolname))


    for index, row in tool_parameters.iterrows():
        content = parameter_definitions.loc[row["use"]]["description"]
        if not pd.isna(parameter_definitions.loc[row["use"]]["code"]):
            if row["tag"].startswith("-"):
                tag=row["tag"] + " "
            else:
                tag=""
            content += "\n\n"+codeblock( tag + parameter_definitions.loc[row["use"]]["code"].replace("tool",toolname))[:-1]
        file.write(param_defin(row["tag"], content, indent=1)[:-1])

    #Issue list

    if not pd.isna(tool_info["issues"]):
        file.write("\n")
        file.write(title("Known Issues", symbol="^"))
        file.write(li(tool_info["issues"].split(";")))

    file.close()
    print(f"Documentation for {toolname} done")


if not tools:
    # render all tools
    tools = as_tools_table.index.values
    for tool in tools:
        rendertool(tool)
elif type(tools) == str:
    # one tool
    tool = tools
    rendertool(tool)
else:
    for tool in tools:
        rendertool(tool)
