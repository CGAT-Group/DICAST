import argparse
import pandas as pd
from rst_render import *

parser = argparse.ArgumentParser(description='Document mapping tools for DICAST in .rst format.')
parser.add_argument('--tools', metavar='-t', type=str, nargs='*',
                    help='Tools to document separated by space. None to document all tools')
parser.add_argument('--inputdir', metavar='-i', type=str, nargs='?', default="./",
                    help='Input directory with documentation tables. Default=directory where the script is executed')
parser.add_argument('--outdir', metavar='-o', type=str, nargs='?', default="./",
                    help='Output directory for .rst files Default=directory where the script is executed')
args = parser.parse_args()

tools = args.tools
inputdir = args.inputdir

mapping_tools_table = pd.read_csv(inputdir + "/mapping_tools.csv", sep="\t", header=0, index_col=0, dtype=str)
compatibility_table = pd.read_csv(inputdir + "/compatibility_table_comma.csv", sep=",", header=0, index_col=0, dtype=str)
required_files = pd.read_csv(inputdir + "/input_files.csv", sep="\t", header=0, index_col=0, dtype=str)
parameter_definitions = pd.read_csv(inputdir + "/parameter_description.csv", sep="\t", header=0, index_col=0, dtype=str)



def rendertool(toolname):
    print(f"Starting documentation for {toolname} ...")
    tool_parameters = pd.read_csv(inputdir + "/tool_specific/" + toolname + ".csv", sep="\t", header=0, dtype=str)
    tool_info = mapping_tools_table.loc[toolname]

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

    # Sidebar

    sidebar_content = ""

    table = small_table([[bold("Toolname"), italic(toolname)],
                         [bold("Version"), italic(version)],
                         [bold("License"), inline_link(licensename, license_link)]])

    sidebar_content += table + "\n"
    sidebar_content += bold("Required Files") + "\n"

    files_index_list = [ref_link(file, f"{file}Mapping") for file in required_files.columns.values if
                        required_files[file].loc[toolname] == "index"]
    files_list = [ref_link(file, f"{file}Mapping") for file in required_files.columns.values if
                  required_files[file].loc[toolname] == "Yes"]

    if len(files_index_list) > 0:
        sidebar_content += indent_content("* " + italic("For indexing only:") + " " + ", ".join(files_index_list), 1)
    sidebar_content += indent_content("* " + italic("For mapping:") + " " + ", ".join(files_list), 1)

    compatible_list = [doc_link("../splicing/"+tool) for tool in compatibility_table.columns.values if
                       compatibility_table[tool].loc[toolname] == "Yes"]
    if len(compatible_list)>0:
        sidebar_content += bold("Compatible splicing tools") + "\n"
        sidebar_content += li(compatible_list, indent=1)

    sidebar_content += bold("Links") + "\n"
    sidebar_content += li(["|tool| `manual`_", "|tool| publication: " + inline_link(tool_info["publication_name"],
                                                                                    tool_info["publication_link"])],
                          indent=1)

    file.write(directive("sidebar", "|tool| Factsheet", sidebar_content, False))

    # Opional: Note
    if not pd.isna(tool_info["note"]):
        file.write(note(tool_info["note"]))


    # Indexing
    file.write(title("Indexing", symbol="^"))

    file.write(note(bold(
        "Indexing might take some time") + " but only has to be run once per fasta file. Make sure to reuse already computed indices if possible."))
    file.write("DICAST will check if :guilabel:`{0}` exists. If there is no index it will be automatically built. "
               "If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.\n"
               "If you want to use your own precomputed index file copy it to :guilabel:`index/{1}-index/` and make sure the index is complete and named appropriately and according to the parameters set in the config files.\n"
               "We recommend including the name of the fasta file in the index name to avoid overwriting. Per default this is already the case and **no parameter changes are needed**.\n\n".format(
        tool_info["indexname"], toolname))

    # Parameters
    file.write(title("Parameters", symbol="^"))
    file.write(
        "These are the default parameters set in the :guilabel:`src/{0}/ENTRYPOINT.sh` script. If you want to change it you can do this in the ENTRYPOINT script directly. Please refer to the |tool| `manual`_.\n\n".format(
            toolname))

    for index, row in tool_parameters.iterrows():
        content = parameter_definitions.loc[row["use"]]["description"] + "\n\n"
        if not pd.isna(parameter_definitions.loc[row["use"]]["code"]):
            if row["tag"].startswith("-"):
                tag=row["tag"] + " "
            else:
                tag=""
            content += codeblock( tag + parameter_definitions.loc[row["use"]]["code"].replace("tool",toolname))
        file.write(param_defin(row["tag"], content, indent=1))

    #Issue list

    if not pd.isna(tool_info["issues"]):
        file.write("\n")
        file.write(title("Known Issues", symbol="^"))
        file.write(li(tool_info["issues"].split(";")))

    file.close()
    print(f"Documentation for {toolname} done")


if not tools:
    # render all tools
    tools = mapping_tools_table.index.values
    for tool in tools:
        rendertool(tool)
elif type(tools) == str:
    # one tool
    tool = tools
    rendertool(tool)
else:
    for tool in tools:
        rendertool(tool)
