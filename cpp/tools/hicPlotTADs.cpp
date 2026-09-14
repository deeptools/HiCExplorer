// hicexplorer/hicPlotTADs.py is pygenometracks.plotTracks.main(args), and it
// stays that delegation (cpp/PLAN.md tier 7: "(a) shell, unchanged").
//
// This binary declares pyGenomeTracks 3.9's parser with hicx::cli, so
// --help-json, the GUI form and the argparse equality of PLAN 10.1 hold for
// the tool (cpp/scripts/tool_specs.py compares it with
// pygenometracks.plotTracks.parse_arguments). It checks the command line as
// argparse does, including opening --tracks and --BED, and hands the tokens to
// plot/hicexplorer_plot/hicPlotTADs.py, which calls plotTracks.main with them.
// There is nothing to compute: pyGenomeTracks reads every track itself.

#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kDescription =
    "Plots genomic tracks on specified region(s). Citations :\nRamirez et al.  High-resolution "
    "TADs reveal DNA sequences underlying genome organization in flies. Nature Communications "
    "(2018) doi:10.1038/s41467-017-02525-w\nLopez-Delisle et al.  pyGenomeTracks: reproducible "
    "plots for multivariate genomic datasets. Bioinformatics (2020) "
    "doi:10.1093/bioinformatics/btaa692";

const char* const kUsage =
    "usage: hicPlotTADs --tracks tracks.ini --region chr1:1000000-4000000 -o image.png\n";

const char* const kHelp =
    "\n"
    "Plots genomic tracks on specified region(s). Citations : Ramirez et al. High-\n"
    "resolution TADs reveal DNA sequences underlying genome organization in flies.\n"
    "Nature Communications (2018) doi:10.1038/s41467-017-02525-w Lopez-Delisle et\n"
    "al. pyGenomeTracks: reproducible plots for multivariate genomic datasets.\n"
    "Bioinformatics (2020) doi:10.1093/bioinformatics/btaa692\n"
    "\n"
    "options:\n"
    "  -h, --help            show this help message and exit\n"
    "  --tracks TRACKS       File containing the instructions to plot the tracks.\n"
    "                        The tracks.ini file can be genarated using the\n"
    "                        `make_tracks_file` program.\n"
    "  --region REGION       Region to plot, the format is chr:start-end\n"
    "  --BED BED             Instead of a region, a file containing the regions to\n"
    "                        plot, in BED format, can be given. If this is the\n"
    "                        case, multiple files will be created. It will use the\n"
    "                        value of --outFileName as a template and put the\n"
    "                        coordinates between the file name and the extension.\n"
    "  --width WIDTH         figure width in centimeters (default is 40)\n"
    "  --plotWidth PLOTWIDTH\n"
    "                        width in centimeters of the plotting (central) part\n"
    "  --height HEIGHT       Figure height in centimeters. If not given, the figure\n"
    "                        height is computed based on the heights of the tracks.\n"
    "                        If given, the track height are proportionally scaled\n"
    "                        to match the desired figure height.\n"
    "  --title TITLE, -t TITLE\n"
    "                        Plot title\n"
    "  --outFileName OUTFILENAME, -out OUTFILENAME\n"
    "                        File name to save the image, file prefix in case\n"
    "                        multiple images are stored\n"
    "  --fontSize FONTSIZE   Font size for the labels of the plot (default is 0.3 *\n"
    "                        figure width)\n"
    "  --dpi DPI             Resolution for the image in case the ouput is a raster\n"
    "                        graphics image (e.g png, jpg) (default is 72)\n"
    "  --trackLabelFraction TRACKLABELFRACTION\n"
    "                        By default the space dedicated to the track labels is\n"
    "                        0.05 of the plot width. This fraction can be changed\n"
    "                        with this parameter if needed.\n"
    "  --trackLabelHAlign {left,right,center}\n"
    "                        By default, the horizontal alignment of the track\n"
    "                        labels is left. This alignemnt can be changed to right\n"
    "                        or center.\n"
    "  --decreasingXAxis     By default, the x-axis is increasing. Use this option\n"
    "                        if you want to see all tracks with a decreasing\n"
    "                        x-axis.\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the command line is checked here and passed to pyGenomeTracks'\n"
    "plotTracks through the hicexplorer_plot drawing layer (HICX_PLOT_PYTHON names\n"
    "the interpreter). The C++-only option --plotData FILE writes the checked\n"
    "command line as JSON to FILE instead of plotting.\n";

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("hicPlotTADs", kDescription);
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& options = parser.group("options");
    options.add({"-h", "--help"}).action(cli::Action::Help).help("show this help message and exit");
    options.add({"--tracks"})
        .file_type("r")
        .required()
        .input({"ini"})
        .help("File containing the instructions to plot the tracks. The tracks.ini file can be "
              "genarated using the `make_tracks_file` program.");
    cli::MutuallyExclusiveGroup& region = parser.mutually_exclusive(options, true);
    region.add({"--region"}).help("Region to plot, the format is chr:start-end");
    region.add({"--BED"})
        .file_type("r")
        .input({"bed"})
        .help("Instead of a region, a file containing the regions to plot, in BED format, can be "
              "given. If this is the case, multiple files will be created. It will use the value "
              "of --outFileName as a template and put the coordinates between the file name and "
              "the extension.");
    cli::MutuallyExclusiveGroup& width = parser.mutually_exclusive(options, false);
    width.add({"--width"})
        .type("float")
        .default_value(std::int64_t{40})
        .help("figure width in centimeters (default is 40)");
    width.add({"--plotWidth"})
        .type("float")
        .help("width in centimeters of the plotting (central) part");
    options.add({"--height"})
        .type("float")
        .help("Figure height in centimeters. If not given, the figure height is computed based on "
              "the heights of the tracks. If given, the track height are proportionally scaled to "
              "match the desired figure height.");
    options.add({"--title", "-t"}).required(false).help("Plot title");
    options.add({"--outFileName", "-out"})
        .required()
        .output({"png", "pdf", "svg"}, "prefix")
        .help("File name to save the image, file prefix in case multiple images are stored");
    options.add({"--fontSize"})
        .type("float")
        .help("Font size for the labels of the plot (default is 0.3 * figure width)");
    options.add({"--dpi"})
        .type("int")
        .default_value(std::int64_t{72})
        .help("Resolution for the image in case the ouput is a raster graphics image (e.g png, "
              "jpg) (default is 72)");
    options.add({"--trackLabelFraction"})
        .default_value(0.05)
        .type("float")
        .help("By default the space dedicated to the track labels is 0.05 of the plot width. This "
              "fraction can be changed with this parameter if needed.");
    options.add({"--trackLabelHAlign"})
        .default_value("left")
        .choices({"left", "right", "center"})
        .help("By default, the horizontal alignment of the track labels is left. This alignemnt "
              "can be changed to right or center.");
    options.add({"--decreasingXAxis"})
        .action(cli::Action::StoreTrue)
        .help("By default, the x-axis is increasing. Use this option if you want to see all "
              "tracks with a decreasing x-axis.");
    options.add({"--version"}).version("3.9");
    options.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The checked command line as JSON, without plotting (cpp/PLAN.md tier 7).")
        .help("Write the checked command line as JSON to this file and do not plot.");

    const cli::Namespace ns = parser.parse(argc, argv);
    if (const int refused = hicx::plot::preflight("hicPlotTADs", !ns.given("plotData")); refused != 0) {
        return refused;
    }

    // The tokens for plotTracks.main, without the C++-only option.
    std::vector<std::string> tokens;
    for (int i = 1; i < argc; ++i) {
        const std::string token = argv[i];
        if (token == "--plotData") {
            ++i;
            continue;
        }
        if (token.rfind("--plotData=", 0) == 0) {
            continue;
        }
        tokens.push_back(token);
    }

    plot::JsonObject data;
    data.add("argv", plot::json_strings(tokens));
    hicx::report_resource_usage("hicPlotTADs");
    return plot::draw("hicPlotTADs", data.str(), ns.opt_str("plotData"));
}
