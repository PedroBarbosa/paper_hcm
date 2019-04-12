import os.path
import argparse
from datasets import preprocess
from filters import filters
from thresholds import threshold_list
from tools import apply_tool_predictions
from plots.performance_comparison import generate_performance_comparison
from plots.heatmap import generate_heatmap

def main():
    """ Generates all the tables and plots for the paper """
    parser = argparse.ArgumentParser(description='Script to trigger the full benchmark analysis')
    parser.add_argument("--out_dir", help='Path to store all the output results. Default: "../../outdir"')
    parser.add_argument("--location", default="HGVSc", choices=("HGVSc","Consequence"),
                        help='VCF field to extract location of the variant')
    referencesetmode = parser.add_argument_group('Tools performance on reference variant sets')
    referencesetmode.add_argument("--limitAnalysis", metavar='dataset', type=lambda s: list(map(str, s.split(","))), default=[],
                        help='Set the analysis to be performed to the given datasets. Default: Run hcm. '
                             'Choices: [hcm]')

    args = parser.parse_args()

    datasets = ['hcm']
    if args.out_dir :
        PAPER_LOCATION = args.out_dir
        if not os.path.isdir(args.out_dir):
            os.mkdir(args.out_dir)
            os.mkdir(os.path.join(args.out_dir,"figures"))
    else:
        PAPER_LOCATION = "../../outdir"

    if all(x in datasets for x in args.limitAnalysis):

        for analysis in datasets:
            data = preprocess(args.location,analysis)
            generate_performance_comparison(data, filters, threshold_list, analysis, PAPER_LOCATION)
            df = apply_tool_predictions(data, threshold_list)
            generate_heatmap(df, filters, threshold_list, analysis, PAPER_LOCATION)

    else:
        print("Please limit your analysis to one (or more) of the following dataset options:\t{}".format(datasets))
        exit(1)


if __name__ == '__main__':
    main()