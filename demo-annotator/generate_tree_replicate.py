import dendropy, random, argparse

import muriqui as m

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="generate random test trees and annotations")

    parser.add_argument("-l", "--label", nargs=1, required=True, help="The label to use for the output files.")

    parser.add_argument("-n", "--ntax", nargs=1, type=int, required=True, help="The number of tips for the tree generator.")

    