# Convert a directory of pdb files into png images using PyMOL
# This script should be run like so: pymol -c code/generate_png_from_pdb.pyy <csv_filepath>
# The CSV file should have a header with columns:
# - approved_symbol - Used to create output png filename.
# - uni_prot_id_supplied_by_uni_prot - Used to determine input pdb.gz filepath.
# The CSV file should be in the same directory as the compressed pdb files.
# The output png file will be created in the same directory as the CSV file.
import os
import os.path
import csv
import sys
import pymol
from pymol import cmd
import gzip
import shutil
PDB_PREFIX = "AF"
PDB_SUFFIX = "F1-model_v2.pdb.gz"


def get_symbol_paths_ary(csv_filepath):
    """
    Given a path to a CSV file containing two columns (approved_symbol, uni_prot_id_supplied_by_uni_prot)
    return an array of (approved_symbol, compressed_pdb_path, pdb_path, png_path)
    for each row of the CSV that has a valid UniProt ID
    """
    dirname = os.path.dirname(csv_filepath)
    symbol_paths_ary = []
    with open(csv_filepath) as csvfile:
       for row in csv.DictReader(csvfile):
            approved_symbol, uniprot_id = row['approved_symbol'], row['uni_prot_id_supplied_by_uni_prot']
            if uniprot_id:
                symbol_paths_ary.append(get_symbol_paths(dirname, approved_symbol, uniprot_id))
    return  symbol_paths_ary


def get_symbol_paths(dirname, approved_symbol, uniprot_id):
    compressed_pdb_path = create_compressed_pdb_path(dirname, uniprot_id)
    outputdir = os.path.dirname(compressed_pdb_path)
    pdb_path = get_extracted_pdb_path(outputdir, approved_symbol)
    png_path = get_png_path(outputdir, approved_symbol)
    return approved_symbol, compressed_pdb_path, pdb_path, png_path


def create_compressed_pdb_path(directory, uniprot_id):
    # Given  directory (eg. /tmp) and a uniprot_id (eg. P20930)
    # return a path to where the compressed pdb file should exist (eg /tmp/AF-P20930-F1-model_v2.pdb.gz)
    return "{}/{}-{}-{}".format(directory, PDB_PREFIX, uniprot_id, PDB_SUFFIX)


def get_extracted_pdb_path(outputdir, approved_symbol):
    return "{}/{}.pdb".format(outputdir, approved_symbol)


def get_png_path(outputdir, approved_symbol):
    return "{}/{}.png".format(outputdir, approved_symbol)


def extract_pdb(compressed_pdb_path, pdb_path):
    with gzip.open(compressed_pdb_path, 'rb') as infile:
        with open(pdb_path, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)


def convert_pdb_to_png(pdb_source_path, png_destination_path):
    """
    Use pymol to load and render pdb_source_path into png_destination_path
    """
    cmd.reinitialize()
    cmd.load(pdb_source_path)
    cmd.remove('solvent')
    cmd.show_as('cartoon')
    cmd.orient()
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_trace_mode', 2)
    cmd.ray(antialias=2)
    cmd.png(png_destination_path, dpi=300)
    print("Created {}.".format(png_destination_path))


def convert_pdb_files_to_pngs(symbol_paths_ary):
    for approved_symbol, compressed_pdb_path, pdb_path, png_path in symbol_paths_ary:
        if os.path.exists(compressed_pdb_path):
            if not os.path.exists(png_path):  # skip if png already exists
                extract_pdb(compressed_pdb_path, pdb_path)
                convert_pdb_to_png(pdb_path, png_path)
        else:
            print("ERROR: File not found: {} for symbol: {}".format(compressed_pdb_path, approved_symbol))


def process_csv_file(csv_filepath):
    print("Reading {}.".format(csv_filepath))
    symbol_paths_ary = get_symbol_paths_ary(csv_filepath)
    print("Finish launching pymol")
    pymol.finish_launching()
    print("Converting pdb files to png")
    convert_pdb_files_to_pngs(symbol_paths_ary)
    print("Done converting pdb files to png")
    cmd.quit()


def main():
    # the script input arguments start at position 3 when run by pymol
    process_csv_file(csv_filepath=sys.argv[3])


# Only run main if we are not running tests
# Standard python __main__ logic doesn't work with pymol
if not os.environ.get('RUNNING_TESTS'):
    main()
