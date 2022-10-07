import unittest
import sys
from unittest.mock import patch, Mock, mock_open

# Mock pymol so we do not need to install it for testing
sys.modules['pymol'] = Mock()

# Setup our path so we can import from code
sys.path.append("code")

from generate_png_from_pdb import get_symbol_paths, get_symbol_paths_ary, convert_pdb_files_to_pngs


class TestFunctions(unittest.TestCase):
    def test_get_symbol_paths(self):
        result = get_symbol_paths(dirname='/tmp', approved_symbol='WASH7P', uniprot_id='123')
        self.assertEqual(result, ('WASH7P',  # approved_symbol
                                  '/tmp/AF-123-F1-model_v2.pdb.gz',  # compressed_pdb_path
                                  '/tmp/WASH7P.pdb',  # pdb_path
                                  '/tmp/WASH7P.png'))  # png_path

    @patch('generate_png_from_pdb.csv')
    def test_get_symbol_paths_ary(self, mock_csv):
        csv_dict_data = [
            {"approved_symbol": "WASH1P", "uni_prot_id_supplied_by_uni_prot": "111"},
            {"approved_symbol": "WASH2P", "uni_prot_id_supplied_by_uni_prot": "222"},
        ]
        mock_csv.DictReader.return_value = csv_dict_data
        with patch('generate_png_from_pdb.open', mock_open(read_data="")):
            result = get_symbol_paths_ary(csv_filepath='/tmp/fakefile.csv')
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0], ('WASH1P',  # approved_symbol
                                     '/tmp/AF-111-F1-model_v2.pdb.gz',  # compressed_pdb_path
                                     '/tmp/WASH1P.pdb',  # pdb_path
                                     '/tmp/WASH1P.png'))
        self.assertEqual(result[1], ('WASH2P',  # approved_symbol
                                     '/tmp/AF-222-F1-model_v2.pdb.gz',  # compressed_pdb_path
                                     '/tmp/WASH2P.pdb',  # pdb_path
                                     '/tmp/WASH2P.png'))

    @patch('generate_png_from_pdb.os.path')
    @patch('generate_png_from_pdb.extract_pdb')
    @patch('generate_png_from_pdb.os.remove')
    @patch('generate_png_from_pdb.cmd')
    def test_convert_pdb_files_to_pngs(self, mock_pymol_cmd, mock_os_remove, mock_extract_pdb, mock_os_path):
        mock_os_path.exists.side_effect = [True, False]  # pdb file exists, but png doesn't
        symbol_paths_ary = [
            ('WASH1P',  # approved_symbol
             '/tmp/AF-111-F1-model_v2.pdb.gz',  # compressed_pdb_path
             '/tmp/WASH1P.pdb',  # pdb_path
             '/tmp/WASH1P.png')

        ]
        with patch('generate_png_from_pdb.gzip.open', mock_open(read_data="")):
            convert_pdb_files_to_pngs(symbol_paths_ary)
        # the pdb gz file should be uncompressed
        mock_extract_pdb.assert_called_with('/tmp/AF-111-F1-model_v2.pdb.gz', '/tmp/WASH1P.pdb')
        # load the pdb file with pymol
        mock_pymol_cmd.load.assert_called_with('/tmp/WASH1P.pdb')
        # save the png with pymol
        mock_pymol_cmd.png.assert_called_with('/tmp/WASH1P.png', dpi=300)
        # no longer removing the uncompressed pdb
        mock_os_remove.assert_not_called()


if __name__ == '__main__':
    unittest.main()
