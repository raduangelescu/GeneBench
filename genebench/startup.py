
import wget
import py7zr
import os
from genebench.utils import Utils
from genebench.geoimporter import GEOImporterConfig


def setup_default_data(config_filename):
    data_section = Utils.get_config(config_filename, 'GEOImporter')
    config = GEOImporterConfig(**data_section)
    default_data_url = 'https://github.com/raduangelescu/GeneBench/raw/main/genebench-data.7z'
    download_file = os.path.join(config.data_path, "genebench-data.7z") 
    wget.download(default_data_url, download_file)
    with py7zr.SevenZipFile(download_file, mode='r') as z:
        z.extractall(config.data_path)