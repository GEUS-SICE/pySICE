
import os
from pysice.sice_io import sice_io
from pysice.sice_lib import process

# this test assumes an OLCI file is available in this directory # can be downloaded from scihub.copernicus.eu, creodias or other dias
datadir = "./S3A_OL_1_EFR____20190101T165746_20190101T165817_20190103T041529_0031_039_368_4680_LN1_O_NT_002.SEN3"
datadir = os.path.join(os.path.dirname(__file__), datadir)
assert os.path.exists(os.path.join(datadir, 'Oa01_radiance.nc'))


def test_sice_io():

    olci_data = sice_io(datadir)

    assert hasattr(olci_data, "open") and callable(olci_data.open)

    olci_data.open()

    assert hasattr(olci_data, "olci_scene")

    assert hasattr(olci_data.olci_scene, "toa") and len(olci_data.olci_scene.toa) == 21


def test_process():

    olci_data = sice_io(datadir)
    olci_data.open()

    snow = process(olci_data.olci_scene)

    assert "isnow" in snow
    assert "diameter" in snow
