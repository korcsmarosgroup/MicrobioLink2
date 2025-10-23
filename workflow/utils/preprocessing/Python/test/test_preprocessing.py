import pytest
from preprocessing import CheckingConfiguration
from pathlib import Path
import os


@pytest.fixture
def config_folder():
    return Path(__file__).resolve().parent



def test_CheckingLogFile(config_folder):

    test_log_file = f"{config_folder}/test_preprocessing.log"

    if os.path.isfile(test_log_file):
        os.remove(test_log_file)

    assert not os.path.exists(test_log_file)


def test_CheckingConfiguration(tmp_path, capsys):

    config_folder = tmp_path

    with pytest.raises(SystemExit) as error:
        CheckingConfiguration(str(config_folder))

    assert error.value.code == 1

    error_message = capsys.readouterr()
    assert f"ERROR MESSAGE: The configuration file does not exists:" in error_message.err
