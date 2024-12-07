from seq import validate_input, read_file, analyze_seq
import pytest 

def test_validate_input():
    # valid input
    assert validate_input(["seq.py", "a_seq.txt"]) == True

    # no file path were given
    with pytest.raises(SystemExit) as excinfo:
        validate_input(["seq.py"])

    assert excinfo.value.code == "Usage: python seq.py filepath1.txt filepath2.txt"

    # some files aren't .txt files
    with pytest.raises(SystemExit) as excinfo:
        validate_input(["seq.py", "invalid"])

    assert excinfo.value.code == "Usage: 'invalid' is not a .txt file"
    
    with pytest.raises(SystemExit) as excinfo:
        validate_input(["seq.py", "a_seq.txt", "invalid"])

    assert excinfo.value.code == "Usage: 'invalid' is not a .txt file"

    # some files don't exist
    with pytest.raises(SystemExit) as excinfo:
        validate_input(["seq.py", "invalid.txt"])

    assert excinfo.value.code == "'invalid.txt' doesn't exists"

    with pytest.raises(SystemExit) as excinfo:
        validate_input(["seq.py","a_seq.txt", "invalid.txt"])

    assert excinfo.value.code == "'invalid.txt' doesn't exists"

def test_read_file():
    assert read_file("a_seq.txt") == "ACCTGXXCXXGTTACTGGGCXTTGTXX"
    assert read_file("b_seq.txt") == "ACCGGGTTTT"


def test_analyze_seq():
    assert analyze_seq("abcdT") == {"A": 1, "T": 1, "C": 1, "G": 0, "Unknown": 2, "Total": 5}
    assert analyze_seq("aAAtTcCgGh") == {"A": 3, "T": 2, "C": 2, "G": 2, "Unknown": 1, "Total": 10}
