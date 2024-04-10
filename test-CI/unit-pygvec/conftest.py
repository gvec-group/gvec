import pytest


def pytest_configure(config):
    """
    Add custom markers to pytest, as if using the `pytest.ini` file.

    Args:
        config (Config): The pytest configuration object.
    """
    # register additional markers
    for marker in [
        "unit: mark test as a unittest - a test that is isolated and tests a specific functionality",
        "pygvec: mark test that are related to the python bindings",
    ]:
        config.addinivalue_line("markers", marker)


def pytest_collection_modifyitems(items):
    """
    Modify collected tests

    Add global markers to each test.

    Args:
        items (List[Item]): The collected pytest testitem objects.
    """
    for item in items:
        if "unit-pygvec" in item.nodeid:
            item.add_marker(pytest.mark.unit)
            item.add_marker(pytest.mark.pygvec)
