from streamlit.testing.v1 import AppTest


def test_app_runs_an_analysis(root):
    at = AppTest.from_file(str(root / "oligotox_app.py"), default_timeout=120)
    at.run()
    assert not at.exception
    at.button[0].click().run()
    assert not at.exception
    assert not [e.value for e in at.error]
