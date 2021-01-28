## Scorecard Tests
The following code is for tests that represent [scorecard tests](https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1435741986).

If you are interested in scheduling jobs through the gitlab schedular, see [gitlab schedular](https://docs.google.com/document/d/1Xd5vF31qNfrbgeFhfu3paXWEBq-2dkpqfHjFCpVAq4o/edit#).

official [pydtest docs](https://docs.pytest.org/en/latest/)

The pytest "fixtures" are kept in the `conftest.py` file.  These are functions that can be re-used by different functions in the pytest scripts

## Create the venv 
To create the venv environment you can do

```
# make sure you have python3
module load python
python --version

# create a venv
python -m venv pytest_venv
source pytest_venvbin/activate
pip install -r requirements.txt
```

This is how I created the requirements file
`pip freeze > requirements.txt`

## Run a test
pytest --verbose --env prod test_jaws_cmds.py



