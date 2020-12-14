# End to End Integration Testing

## Summary

official [pydtest docs](https://docs.pytest.org/en/latest/)

The pytest "fixtures" are kept in the `conftest.py` file.  These are functions that can be re-used by different functions in the main_runner.py. 

In the "TestsWDLs" directory, there are wdls with names like "testcase22.wdl". This name corresponds to the testcase in the "score-card" tests at [google doc](https://docs.google.com/document/d/1nXuPDVZ3dXl0AetyU5Imdbi0Gvc5sUhAR0OfYxss2uI/edit#heading=h.rmy1jmsa0m7n) and [google sheet](https://docs.google.com/spreadsheets/d/1eBWvk4FSPpbFclTuzu0o77aPAxcZ78C_mVKCnHoMMAo/edit#gid=1883830451)

## Run Tests

```
. ~/jaws-dev/bin/activate
pytest main_runner.py
```


