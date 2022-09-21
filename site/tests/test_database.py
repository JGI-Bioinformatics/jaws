# import sqlalchemy
# import jaws_site
#
#
# def test_create_engine_exception(monkeypatch, configuration):
#     def mock_create_engine():
#         raise Exception
#
#     monkeypatch.setattr(sqlalchemy, "create_engine", mock_create_engine)
#
#     from jaws_site.database import Base
#
#     assert Base is not None
