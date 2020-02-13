import os
import connexion
from flask_sqlalchemy import SQLAlchemy

# create Connexion app
basedir = os.path.abspath(os.path.dirname(__file__))
connex = connexion.App(__name__, specification_dir=basedir)

# get underlying Flask instance
app = connex.app

# configure Flask
app.config["SQLALCHEMY_ECHO"] = False
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
app.config.from_envvar('JAWS_FLASK_CONFIG')

# init ORM
db = SQLAlchemy(app)
