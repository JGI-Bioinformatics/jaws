from flask_sqlalchemy import SQLAlchemy
from database import metadata

db = SQLAlchemy(metadata=metadata)
