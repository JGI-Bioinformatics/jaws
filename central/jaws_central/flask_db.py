from flask_sqlalchemy import SQLAlchemy
from jaws_central.database import metadata

db = SQLAlchemy(metadata=metadata)
